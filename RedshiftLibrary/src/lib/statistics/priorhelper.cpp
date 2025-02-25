// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
//
// https://www.lam.fr/
//
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
//
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use,
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info".
//
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability.
//
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or
// data to be ensured and,  more generally, to use and operate it in the
// same conditions as regards security.
//
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include <fstream>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/replace.hpp>

#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/statistics/priorhelper.h"

using namespace NSEpic;
using namespace std;
using namespace boost;

void CPriorHelper::Init(std::string priorDirPath, Int32 type) {
  bfs::path rootFolder(priorDirPath.c_str());
  if (!bfs::exists(rootFolder)) {
    if (!rootFolder.string().empty()) {
      Flag.warning(WarningCode::INVALID_FOLDER_PATH,
                   Formatter() << "    CPriorHelper::" << __func__
                               << ": rootFolder path does not exist: "
                               << rootFolder.string());
      Log.LogDetail("    CPriorHelper: priors won't be used");
    }
    mInitFailed = true;
    return;
  }
  m_type = type;
  std::string part;
  if (m_type == 0)
    part = "continuum";
  else if (m_type == 1)
    part = "lines";
  TStringList EZTfilesPathList;
  bfs::directory_iterator end_itr;
  std::string ezt_path = "prior_" + part + "_hist_Ebmvr_Z";
  for (bfs::directory_iterator itr(rootFolder / ezt_path.c_str());
       itr != end_itr; ++itr)
    if (!is_directory(itr->status()))
      EZTfilesPathList.push_back(itr->path().c_str());

  TStringList AGaussMeanfilesPathList =
      loadfilesPathList(EZTfilesPathList, rootFolder, part, "mean");
  TStringList AGaussSigmafilesPathList =
      loadfilesPathList(EZTfilesPathList, rootFolder, part, "sigma");

  // allocate the data buffer
  setSize(EZTfilesPathList.size());

  // set the template names
  for (Int32 k = 0; k < ssize(EZTfilesPathList); k++) {
    bfs::path fPath = EZTfilesPathList[k];
    std::string fNameStr = fPath.filename().c_str();
    setTNameData(k, fNameStr);
  }

  // read the EZT data from files
  for (Int32 k = 0; k < ssize(EZTfilesPathList); k++) {
    bfs::path fPath = EZTfilesPathList[k];
    std::string fPathStr = (fPath).string();

    std::vector<TFloat64List> read_buffer;
    loadFileEZ(fPathStr.c_str(), read_buffer);
    if (mInitFailed)
      return;
    setEZTData(k, read_buffer);
    mInitFailed = false;
  }

  // read the AGaussMean data from files
  loadFromFileList("mean", AGaussMeanfilesPathList);
  // read the AGaussSigma data from files
  loadFromFileList("sigma", AGaussSigmafilesPathList);

  std::string z_dirpath = "prior_" + part + "_hist_Z";
  std::string z_filepath = "prior_" + part + "_hist_Z.txt";
  bfs::path pz_fpath = rootFolder / z_dirpath.c_str() / z_filepath.c_str();
  if (!bfs::exists(pz_fpath))
    THROWG(ErrorCode::INVALID_FILEPATH,
           Formatter() << "Pz path does not exist: " << pz_fpath.string());

  TFloat64List read_buffer;
  loadFileZ(pz_fpath.string().c_str(), read_buffer);
  if (mInitFailed)
    return;
  setPzData(read_buffer);
  mInitFailed = false;
}

void CPriorHelper::SetBetaA(Float64 beta) { m_betaA = beta; }
void CPriorHelper::SetBetaTE(Float64 beta) { m_betaTE = beta; }
void CPriorHelper::SetBetaZ(Float64 beta) { m_betaZ = beta; }

void CPriorHelper::setSize(Int32 size) {
  m_data.clear();
  SPriorTZE _tze;
  std::vector<SPriorTZE> _elist(m_nEbv, _tze);
  TPriorZEList _zelist(m_nZ, _elist);
  m_data.resize(size, _zelist);
  m_data_pz.clear();
  m_data_pz.resize(m_nZ, 0.0);

  m_tplnames.clear();
  m_tplnames.resize(size, "undefined");
  return;
}

const TStringList CPriorHelper::loadfilesPathList(
    const TStringList &EZTfilesPathList, const bfs::path &rootFolder,
    const std::string &part, const std::string &sigma_mean) const {
  if (sigma_mean != "mean" && sigma_mean != "sigma")
    THROWG(ErrorCode::INTERNAL_ERROR, "Invalid argument");
  if (part != "continuum" && part != "lines")
    THROWG(ErrorCode::INTERNAL_ERROR, "Invalid argument");

  std::string ezt_tag = "prior_" + part + "_hist_Ebmvc_Z_Tplc_";
  std::string a_tag = "prior_" + part + "_gauss" + sigma_mean + "_Ac_Z_Tplc_";
  std::string a_dirpath = "prior_" + part + "_gauss" + sigma_mean + "_Ac_Z";

  TStringList AGaussfilesPathList;
  AGaussfilesPathList.reserve(EZTfilesPathList.size());
  for (bfs::path fPath : EZTfilesPathList) {
    std::string fNameStr = fPath.filename().c_str();

    boost::replace_all(fNameStr, ezt_tag.c_str(), a_tag.c_str());

    bfs::path agaussfpath =
        rootFolder / a_dirpath.c_str() / bfs::path(fNameStr);
    if (!bfs::exists(agaussfpath)) {
      THROWG(ErrorCode::INVALID_FILEPATH,
             Formatter() << "path does not exist: " << agaussfpath.string());
    }
    AGaussfilesPathList.push_back(agaussfpath.string());
  }
  return AGaussfilesPathList;
}

void CPriorHelper::loadFromFileList(const std::string &sigma_mean,
                                    const TStringList &filePathList) {
  // read the AGaussSigma data from files
  for (Int32 k = 0; k < ssize(filePathList); k++) {
    bfs::path fPath = filePathList[k];
    std::string fPathStr = (fPath).string();

    std::vector<TFloat64List> read_buffer;
    loadFileEZ(fPathStr.c_str(), read_buffer);
    if (mInitFailed)
      return;
    if (sigma_mean == "sigma")
      setAGausssigmaData(k, read_buffer);
    else if (sigma_mean == "mean")
      setAGaussmeanData(k, read_buffer);

    mInitFailed = false;
  }
  return;
}

void CPriorHelper::setTNameData(Int32 k, std::string tname) {
  if (k >= ssize(m_tplnames))
    THROWG(ErrorCode::INTERNAL_ERROR, Formatter()
                                          << "Out-of-bound index k=" << k);

  m_tplnames[k] = tname;
}

void CPriorHelper::setEZTData(Int32 k,
                              const std::vector<TFloat64List> &ezt_data) {
  if (k >= ssize(m_data))
    THROWG(ErrorCode::INTERNAL_ERROR, Formatter()
                                          << "Out-of-bound index k=" << k);

  for (Int32 kz = 0; kz < m_nZ; kz++)
    for (Int32 ke = 0; ke < m_nEbv; ke++)
      m_data[k][kz][ke].priorTZE = ezt_data[kz][ke];
}

void CPriorHelper::setAGaussmeanData(
    Int32 k, const std::vector<TFloat64List> &agaussmean_data) {
  if (k >= ssize(m_data))
    THROWG(ErrorCode::INTERNAL_ERROR, Formatter()
                                          << "Out-of-bound index k=" << k);

  for (Int32 kz = 0; kz < m_nZ; kz++)
    for (Int32 ke = 0; ke < m_nEbv; ke++)
      m_data[k][kz][ke].A_mean = agaussmean_data[kz][ke];
}

void CPriorHelper::setAGausssigmaData(
    Int32 k, const std::vector<TFloat64List> &agausssigma_data) {
  if (k >= ssize(m_data))
    THROWG(ErrorCode::INTERNAL_ERROR, Formatter()
                                          << "Out-of-bound index k=" << k);

  for (Int32 kz = 0; kz < m_nZ; kz++)
    for (Int32 ke = 0; ke < m_nEbv; ke++)
      m_data[k][kz][ke].A_sigma = agausssigma_data[kz][ke];
}

void CPriorHelper::setPzData(const TFloat64List &z_data) {
  if (z_data.size() != m_data_pz.size())
    THROWG(ErrorCode::INTERNAL_ERROR, " Data and prior sizes do not match");

  std::copy(z_data.cbegin(), z_data.cbegin() + m_nZ, m_data_pz.begin());
}

void CPriorHelper::loadFileEZ(const char *filePath,
                              std::vector<TFloat64List> &data) {
  Log.LogDetail(Formatter()
                << "CPriorHelper: start load prior file: " << filePath);
  std::ifstream file;
  file.open(filePath, std::ifstream::in);
  bool fileOpenFailed = file.rdstate() & std::ios_base::failbit;
  if (fileOpenFailed)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Unable to load the prior data from file: "
                       << filePath);

  Int32 nlinesRead = 0;
  std::string line;
  // Read file line by line
  while (getline(file, line)) {
    if (boost::starts_with(line, "#"))
      continue;

    TFloat64List lineVals(m_nEbv);
    std::istringstream iss(line);
    for (auto &x : lineVals) {
      iss >> x;
      if (!isValidFloat(x) || x <= 0.)
        x = m_priorminval;
    }

    nlinesRead++;
    data.push_back(lineVals);
  }
  file.close();
}

void CPriorHelper::loadFileZ(const char *filePath, TFloat64List &data) {
  Log.LogDetail(Formatter()
                << "    CPriorHelper: start load prior file: " << filePath);
  std::ifstream file;
  file.open(filePath, std::ifstream::in);
  bool fileOpenFailed = file.rdstate() & std::ios_base::failbit;
  if (fileOpenFailed)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Unable to load the prior data from file: "
                       << filePath);
  Int32 nlinesRead = 0;
  std::string line;
  // Read file line by line
  while (getline(file, line)) {
    if (!boost::starts_with(line, "#")) {
      std::istringstream iss(line);
      Float64 x;
      iss >> x;
      if (!isValidFloat(x) || x < 0.) {
        x = m_priorminval;
      }

      nlinesRead++;
      data.push_back(x);
    }
  }
  file.close();

  if (nlinesRead != m_nZ)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "read n=" << nlinesRead << " lines");
}

/**
 * @brief CPriorHelper::GetTplPriorData
 * @param tplname
 * @param redshifts
 * @param zePriorData
 * @param outsideZRangeExtensionMode: 0=extend 0 value for z<m_z0, and n-1
 * value for z>m_z0+m_dZ*m_nZ, 1=return error if z outside prior range
 * @return
 */
void CPriorHelper::GetTplPriorData(const std::string &tplname,
                                   const TRedshiftList &redshifts,
                                   TPriorZEList &zePriorData,
                                   Int32 outsideZRangeExtensionMode) const {
  if (m_betaA <= 0.0 && m_betaTE <= 0.0 && m_betaZ <= 0.0)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "beta coeffs are all zero (betaA=" << m_betaA
                       << ", betaTE=" << m_betaTE << ", betaZ=" << m_betaZ
                       << ")");

  if (mInitFailed) {
    Log.LogDetail(
        "CPriorHelper::GetTplPriorData: failed, unable to provide priors. "
        "Priors won't be used.");
    zePriorData.clear();
    return;
  }

  // find idx for tplname
  Int32 idx = getTemplateIndex(tplname);

  zePriorData.clear();
  zePriorData.reserve(redshifts.size());
  for (Float64 z : redshifts) {
    Int32 idz = getRedshiftIndex(z, outsideZRangeExtensionMode);
    TPriorEList dataz = m_data[idx][idz];
    fillPriorDataPerZ(dataz, idz);
    zePriorData.push_back(std::move(dataz));
  }
}

void CPriorHelper::fillPriorDataPerZ(TPriorEList &dataz, Int32 idz) const {

  for (auto &dataz_i : dataz) {
    if (dataz_i.priorTZE <= 0.0)
      THROWG(ErrorCode::INTERNAL_ERROR, "P_TZE cannot be null");

    Float64 logPA = dataz_i.A_sigma > 0.0
                        ? -0.5 * log(2 * M_PI) - log(dataz_i.A_sigma)
                        : log(1. / m_deltaA);

    Float64 logPTE = log(dataz_i.priorTZE);

    if (!isValidFloat(logPTE))
      THROWG(ErrorCode::INTERNAL_ERROR, "logP_TZE is NAN or inf, or invalid");

    Float64 logPZ = m_data_pz[idz] > 0.0 ? log(m_data_pz[idz] / m_dz) : 0.;

    dataz_i.logprior_precompA = logPA;
    dataz_i.logprior_precompTE = logPTE;
    dataz_i.logprior_precompZ = logPZ;
    dataz_i.betaA = m_betaA;
    dataz_i.betaTE = m_betaTE;
    dataz_i.betaZ = m_betaZ;
  }
}

CPriorHelper::SPriorTZE
CPriorHelper::GetTZEPriorData(const std::string &tplname, Int32 EBVIndexfilter,
                              Float64 redshift,
                              Int32 outsideZRangeExtensionMode) const {
  if (EBVIndexfilter < 0 || EBVIndexfilter > m_nEbv - 1)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Bad EBV index requested =" << EBVIndexfilter
                       << " nEBV=" << m_nEbv);

  TFloat64List redshifts(1, redshift);
  TPriorZEList zePriorData;
  GetTplPriorData(tplname, redshifts, zePriorData, outsideZRangeExtensionMode);

  return zePriorData[0][EBVIndexfilter];
}

Int32 CPriorHelper::getTemplateIndex(const std::string &tplname) const {

  auto itr = std::find(m_tplnames.begin(), m_tplnames.end(), tplname);
  if (itr == m_tplnames.end())
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "unable to match this tplname in prior names list : "
                       << tplname);

  return itr - m_tplnames.begin();
}

Int32 CPriorHelper::getRedshiftIndex(Float64 redshift,
                                     Int32 outsideZRangeExtensionMode) const {
  Int32 idz = Int32((redshift - m_z0) / m_dz);
  if (outsideZRangeExtensionMode == 1) {
    if (idz < 0 || idz >= m_nZ)
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "unable to match this redshift in prior list:"
                         << redshift);
  }
  if (outsideZRangeExtensionMode == 0) {
    idz = std::max(idz, 0);
    idz = std::min(idz, m_nZ);
  }
  return idz;
}

bool CPriorHelper::isValidFloat(Float64 val) const {
  bool invalid = std::isnan(val) || val != val || std::isinf(val);
  return !invalid;
}
