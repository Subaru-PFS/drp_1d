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
#include <algorithm> // std::sort
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/line/lineprofile.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/log/log.h"

using namespace NSEpic;
using namespace std;
using namespace boost;
using namespace boost::filesystem;

template <typename TLine>
typename CLineCatalogBase<TLine>::TLineMap
CLineCatalogBase<TLine>::GetFilteredList(CLine::EType typeFilter,
                                         CLine::EForce forceFilter) const {

  TLineMap filteredList;

  for (const auto &[id, line] : m_List) {
    if ((typeFilter == CLine::EType::nType_All ||
         typeFilter == line.GetType()) &&
        (forceFilter == CLine::EForce::nForce_All ||
         forceFilter == line.GetForce())) {
      filteredList[id] = line;
    }
  }

  return filteredList;
}

template <typename TLine>
typename CLineCatalogBase<TLine>::TLineMap
CLineCatalogBase<TLine>::GetFilteredList(const std::string &typeFilter,
                                         const std::string &forceFilter) const {
  const auto &etypeFilter = CLine::string2Type(typeFilter);
  const auto &eforceFilter = CLine::string2Force(forceFilter);

  return GetFilteredList(etypeFilter, eforceFilter);
}

template <typename TLine>
std::map<std::string, typename CLineCatalogBase<TLine>::TLineVector>
CLineCatalogBase<TLine>::ConvertToGroupList(
    const CLineCatalogBase<TLine>::TLineMap &filteredList) {

  std::map<std::string, TLineVector> fullList;

  for (const auto &[_, line] : filteredList) {
    auto group_name = line.GetGroupName();
    if (group_name == undefStr)
      // non grouped lines are added in dedicated maps (one element)
      group_name = "single_" + line.GetStrID();
    fullList[group_name].push_back(line);
  }

  return fullList;
}

void CLineCatalog::AddLineFromParams(
    const std::string &name, Float64 position, const std::string &type,
    const std::string &force, const std::string &profileName,
    const TAsymParams &asymParams, const std::string &groupName,
    Float64 nominalAmplitude, const std::string &velocityGroup,
    Float64 velocityOffset, bool enableVelocityFit, Int32 id,
    const std::string &str_id,
    const std::shared_ptr<CSpectrumFluxCorrectionMeiksin> &igmcorrection) {

  const auto &etype = CLine::string2Type(type);
  const auto &eforce = CLine::string2Force(force);

  TAsymParams _asymParams = {1., 4.5, 0.};
  TAsymParams _asymFitParams = {2., 2., 0.};
  std::unique_ptr<CLineProfile> profile;

  if (profileName.find("ASYMFIXED") != std::string::npos)
    profile = std::unique_ptr<CLineProfileASYM>(
        new CLineProfileASYM(m_nSigmaSupport, asymParams, "mean"));
  else if (profileName == "SYM")
    profile =
        std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM(m_nSigmaSupport));
  else if (profileName == "LOR")
    profile =
        std::unique_ptr<CLineProfileLOR>(new CLineProfileLOR(m_nSigmaSupport));
  else if (profileName == "ASYM")
    profile = std::unique_ptr<CLineProfileASYM>(
        new CLineProfileASYM(m_nSigmaSupport, _asymParams, "none"));
  else if (profileName == "ASYMFIT")
    profile = std::unique_ptr<CLineProfileASYMFIT>(
        new CLineProfileASYMFIT(m_nSigmaSupport, _asymFitParams, "mean"));
  else if (profileName == "SYMIGM")
    profile = std::unique_ptr<CLineProfileSYMIGM>(
        new CLineProfileSYMIGM(igmcorrection, m_nSigmaSupport));
  else {
    THROWG(ErrorCode::INTERNAL_ERROR, Formatter()
                                          << "Profile name " << profileName
                                          << " is no recognized.");
  }

  Add(CLine(name, position, etype, std::move(profile), eforce, velocityOffset,
            enableVelocityFit, groupName, nominalAmplitude, velocityGroup, id,
            str_id));
}

template <typename TLine>
void CLineCatalogBase<TLine>::setLineAmplitude(Int32 id,
                                               Float64 nominalAmplitude) {

  const auto &search = m_List.find(id);
  if (search != m_List.end())
    search->second.setNominalAmplitude(nominalAmplitude);
  else
    THROWG(ErrorCode::INTERNAL_ERROR, Formatter()
                                          << " Line with id " << id
                                          << " does not exist in catalog");
}

/**
 * @brief only called for tplRatio catalogs
 * hereinafter we consider that only one ASYMFIT profile exists in linecatalog
 * and that it corresponds to Lya
 * @param profile
 * @param params
 */
template <typename TLine>
void CLineCatalogBase<TLine>::setAsymProfileAndParams(
    const std::string &profile, TAsymParams params) {

  for (auto &it : m_List) {
    auto &line = it.second;
    if (line.GetProfile()->isAsymFit() || line.GetProfile()->isAsymFixed())
      line.setProfileAndParams(profile, params, m_nSigmaSupport);
  }
}

template <typename TLine>
void CLineCatalogBase<TLine>::convertLineProfiles2SYMIGM(
    const std::shared_ptr<CSpectrumFluxCorrectionMeiksin> &igmcorrection) {

  for (auto &it : m_List) {
    auto &line = it.second;
    if (!line.IsEmission())
      continue;
    if (line.GetName() == linetags::lya_em ||
        line.GetPosition() < RESTLAMBDA_LYA)
      line.setProfileAndParams("SYMIGM", TAsymParams(), m_nSigmaSupport,
                               igmcorrection);
  }
}

void CLineDetectedCatalog::Sort() {

  // order the result by value instead of key
  typename CLineDetectedCatalog::TLineVector buff;
  buff.reserve(m_List.size());
  for (auto &[id, line] : m_List)
    buff.push_back(std::move(line));

  sort(buff.begin(), buff.end());

  m_List.clear();
  for (Int32 i = 0, end = buff.size(); i < end; ++i)
    m_List[i] = std::move(buff[i]);
}
