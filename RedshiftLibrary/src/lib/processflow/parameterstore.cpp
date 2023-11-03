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
#include "RedshiftLibrary/processflow/parameterstore.h"
#include <cfloat>
namespace bpt = boost::property_tree;

namespace NSEpic {
CParameterStore::CParameterStore(const TScopeStack &stack)
    : CScopeStore(stack) {}

void CParameterStore::Set(const std::string &name, const TFloat64List &v) {
  boost::optional<bpt::ptree &> property =
      m_PropertyTree.get_child_optional(name);

  bpt::ptree array;

  for (Int32 i = 0; i < v.size(); i++) {
    bpt::ptree item;
    item.put("", v[i]);

    array.push_back(std::make_pair("", item));
  }

  m_PropertyTree.put_child(name, array);
}

void CParameterStore::Set(const std::string &name, const TStringList &v) {
  boost::optional<bpt::ptree &> property =
      m_PropertyTree.get_child_optional(name);

  bpt::ptree array;

  for (Int32 i = 0; i < v.size(); i++) {
    bpt::ptree item;
    item.put("", v[i]);

    array.push_back(std::make_pair("", item));
  }

  m_PropertyTree.put_child(name, array);
}

void CParameterStore::Set(const std::string &name, const TFloat64Range &v) {
  TFloat64List list(2);

  list[0] = v.GetBegin();
  list[1] = v.GetEnd();

  Set(name, list);
}

void CParameterStore::Set(const std::string &name, const TInt64List &v) {
  boost::optional<bpt::ptree &> property =
      m_PropertyTree.get_child_optional(name);

  bpt::ptree array;

  for (Int32 i = 0; i < v.size(); i++) {
    bpt::ptree item;
    item.put("", v[i]);

    array.push_back(std::make_pair("", item));
  }

  m_PropertyTree.put_child(name, array);
}

void CParameterStore::Set(const std::string &name, const TBoolList &v) {
  boost::optional<bpt::ptree &> property =
      m_PropertyTree.get_child_optional(name);

  bpt::ptree array;

  for (Int32 i = 0; i < v.size(); i++) {
    bpt::ptree item;
    item.put("", v[i]);

    array.push_back(std::make_pair("", item));
  }

  m_PropertyTree.put_child(name, array);
}

void CParameterStore::Set(const std::string &name, Float64 v) {
  boost::optional<Float64> property =
      m_PropertyTree.get_optional<Float64>(name);

  m_PropertyTree.put(name, v);
}

void CParameterStore::Set(const std::string &name, Int64 v) {
  boost::optional<Int64> property = m_PropertyTree.get_optional<Int64>(name);

  m_PropertyTree.put(name, v);
}

void CParameterStore::Set(const std::string &name, bool v) {
  boost::optional<bool> property = m_PropertyTree.get_optional<bool>(name);

  m_PropertyTree.put(name, v);
}

void CParameterStore::Set(const std::string &name, const std::string &v) {
  boost::optional<std::string> property =
      m_PropertyTree.get_optional<std::string>(name);

  m_PropertyTree.put(name, v);
}

void CParameterStore::Save(const std::string &path) const {
  bpt::json_parser::write_json(path, m_PropertyTree);
}

void CParameterStore::FromString(const std::string &json) {
  std::istringstream jsonstream(json);
  bpt::json_parser::read_json(jsonstream, m_PropertyTree);
}

bool CParameterStore::HasTplIsmExtinction(const std::string &spectrumModel) const {
  bool extinction = false;
  const std::string methodScope = spectrumModel + ".method";
  if (!Has<std::string>(methodScope))
    return false;
  const auto &method = Get<std::string>(methodScope);
  if (method == "templateFittingSolver" || method == "tplCombinationSolver") {
    const std::string scopeStr = spectrumModel + "." + method + ".ismFit";
    if (Has<bool>(scopeStr))
      extinction = Get<bool>(scopeStr);
  } else if (method == "lineModelSolver") {
    // two parameters play here: continuumfit.ismfit and tplRatioIsmFit
    //
    const std::string scopeStr =
        spectrumModel + "." + method + ".lineModel.continuumFit.ismFit";
    if (Has<bool>(scopeStr))
      extinction = Get<bool>(scopeStr);

    const std::string scopeStr_tplratio =
        spectrumModel + "." + method + ".lineModel.tplRatioIsmFit";
    if (Has<bool>(scopeStr_tplratio))
      extinction |= Get<bool>(scopeStr_tplratio);
  }
  return extinction;
}

bool CParameterStore::HasTplIgmExtinction(const std::string &spectrumModel) const {
  bool extinction = false;
  const std::string methodScope = spectrumModel + ".method";
  if (!Has<std::string>(methodScope))
    return false;
  const auto &method = Get<std::string>(methodScope);
  std::string scopeStr = spectrumModel + "." + method;
  if (method == "templateFittingSolver" || method == "tplCombinationSolver")
    scopeStr += ".igmFit";
  else if (method == "lineModelSolver")
    scopeStr += ".lineModel.continuumFit.igmFit";

  if (Has<bool>(scopeStr))
    extinction = Get<bool>(scopeStr);

  return extinction;
}

bool CParameterStore::HasFFTProcessing(const std::string &spectrumModel) const {
  bool fft_processing = false;
  const std::string methodScope = spectrumModel + ".method";
  if (!Has<std::string>(methodScope))
    return false;
  const auto &method = Get<std::string>(methodScope);
  if (method == "templateFittingSolver") {
    const std::string scopeStr =
        spectrumModel + ".templateFittingSolver.fftProcessing";
    if (Has<bool>(scopeStr))
      fft_processing = Get<bool>(scopeStr);
  }
  if (method == "lineModelSolver") {
    const std::string scopeStr =
        spectrumModel + ".lineModelSolver.lineModel.continuumFit.fftProcessing";
    if (Has<bool>(scopeStr))
      fft_processing = Get<bool>(scopeStr);
  }

  return fft_processing;
}

bool CParameterStore::HasToOrthogonalizeTemplates(
    const std::string &spectrumModel) const {

  const std::string methodScope = spectrumModel + ".method";
  bool orthogonalize = Get<std::string>(methodScope) == "lineModelSolver";
  if (orthogonalize) {
    std::string continuumComponent = Get<std::string>(
        spectrumModel + ".lineModelSolver.lineModel.continuumComponent");
    orthogonalize &=
        (continuumComponent == "tplFit" || continuumComponent == "tplFitAuto");
  }
  return orthogonalize;
}

bool CParameterStore::EnableTemplateOrthogonalization(
    const std::string &spectrumModel) const {
  bool enableOrtho = HasToOrthogonalizeTemplates(spectrumModel);
  if (enableOrtho) {
    enableOrtho &=
        !Get<bool>(spectrumModel +
                   ".lineModelSolver.lineModel.continuumFit.ignoreLineSupport");
  }
  return enableOrtho;
}

bool CParameterStore::hasToLogRebin(
    const TStringList &categories,
    std::map<std::string, bool> &fft_processing) const {
  bool usefftSpectrum = false;
  for (std::string cat : categories) {
    fft_processing[cat] = HasFFTProcessing(cat);
    if (fft_processing[cat])
      usefftSpectrum = true;
  }
  return usefftSpectrum;
}

Float64 CParameterStore::getMinZStepForFFTProcessing(
    std::map<std::string, bool> fftprocessing) const {
  Float64 logGridStep = DBL_MAX;
  for (const auto &p : fftprocessing) {
    if (!p.second)
      continue;
    Float64 redshift_step = Get<Float64>(p.first + ".redshiftStep");
    logGridStep = std::min(redshift_step, logGridStep);
  }
  return logGridStep;
}
} // namespace NSEpic
