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
#include "RedshiftLibrary/processflow/resultstore.h"

#include "RedshiftLibrary/linemodel/linemodelextremaresult.h"
#include "RedshiftLibrary/method/classificationresult.h"
#include "RedshiftLibrary/method/linemodelsolveresult.h"
#include "RedshiftLibrary/method/reliabilityresult.h"
#include "RedshiftLibrary/operator/extremaresult.h"
#include "RedshiftLibrary/operator/flagResult.h"
#include "RedshiftLibrary/operator/logZPdfResult.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/operator/tplCombinationExtremaResult.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/statistics/pdfcandidateszresult.h"

#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/log/log.h"

using namespace NSEpic;

COperatorResultStore::COperatorResultStore(
    const std::shared_ptr<const CScopeStack> &scope)
    : CScopeStore(scope) {}

void COperatorResultStore::StoreResult(
    TResultsMap &map, const std::string &path, const std::string &name,
    std::shared_ptr<const COperatorResult> result) {
  std::string scopedName;
  if (!path.empty()) {
    scopedName = path;
    scopedName.append(".");
  }
  scopedName.append(name);

  TResultsMap::iterator it = map.find(scopedName);
  if (it != map.end()) {
    THROWG(ErrorCode::INTERNAL_ERROR,
           "Can not store results: result already exists");
  }
  map[scopedName] = result;
}

void COperatorResultStore::StorePerTemplateResult(
    const std::shared_ptr<const CTemplate> &t, const std::string &path,
    const std::string &name, std::shared_ptr<const COperatorResult> result) {
  TPerTemplateResultsMap::iterator it = m_PerTemplateResults.find(t->GetName());
  if (it == m_PerTemplateResults.end()) {
    m_PerTemplateResults[t->GetName()] = TResultsMap();
  }

  StoreResult(m_PerTemplateResults[t->GetName()], path, name, result);
}

void COperatorResultStore::StoreGlobalResult(
    const std::string &path, const std::string &name,
    std::shared_ptr<const COperatorResult> result) {
  StoreResult(m_GlobalResults, path, name, result);
}

std::weak_ptr<const COperatorResult> COperatorResultStore::GetPerTemplateResult(
    const std::shared_ptr<const CTemplate> &t, const std::string &name) const {
  TPerTemplateResultsMap::const_iterator it1 =
      m_PerTemplateResults.find(t->GetName());
  if (it1 != m_PerTemplateResults.end()) {
    const TResultsMap &resultsMap = (*it1).second;
    TResultsMap::const_iterator it2 = resultsMap.find(name);
    if (it2 != resultsMap.end()) {
      return (*it2).second;
    }
  }

  THROWG(ErrorCode::INTERNAL_ERROR,
         Formatter() << "Not found result for template " << name);
  return std::weak_ptr<const COperatorResult>();
}

std::weak_ptr<const COperatorResult>
COperatorResultStore::GetScopedPerTemplateResult(
    const std::shared_ptr<const CTemplate> &t, const std::string &name) const {
  return GetPerTemplateResult(t, GetScopedName(name));
}

TResultsMap COperatorResultStore::GetPerTemplateResult(
    const std::string &resultName) const {
  TResultsMap map;

  for (auto &[tplName, resultsMap] : m_PerTemplateResults) {
    auto foundResult = resultsMap.find(resultName);
    if (foundResult != resultsMap.end())
      map[tplName] = (*foundResult).second;
  }
  return map;
}

TResultsMap COperatorResultStore::GetScopedPerTemplateResult(
    const std::string &name) const {
  return GetPerTemplateResult(GetScopedName(name));
}

/**
 * /brief Returns the best global result, if there is one.
 *
 * The somewhat strange syntax of this method is due to the usage of std::map,
 * which will yield an iterator when accessing a member using .find(). This
 * method will find the global result entry with key "name", and if it exists,
 * will return its .second member (which will be a COperatorResult or derived
 * object). Otherwise, will return a pointer to an empty COperatorResult.
 */

std::weak_ptr<const COperatorResult>
COperatorResultStore::GetGlobalResult(const std::string &name) const {
  TResultsMap::const_iterator it = m_GlobalResults.find(name);
  if (it != m_GlobalResults.end()) {
    return (*it).second;
  } else
    THROWG(ErrorCode::UNKNOWN_ATTRIBUTE,
           Formatter() << "Unknown global result:" << name);
}

std::weak_ptr<const COperatorResult>
COperatorResultStore::GetScopedGlobalResult(const std::string &name) const {
  return GetGlobalResult(GetScopedName(name));
}

std::string COperatorResultStore::buildFullname(
    const std::string &spectrumModel, const std::string &stage,
    const std::string &method, const std::string &name) {

  std::ostringstream oss;
  if (!spectrumModel.empty())
    oss << spectrumModel << ".";
  if (!stage.empty())
    oss << stage << ".";
  if (!method.empty())
    oss << method << ".";
  oss << name;
  return oss.str();
}

std::weak_ptr<const COperatorResult> COperatorResultStore::GetGlobalResult(
    const std::string &spectrumModel, const std::string &stage,
    const std::string &method, const std::string &name) const {
  return GetGlobalResult(buildFullname(spectrumModel, stage, method, name));
}

std::shared_ptr<const CClassificationResult>
COperatorResultStore::GetClassificationResult(const std::string &spectrumModel,
                                              const std::string &stage,
                                              const std::string &method,
                                              const std::string &name) const {
  return std::dynamic_pointer_cast<const CClassificationResult>(
      GetGlobalResult(buildFullname(spectrumModel, stage, method, name))
          .lock());
}

std::shared_ptr<const CReliabilityResult>
COperatorResultStore::GetReliabilityResult(const std::string &spectrumModel,
                                           const std::string &stage,
                                           const std::string &method,
                                           const std::string &name) const {
  return std::dynamic_pointer_cast<const CReliabilityResult>(
      GetGlobalResult(buildFullname(spectrumModel, stage, method, name))
          .lock());
}

std::shared_ptr<const CLogZPdfResult> COperatorResultStore::GetLogZPdfResult(
    const std::string &spectrumModel, const std::string &stage,
    const std::string &method, const std::string &name) const {
  return std::dynamic_pointer_cast<const CLogZPdfResult>(
      GetGlobalResult(buildFullname(spectrumModel, stage, method, name))
          .lock());
}

std::shared_ptr<const CFlagLogResult> COperatorResultStore::GetFlagLogResult(
    const std::string &spectrumModel, const std::string &stage,
    const std::string &method, const std::string &name) const {
  std::ostringstream oss;
  return std::dynamic_pointer_cast<const CFlagLogResult>(
      GetGlobalResult(buildFullname(spectrumModel, stage, method, name))
          .lock());
}

std::shared_ptr<const TLineModelResult>
COperatorResultStore::GetLineModelResult(
    const std::string &spectrumModel, const std::string &stage,
    const std::string &method, const std::string &name,
    const std::string &dataset, const int &rank,
    bool firstpassCorrespondingResult) const

{
  std::shared_ptr<const COperatorResult> cop =
      GetGlobalResult(spectrumModel, stage, method, name)
          .lock()
          ->getCandidate(rank, dataset, firstpassCorrespondingResult);

  std::shared_ptr<const TLineModelResult> tlm =
      std::dynamic_pointer_cast<const TLineModelResult>(cop);
  if (tlm == nullptr && cop != nullptr)
    THROWG(ErrorCode::INTERNAL_ERROR, "tlm is nullptr from GetLineModelResult");
  return tlm;
}

std::shared_ptr<const TTplCombinationResult>
COperatorResultStore::GetTplCombinationResult(const std::string &spectrumModel,
                                              const std::string &stage,
                                              const std::string &method,
                                              const std::string &name,
                                              const std::string &dataset,
                                              const int &rank) const

{
  std::shared_ptr<const COperatorResult> cop =
      GetGlobalResult(spectrumModel, stage, method, name)
          .lock()
          ->getCandidate(rank, dataset);
  std::shared_ptr<const TTplCombinationResult> ttc =
      std::dynamic_pointer_cast<const TTplCombinationResult>(cop);
  return ttc;
}

std::shared_ptr<const TExtremaResult> COperatorResultStore::GetExtremaResult(
    const std::string &spectrumModel, const std::string &stage,
    const std::string &method, const std::string &name,
    const std::string &dataset, const int &rank,
    bool firstpassCorrespondingResult) const

{
  std::shared_ptr<const COperatorResult> cop =
      GetGlobalResult(spectrumModel, stage, method, name)
          .lock()
          ->getCandidate(rank, dataset, firstpassCorrespondingResult);
  std::shared_ptr<const TExtremaResult> tlm =
      std::dynamic_pointer_cast<const TExtremaResult>(cop);
  return tlm;
}

std::shared_ptr<const CLineModelSolution>
COperatorResultStore::GetLineModelSolution(const std::string &spectrumModel,
                                           const std::string &stage,
                                           const std::string &method,
                                           const std::string &name) const {
  return std::dynamic_pointer_cast<const CLineModelSolution>(
      GetGlobalResult(spectrumModel, stage, method, name).lock());
}

std::shared_ptr<const CModelSpectrumResult>
COperatorResultStore::GetModelSpectrumResult(const std::string &spectrumModel,
                                             const std::string &stage,
                                             const std::string &method,
                                             const std::string &name,
                                             const std::string &dataset,
                                             const int &rank) const

{
  return std::dynamic_pointer_cast<const CModelSpectrumResult>(
      GetGlobalResult(spectrumModel, stage, method, name)
          .lock()
          ->getCandidate(rank, dataset));
}

std::shared_ptr<const CModelPhotValueResult>
COperatorResultStore::GetModelPhotValueResult(const std::string &spectrumModel,
                                              const std::string &stage,
                                              const std::string &method,
                                              const std::string &name,
                                              const std::string &dataset,
                                              const int &rank) const

{
  return std::dynamic_pointer_cast<const CModelPhotValueResult>(
      GetGlobalResult(spectrumModel, stage, method, name)
          .lock()
          ->getCandidate(rank, dataset));
}

std::shared_ptr<const CLineModelSolution>
COperatorResultStore::GetLineModelSolution(const std::string &spectrumModel,
                                           const std::string &stage,
                                           const std::string &method,
                                           const std::string &name,
                                           const std::string &dataset,
                                           const int &rank) const

{
  return std::dynamic_pointer_cast<const CLineModelSolution>(
      GetGlobalResult(spectrumModel, stage, method, name)
          .lock()
          ->getCandidate(rank, dataset));
}

std::shared_ptr<const CModelSpectrumResult>
COperatorResultStore::GetModelSpectrumResult(const std::string &spectrumModel,
                                             const std::string &stage,
                                             const std::string &method,
                                             const std::string &name) const

{
  return std::dynamic_pointer_cast<const CModelSpectrumResult>(
      GetGlobalResult(spectrumModel, stage, method, name).lock());
}

const std::string &COperatorResultStore::GetGlobalResultType(
    const std::string &spectrumModel, const std::string &stage,
    const std::string &method, const std::string &name) const {
  return GetGlobalResult(spectrumModel, stage, method, name).lock()->getType();
}

const std::string &COperatorResultStore::GetCandidateResultType(
    const std::string &spectrumModel, const std::string &stage,
    const std::string &method, const std::string &name,
    const std::string &dataset) const {
  return GetGlobalResult(spectrumModel, stage, method, name)
      .lock()
      ->getCandidateDatasetType(dataset);
}

bool COperatorResultStore::HasCandidateDataset(
    const std::string &spectrumModel, const std::string &stage,
    const std::string &method, const std::string &name,
    const std::string &dataset) const {
  if (HasDataset(spectrumModel, stage, method, name)) {
    return GetGlobalResult(spectrumModel, stage, method, name)
        .lock()
        ->HasCandidateDataset(dataset);
  } else
    return false;
}

bool COperatorResultStore::HasDataset(const std::string &spectrumModel,
                                      const std::string &stage,
                                      const std::string &method,
                                      const std::string &name) const {
  TResultsMap::const_iterator it =
      m_GlobalResults.find(buildFullname(spectrumModel, stage, method, name));
  return (it != m_GlobalResults.end());
}

bool COperatorResultStore::hasContextWarningFlag() const {

  TResultsMap::const_iterator it = m_GlobalResults.find("context_warningFlag");
  return (it != m_GlobalResults.end());
}

bool COperatorResultStore::hasInitWarningFlag() const {

  TResultsMap::const_iterator it = m_GlobalResults.find("init_warningFlag");
  return (it != m_GlobalResults.end());
}

bool COperatorResultStore::hasCurrentMethodWarningFlag() const {
  TResultsMap::const_iterator it =
      m_GlobalResults.find(GetScopedNameAt("warningFlag", ScopeType::METHOD));
  return (it != m_GlobalResults.end());
}

int COperatorResultStore::getNbRedshiftCandidates(
    const std::string &spectrumModel, const std::string &stage,
    const std::string &method) const {
  std::shared_ptr<const COperatorResult> cor =
      GetGlobalResult(
          buildFullname(spectrumModel, stage, method, "extrema_results"))
          .lock();
  std::string type = cor->getType();
  if (type == "PdfCandidatesZResult")
    return std::dynamic_pointer_cast<const PdfCandidatesZResult>(cor)->size();
  if (type == "ExtremaResult")
    return std::dynamic_pointer_cast<const ExtremaResult>(cor)->size();
  if (type == "LineModelExtremaResult")
    return std::dynamic_pointer_cast<const LineModelExtremaResult>(cor)->size();
  if (type == "TplCombinationExtremaResult")
    return std::dynamic_pointer_cast<const TplCombinationExtremaResult>(cor)
        ->size();

  return 0;
}

std::shared_ptr<const CLineModelSolveResult>
COperatorResultStore::GetLineModelSolveResult(const std::string &spectrumModel,
                                              const std::string &stage,
                                              const std::string &method,
                                              const std::string &name) const {

  return std::dynamic_pointer_cast<const CLineModelSolveResult>(
      GetGlobalResult(spectrumModel, stage, method, name).lock());
}

void COperatorResultStore::StoreScopedPerTemplateResult(
    const std::shared_ptr<const CTemplate> &t, const std::string &name,
    std::shared_ptr<const COperatorResult> result) {
  StorePerTemplateResult(t, GetCurrentScopeName(), name, result);
}

void COperatorResultStore::StoreScopedGlobalResult(
    const std::string &name, std::shared_ptr<const COperatorResult> result) {
  StoreGlobalResult(GetCurrentScopeName(), name, result);
}

void COperatorResultStore::StoreGlobalResult(
    const std::string &name, std::shared_ptr<const COperatorResult> result) {
  StoreGlobalResult("", name, result);
}

void COperatorResultStore::StoreScopedFlagResult(const std::string &name) {
  StoreScopedGlobalResult(name, std::make_shared<const CFlagLogResult>(
                                    Flag.getBitMask(), Flag.getListMessages()));
}

std::weak_ptr<const COperatorResult>
COperatorResultStore::GetSolveResult(const std::string &spectrumModel,
                                     const std::string &stage,
                                     const std::string &method) const {
  return GetGlobalResult(spectrumModel, stage, method, "solveResult");
}

bool COperatorResultStore::hasSolveResult(const std::string &spectrumModel,
                                          const std::string &stage,
                                          const std::string &method) const {
  return HasDataset(spectrumModel, stage, method, "solveResult");
}
