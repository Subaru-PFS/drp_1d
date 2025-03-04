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
#ifndef _REDSHIFT_PROCESSFLOW_OPERATORRESULTSTORE_
#define _REDSHIFT_PROCESSFLOW_OPERATORRESULTSTORE_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/processflow/result.h"
#include "RedshiftLibrary/processflow/scopestore.h"

namespace ResultStore { // boost_test_suite
// all boost_auto_test_case that use private method
class CreateResultStorage_test;
class StoreResult_test;
class StoreTemplateMethods_test;
class GetMethods_test;
} // namespace ResultStore

namespace NSEpic {

class CTemplate;
class CClassificationResult;
class CReliabilityResult;
class CLogZPdfResult;
class TLineModelResult;
class TTplCombinationResult;
class TExtremaResult;
class CModelSpectrumResult;
class CModelPhotValueResult;
class CLineModelSolution;
class CFlagLogResult;
template <class T = TLineModelResult> class CLineModelExtremaResult;
template <class T = TTplCombinationResult> class CTplCombinationExtremaResult;
//  class CLineModelExtremaResult<TLineModelResult>;
// class LineModelExtremaResult;
/**
 * \ingroup Redshift
 */
using TResultsMap =
    std::map<std::string, std::shared_ptr<const COperatorResult>>;
using TPerTemplateResultsMap = std::map<std::string, TResultsMap>;

class COperatorResultStore : public CScopeStore {

public:
  static std::string buildFullname(const std::string &spectrumModel,
                                   const std::string &stage,
                                   const std::string &method,
                                   const std::string &name);

  COperatorResultStore(const std::shared_ptr<const CScopeStack> &scopeStack);

  void StorePerTemplateResult(const std::shared_ptr<const CTemplate> &t,
                              const std::string &path, const std::string &name,
                              std::shared_ptr<const COperatorResult> result);
  void StoreGlobalResult(const std::string &path, const std::string &name,
                         std::shared_ptr<const COperatorResult> result);
  std::weak_ptr<const COperatorResult>
  GetScopedPerTemplateResult(const std::shared_ptr<const CTemplate> &t,
                             const std::string &name) const;
  TResultsMap GetScopedPerTemplateResult(const std::string &name) const;
  std::weak_ptr<const COperatorResult>
  GetScopedGlobalResult(const std::string &name) const;

  const std::string &GetGlobalResultType(const std::string &spectrumModel,
                                         const std::string &stage,
                                         const std::string &method,
                                         const std::string &name) const;

  const std::string &GetCandidateResultType(const std::string &spectrumModel,
                                            const std::string &stage,
                                            const std::string &method,
                                            const std::string &name,
                                            const std::string &dataset) const;

  bool HasCandidateDataset(const std::string &spectrumModel,
                           const std::string &stage, const std::string &method,
                           const std::string &name,
                           const std::string &dataset) const;

  bool HasDataset(const std::string &spectrumModel, const std::string &stage,
                  const std::string &method, const std::string &name) const;
  bool hasContextWarningFlag() const;
  bool hasInitWarningFlag() const;
  bool hasCurrentMethodWarningFlag() const;

  std::shared_ptr<const CClassificationResult>
  GetClassificationResult(const std::string &spectrumModel,
                          const std::string &stage, const std::string &method,
                          const std::string &name) const;

  std::shared_ptr<const CReliabilityResult>
  GetReliabilityResult(const std::string &spectrumModel,
                       const std::string &stage, const std::string &method,
                       const std::string &name) const;

  std::shared_ptr<const CLogZPdfResult>
  GetLogZPdfResult(const std::string &spectrumModel, const std::string &stage,
                   const std::string &method, const std::string &name) const;

  std::shared_ptr<const CFlagLogResult>
  GetFlagLogResult(const std::string &spectrumModel, const std::string &stage,
                   const std::string &method, const std::string &name) const;

  std::shared_ptr<const TLineModelResult>
  GetLineModelResult(const std::string &spectrumModel, const std::string &stage,
                     const std::string &method, const std::string &name,
                     const std::string &dataset, const int &rank,
                     bool firstpassCorrespondingResult = false) const;

  std::shared_ptr<const TTplCombinationResult>
  GetTplCombinationResult(const std::string &spectrumModel,
                          const std::string &stage, const std::string &method,
                          const std::string &name, const std::string &dataset,
                          const int &rank) const;
  std::shared_ptr<const TExtremaResult>
  GetExtremaResult(const std::string &spectrumModel, const std::string &stage,
                   const std::string &method, const std::string &name,
                   const std::string &dataset, const int &rank,
                   bool firstpassCorrespondingResult = false) const;

  std::shared_ptr<const CLineModelSolution>
  GetLineModelSolution(const std::string &spectrumModel,
                       const std::string &stage, const std::string &method,
                       const std::string &name) const;

  std::shared_ptr<const CLineModelSolution>
  GetLineModelSolution(const std::string &spectrumModel,
                       const std::string &stage, const std::string &method,
                       const std::string &name, const std::string &dataset,
                       const int &rank) const;

  std::shared_ptr<const CModelSpectrumResult>
  GetModelSpectrumResult(const std::string &spectrumModel,
                         const std::string &stage, const std::string &method,
                         const std::string &name, const std::string &dataset,
                         const int &rank) const;
  std::shared_ptr<const CModelSpectrumResult>
  GetModelSpectrumResult(const std::string &spectrumModel,
                         const std::string &stage, const std::string &method,
                         const std::string &name) const;

  std::shared_ptr<const CModelPhotValueResult>
  GetModelPhotValueResult(const std::string &spectrumModel,
                          const std::string &stage, const std::string &method,
                          const std::string &name, const std::string &dataset,
                          const int &rank) const;
  int getNbRedshiftCandidates(const std::string &spectrumModel,
                              const std::string &stage,
                              const std::string &method) const;

  std::weak_ptr<const COperatorResult>
  GetSolveResult(const std::string &spectrumModel, const std::string &stage,
                 const std::string &method) const;
  bool hasSolveResult(const std::string &spectrumModel,
                      const std::string &stage,
                      const std::string &method) const;
  // From DataStore, above should be removed and integrated into these
  void StoreGlobalResult(const std::string &name,
                         std::shared_ptr<const COperatorResult> result);

  void
  StoreScopedPerTemplateResult(const std::shared_ptr<const CTemplate> &t,
                               const std::string &name,
                               std::shared_ptr<const COperatorResult> result);
  void StoreScopedGlobalResult(const std::string &name,
                               std::shared_ptr<const COperatorResult> result);
  void StoreScopedFlagResult(const std::string &name);

  void reset() {
    m_GlobalResults.clear();
    m_PerTemplateResults.clear();
  }

private:
  friend class ResultStore::StoreTemplateMethods_test;
  friend class ResultStore::GetMethods_test;

  std::weak_ptr<const COperatorResult>
  GetPerTemplateResult(const std::shared_ptr<const CTemplate> &t,
                       const std::string &name) const;
  TResultsMap GetPerTemplateResult(const std::string &name) const;
  std::weak_ptr<const COperatorResult>
  GetGlobalResult(const std::string &name) const;
  std::weak_ptr<const COperatorResult>
  GetGlobalResult(const std::string &spectrumModel, const std::string &stage,
                  const std::string &method, const std::string &name) const;

protected:
  friend class ResultStore::CreateResultStorage_test;
  friend class ResultStore::StoreResult_test;

  void StoreResult(TResultsMap &map, const std::string &path,
                   const std::string &name,
                   std::shared_ptr<const COperatorResult> result);

  TPerTemplateResultsMap m_PerTemplateResults;
  TResultsMap m_GlobalResults;
};

} // namespace NSEpic

#endif
