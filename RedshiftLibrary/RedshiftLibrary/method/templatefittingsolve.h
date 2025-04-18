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
#ifndef _REDSHIFT_METHOD_TEMPLATEFITTINGSOLVE_
#define _REDSHIFT_METHOD_TEMPLATEFITTINGSOLVE_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/method/templatefittingsolveresult.h"
#include "RedshiftLibrary/method/twopasssolve.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/operator/templatefittingBase.h"
#include "RedshiftLibrary/operator/twopass.h"
#include "RedshiftLibrary/processflow/inputcontext.h"
#include "RedshiftLibrary/processflow/resultstore.h"

namespace templateFitting_test {
class fitQuality_test;
}

namespace templateFittingSolve_test {
class computeNoFFT_test;
}
namespace NSEpic {

class CSpectrum;
class CTemplateCatalog;

using TTemplateFittingResultMap =
    std::map<std::string, std::shared_ptr<CTemplateFittingResult>>;
using TConstTemplateFittingResultMap =
    std::map<std::string, std::shared_ptr<const CTemplateFittingResult>>;

/**
 * \ingroup Redshift
 */
class CTemplateFittingSolve : public CTwoPassSolve {

public:
  enum class EType {
    raw,
    continuumOnly,
    noContinuum,
    all,
  };

  static const std::unordered_map<EType, CSpectrum::EType>
      fittingTypeToSpectrumType;

  CTemplateFittingSolve();

private:
  friend class templateFitting_test::fitQuality_test;
  friend class templateFittingSolve_test::computeNoFFT_test;
  void PopulateParameters(
      const std::shared_ptr<const CParameterStore> &parameterStore);
  void InitFittingOperator();
  void LogParameters();
  void CheckTemplateCatalog();
  std::string getResultName() const;
  EType getFitTypeFromParam(const std::string &component);

  std::shared_ptr<CSolveResult> compute() override;

  std::shared_ptr<CTemplateFittingSolveResult> computeSinglePass();

  std::shared_ptr<CTemplateFittingSolveResult> computeTwoPass();
  void computeFirstPass();

  std::shared_ptr<const ExtremaResult>
  computeResults(COperatorPdfz &pdfz, const TZGridListParams &zgridParams = {});
  void storeFirstPassResults(
      const COperatorPdfz &pdfz,
      std::shared_ptr<const ExtremaResult> const &extremaResult);
  void storeResults(const COperatorPdfz &pdfz,
                    std::shared_ptr<const ExtremaResult> const &extremaResult);
  void computeSecondPass(std::shared_ptr<const ExtremaResult> extremaResult);

  std::shared_ptr<CTemplateFittingResult>
  Solve(std::shared_ptr<COperatorResultStore> resultStore,
        const std::shared_ptr<const CTemplate> &tpl, Int32 FitEbmvIdx = allIdx,
        Int32 FitMeiksinIdx = allIdx, std::string parentId = "",
        Int32 candidateIdx = undefIdx,
        std::shared_ptr<CTemplateFittingResult> const &result = nullptr);

  ChisquareArray BuildChisquareArray(const std::string &resultName,
                                     TZGridListParams zgridParams = {}) const;

  std::shared_ptr<ExtremaResult>
  buildExtremaResults(const std::string &resultName,
                      const TCandidateZbyRank &ranked_zCandidates);
  void initSkipSecondPass() override;
  void initTwoPassZStepFactor() override;
  TConstTemplateFittingResultMap
  getPerTemplateResultMap(const std::string &resultName) const;
  TTemplateFittingResultMap
  getPerTemplateResultMapCopy(const std::string &resultName) const;

  std::shared_ptr<COperatorTemplateFittingBase> m_templateFittingOperator;
  std::string m_opt_pdfcombination;
  Float64 m_redshiftSeparation;
  Int32 m_opt_maxCandidate;
  Float64 m_overlapThreshold;
  EType m_spectrumType = EType::raw;
  std::string m_interpolation;
  bool m_extinction;
  bool m_dustFit;
  bool m_fftProcessing;
  bool m_usePhotometry = false;
  Float64 m_photometryWeight = NAN;
  EContinuumFit m_secondPassContinuumFit = EContinuumFit::undefined;
  Float64 m_secondPass_halfwindowsize = NAN;
  Int32 m_opt_extremacount;
  Float64 m_opt_candidatesLogprobaCutThreshold = 0.0;
  bool m_isFirstPass = true;
};

} // namespace NSEpic

#endif
