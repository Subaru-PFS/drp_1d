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
#ifndef _REDSHIFT_OPERATOR_LINEMODEL_
#define _REDSHIFT_OPERATOR_LINEMODEL_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/linemodel/linemodelfitting.h"
#include "RedshiftLibrary/linemodel/modelrulesresult.h"
#include "RedshiftLibrary/linemodel/velocityfit.h"
#include "RedshiftLibrary/method/twopasssolve.h"
#include "RedshiftLibrary/operator/continuumfitting.h"
#include "RedshiftLibrary/operator/linemodelresult.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/operator/operator.h"
#include "RedshiftLibrary/operator/passextremaresult.h"
#include "RedshiftLibrary/operator/twopass.h"
#include "RedshiftLibrary/photometry/photometricband.h"
#include "RedshiftLibrary/processflow/resultstore.h"

#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/linemodel/linemodelextremaresult.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
namespace Linemodel {
class spanRedshift_test;
class checkSecondPassWindowSize_test;
} // namespace Linemodel

namespace NSEpic {

class CInputContext;
class CTemplateFittingResult;
class COperatorTemplateFittingBase;

/**
 * \ingroup Redshift
 */
class COperatorLineModel : public COperatorTwoPass<TLineModelResult> {

public:
  void Init(const TFloat64List &redshifts, Float64 finestep,
            const bool zLogSampling);

  std::shared_ptr<COperatorResult> getResult();

  void ComputeFirstPass();

  void SetFirstPassCandidates(const TCandidateZbyRank &candidatesz);

  void ComputeSecondPass();

  CLineModelSolution
  computeForLineMeas(std::shared_ptr<const CInputContext> context,
                     const TFloat64List &redshiftsGrid, Float64 &bestZ);

  std::shared_ptr<LineModelExtremaResult> getFirstPassExtremaResults() {
    return m_firstpass_extremaResult;
  }

  std::shared_ptr<LineModelExtremaResult>
  buildExtremaResults(const TCandidateZbyRank &zCandidates,
                      const std::string &opt_continuumreest = "no");

  const bool m_enableWidthFitByGroups = false;
  // m_enableWidthFitByGroups: enable/disable fit by groups. Once enabled, the
  // velocity fitting groups are defined in the line catalog from v4.0 on.

  Int32 m_maxModelSaveCount = 20;
  Float64 m_secondPass_halfwindowsize; // = 0.005;
  std::string m_tplCategory;

  bool m_opt_tplfit_fftprocessing =
      false; // we cant set it as the default since not taken into account when
             // deciding on rebinning
  bool m_opt_tplfit_fftprocessing_secondpass = false; // true;
  bool m_opt_tplfit_use_photometry = false;
  bool m_opt_tplfit_dustFit = true;
  bool m_opt_tplfit_extinction = true;
  Int32 m_opt_fitcontinuum_maxN = 2;
  bool m_opt_tplfit_ignoreLinesSupport =
      false; // default: false, as ortho templates store makes this un-necessary
  bool m_opt_firstpass_multiplecontinuumfit_disable = true;
  bool m_opt_skipsecondpass = false;
  std::string m_opt_firstpass_fittingmethod;
  std::string m_opt_secondpasslcfittingmethod = undefStr;
  TContinuumComponent m_opt_continuumcomponent;
  Float64 m_opt_continuum_neg_amp_threshold = -INFINITY;
  Float64 m_opt_continuum_null_amp_threshold = 0;
  Float64 m_opt_continuum_bad_chi2_threshold = 100;

  EContinuumFit m_continnuum_fit_option = EContinuumFit::retryAll;

  CLineModelSolution
  fitWidthByGroups(std::shared_ptr<const CInputContext> context,
                   Float64 redshift);

  void setHapriorOption(Int32 opt);
  const CSpectrum &
  getFittedModelWithoutcontinuum(const CLineModelSolution &bestModelSolution);

private:
  friend class Linemodel::spanRedshift_test;
  friend class Linemodel::checkSecondPassWindowSize_test;

  std::shared_ptr<CContinuumFitStore>
  PrecomputeContinuumFit(const TFloat64List &redshifts,
                         Int32 candidateIdx = -1);
  void EstimateSecondPassParameters();

  void RecomputeAroundCandidates(
      const std::string &opt_continuumreest, const EContinuumFit tplfit_option,
      const bool overrideRecomputeOnlyOnTheCandidate = false);

  void evaluateAndUpdateContinuumComponent(
      const std::shared_ptr<CContinuumFitStore> &continuumFitStore);

  bool
  updateContinuumComponentIfBadChi2(CContinuumModelSolution const &fitValues);

  bool
  updateContinuumComponentIfNegative(Float64 max_fitamplitudeSigma,
                                     CContinuumModelSolution const &fitValues);

  bool updateContinuumComponentIfNotSignificant(
      Float64 max_fitamplitudeSigma, CContinuumModelSolution const &fitValues);

  std::shared_ptr<CLineModelResult> m_result;
  std::shared_ptr<CLineModelFitting> m_fittingManager;
  TFloat64List m_Redshifts; // coarse grid
  Float64 m_fineStep = NAN;
  bool m_zLogSampling = false;
  Int32 m_estimateLeastSquareFast = 0;
  void fitVelocity(Int32 Zidx, Int32 candidateIdx, Int32 contreest_iterations);

  TFloat64List SpanRedshiftWindow(Float64 z) const;

  Float64 FitBayesWidth(const CSpectrumSpectralAxis &spectralAxis,
                        const CSpectrumFluxAxis &fluxAxis, Float64 z,
                        Int32 start, Int32 end);

  bool AllAmplitudesAreZero(const TBoolList &amplitudesZero);

  bool isfftprocessingActive(Int32 redshiftsTplFitCount);
  void fitContinuumTemplates(
      Int32 candidateIdx, const TFloat64List &redshiftsTplFit,
      std::vector<std::shared_ptr<COperatorResult>> &chisquareResultsAllTpl,
      TStringList &chisquareResultsTplName);
  void getContinuumInfoFromFirstpassFitStore(Int32 candidateIdx,
                                             TInt32List &meiksinIndices,
                                             TInt32List &ebmvIndices,
                                             TTemplateConstRefList &tplList,
                                             bool fft) const;
  void makeContinuumFittingOperator(const TFloat64List &redshifts);
  bool mustReFit(EContinuumFit continuumFit) {
    return continuumFit == EContinuumFit::retryAll ||
           continuumFit == EContinuumFit::reFitFirstPass;
  }
  std::shared_ptr<COperatorContinuumFitting> m_continuumFittingOperator;

  std::shared_ptr<CPriorHelper> m_phelperContinuum;
  std::shared_ptr<CContinuumFitStore> m_tplfitStore_firstpass;
  std::vector<std::shared_ptr<CContinuumFitStore>> m_tplfitStore_secondpass;
  std::vector<VelocityFitSolution> m_velocitySolutions;
};

} // namespace NSEpic

#endif
