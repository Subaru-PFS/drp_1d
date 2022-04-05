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
#include "RedshiftLibrary/linemodel/multirollmodel.h"
#include "RedshiftLibrary/operator/linemodelpassextremaresult.h"
#include "RedshiftLibrary/operator/linemodelresult.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/operator/operator.h"
#include "RedshiftLibrary/operator/spectraFluxResult.h"
#include "RedshiftLibrary/photometry/photometricband.h"
#include "RedshiftLibrary/processflow/resultstore.h"

#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

#include "RedshiftLibrary/linemodel/linemodelextremaresult.h"

namespace NSEpic {

class CInputContext;
/**
 * \ingroup Redshift
 */
class COperatorLineModel {

public:
  Int32 Init(const CSpectrum &spectrum, const TFloat64List &redshifts,
             const CLineCatalog::TLineVector restLineCatalog,
             const TStringList &tplCategoryList,
             const std::string &opt_continuumcomponent,
             const Float64 nsigmasupport, const Float64 halfwdwsize = NAN,
             const Float64 radius = NAN);

  std::shared_ptr<COperatorResult> getResult();

  std::shared_ptr<CTemplatesFitStore> PrecomputeContinuumFit(
      const CSpectrum &spectrum, const CSpectrum &logSampledSpectrum,
      const CTemplateCatalog &tplCatalog, const TFloat64Range &lambdaRange,
      const TFloat64List &redshifts,
      const std::shared_ptr<const CPhotBandCatalog> &photBandCat,
      const Float64 photometry_weight, bool ignoreLinesSupport = false,
      Int32 candidateIdx = -1);

  Int32 ComputeFirstPass(
      const CSpectrum &spectrum, const CSpectrum &logSampledSpc,
      const CTemplateCatalog &tplCatalog,
      const CLineCatalogsTplShape &tplRatioCatalog,
      const TFloat64Range &lambdaRange,
      const std::shared_ptr<const CPhotBandCatalog> &photBandCat,
      const Float64 photo_weight, const std::string &opt_fittingmethod,
      const std::string &opt_lineWidthType, const Float64 opt_velocityEmission,
      const Float64 opt_velocityAbsorption,
      const std::string &opt_continuumreest = "no",
      const std::string &opt_rules = "all",
      const bool &opt_velocityFitting = false,
      const Int32 &opt_twosteplargegridstep_ratio = 10,
      const string &opt_twosteplargegridsampling = "log",
      const std::string &opt_rigidity = "rules",
      const Float64 opt_haprior = -1.);

  void CreateRedshiftLargeGrid(Int32 ratio, TFloat64List &largeGridRedshifts);
  Int32 SetFirstPassCandidates(const TCandidateZbyRank &candidatesz);

  void Combine_firstpass_candidates(
      std::shared_ptr<const CLineModelPassExtremaResult> results_b);

  Int32 ComputeSecondPass(
      const CSpectrum &spectrum, const CSpectrum &logSampledSpectrum,
      const CTemplateCatalog &tplCatalog, const TFloat64Range &lambdaRange,
      const std::shared_ptr<const CPhotBandCatalog> &photBandCat,
      const Float64 photo_weight, const std::string &opt_fittingmethod,
      const std::string &opt_lineWidthType, const Float64 opt_velocityEmission,
      const Float64 opt_velocityAbsorption,
      const std::string &opt_continuumreest = "no",
      const std::string &opt_rules = "all",
      const bool &opt_velocityFitting = false,
      const std::string &opt_rigidity = "rules",
      const Float64 &opt_emvelocityfitmin = 20.,
      const Float64 &opt_emvelocityfitmax = 500.,
      const Float64 &opt_emvelocityfitstep = 20.,
      const Float64 &opt_absvelocityfitmin = 150.,
      const Float64 &opt_absvelocityfitmax = 500.,
      const Float64 &opt_absvelocityfitstep = 20.,
      const string &opt_continuumfit_method = "fromfirstpass");

  Int32 EstimateSecondPassParameters(
      const CSpectrum &spectrum, const TFloat64Range &lambdaRange,
      const std::string &opt_continuumreest, const string &opt_fittingmethod,
      const std::string &opt_rigidity, const bool &opt_velocityFitting,
      const Float64 &opt_emvelocityfitmin, const Float64 &opt_emvelocityfitmax,
      const Float64 &opt_emvelocityfitstep,
      const Float64 &opt_absvelocityfitmin,
      const Float64 &opt_absvelocityfitmax,
      const Float64 &opt_absvelocityfitstep);

  Int32 RecomputeAroundCandidates(
      const TFloat64Range &lambdaRange, const std::string &opt_continuumreest,
      const Int32 tplfit_option,
      const bool overrideRecomputeOnlyOnTheCandidate = false);
  CLineModelSolution
  computeForLineMeas(std::shared_ptr<const CInputContext> context,
                     const TFloat64List &redshiftsGrid, Float64 &bestZ);

  std::shared_ptr<const LineModelExtremaResult>
  saveFirstPassExtremaResults(const TCandidateZbyRank &zCandidates);
  std::shared_ptr<LineModelExtremaResult>
  SaveExtremaResults(const CSpectrum &spectrum,
                     const TFloat64Range &lambdaRange,
                     const TCandidateZbyRank &zCandidates,
                     const std::string &opt_continuumreest = "no");

  void InitTplratioPriors();

  bool m_enableWidthFitByGroups = false;

  Float64 m_linesmodel_nsigmasupport;

  Int32 m_maxModelSaveCount = 20;
  Float64 m_secondPass_halfwindowsize; // = 0.005;
  Float64 m_extremaRedshiftSeparation;
  TStringList m_tplCategoryList;

  bool m_enableLoadContTemplate = false;
  Int32 m_iRollContaminated = -1;
  Float64 m_contLambdaOffset = 0;
  std::shared_ptr<CTemplate> m_tplContaminant = NULL;
  Int32 initContaminant(std::shared_ptr<CModelSpectrumResult> contModelSpectrum,
                        Int32 iRollContaminated, Float64 contLambdaOffset);
  std::shared_ptr<CModelSpectrumResult> GetContaminantSpectrumResult();
  std::shared_ptr<CModelSpectrumResult> m_savedContaminantSpectrumResult;

  bool m_opt_tplfit_fftprocessing =
      false; // we cant set it as the default since not taken into account when
             // deiding on rebinning
  bool m_opt_tplfit_fftprocessing_secondpass = false; // true;
  bool m_opt_tplfit_use_photometry = false;
  Int32 m_opt_tplfit_dustFit = 1;
  Int32 m_opt_tplfit_extinction = 1;
  Int32 m_opt_fitcontinuum_maxN = 2;
  bool m_opt_tplfit_ignoreLinesSupport =
      false; // default: false, as ortho templates store makes this un-necessary
  Float64 m_opt_tplfit_continuumprior_betaA = 1.0;
  Float64 m_opt_tplfit_continuumprior_betaTE = 1.0;
  Float64 m_opt_tplfit_continuumprior_betaZ = 1.0;
  std::string m_opt_tplfit_continuumprior_dirpath = "";

  Int32 m_opt_tplratio_ismFit = 1;
  Int32 m_opt_firstpass_tplratio_ismFit = 0;
  Int32 m_opt_firstpass_multiplecontinuumfit_disable = 1;
  std::string m_opt_firstpass_fittingmethod;
  std::string m_opt_secondpasslcfittingmethod = "-1";
  Float64 m_opt_tplratio_prior_betaA = 1.0;
  Float64 m_opt_tplratio_prior_betaTE = 1.0;
  Float64 m_opt_tplratio_prior_betaZ = 1.0;
  std::string m_opt_tplratio_prior_dirpath = "";
  std::string m_opt_continuumcomponent;
  Float64 m_opt_continuum_neg_amp_threshold = -INFINITY;

  bool m_opt_lya_forcefit;
  bool m_opt_lya_forcedisablefit;
  Float64 m_opt_lya_fit_asym_min;
  Float64 m_opt_lya_fit_asym_max;
  Float64 m_opt_lya_fit_asym_step;
  Float64 m_opt_lya_fit_width_min;
  Float64 m_opt_lya_fit_width_max;
  Float64 m_opt_lya_fit_width_step;
  Float64 m_opt_lya_fit_delta_min;
  Float64 m_opt_lya_fit_delta_max;
  Float64 m_opt_lya_fit_delta_step;

  bool m_opt_enableImproveBalmerFit;
  Int32 m_continnuum_fit_option = 0; // default to "retryall" templates
  // candidates
  std::shared_ptr<CLineModelPassExtremaResult> m_firstpass_extremaResult;
  CLineModelPassExtremaResult m_secondpass_parameters_extremaResult;

  CLineModelSolution
  fitWidthByGroups(std::shared_ptr<const CInputContext> context,
                   Float64 redshift);
  void fitVelocityByGroups(TFloat64List velfitlist, TFloat64List zfitlist,
                           Int32 lineType);

  void setHapriorOption(Int32 opt);
  const CSpectrum &
  getFittedModelWithoutcontinuum(Float64 z,
                                 const CLineModelSolution &bestModelSolution);

private:
  std::shared_ptr<CLineModelResult> m_result;
  std::shared_ptr<CLineModelFitting> m_model;
  CLineCatalog::TLineVector m_RestLineList;
  TFloat64List m_sortedRedshifts;
  Int32 m_enableFastFitLargeGrid = 0;
  Int32 m_estimateLeastSquareFast = 0;
  Float64 m_extremaCount;
  Float64 m_Zlinemeasref;

  TFloat64List SpanRedshiftWindow(Float64 z) const;

  Float64 FitBayesWidth(const CSpectrumSpectralAxis &spectralAxis,
                        const CSpectrumFluxAxis &fluxAxis, Float64 z,
                        Int32 start, Int32 end);

  bool AllAmplitudesAreZero(const TBoolList &amplitudesZero, Int32 nbZ);

  Int32 interpolateLargeGridOnFineGrid(const TFloat64List &redshiftsLargeGrid,
                                       const TFloat64List &redshiftsFineGrid,
                                       const TFloat64List &meritLargeGrid,
                                       TFloat64List &meritFineGrid) const;

  std::shared_ptr<COperatorTemplateFittingBase> m_templateFittingOperator;

  // lmfit

  bool mlmfit_modelInfoSave = false;
  std::vector<std::shared_ptr<CModelSpectrumResult>>
      mlmfit_savedModelSpectrumResults_lmfit;
  std::vector<std::shared_ptr<CLineModelSolution>>
      mlmfit_savedModelFittingResults_lmfit;
  std::vector<std::shared_ptr<CModelRulesResult>>
      mlmfit_savedModelRulesResults_lmfit;
  std::vector<std::shared_ptr<CSpectraFluxResult>>
      mlmfit_savedBaselineResult_lmfit;

  std::shared_ptr<CPriorHelper> m_phelperContinuum;
  std::shared_ptr<CTemplatesFitStore> m_tplfitStore_firstpass;
  std::vector<std::shared_ptr<CTemplatesFitStore>> m_tplfitStore_secondpass;
};

} // namespace NSEpic

#endif
