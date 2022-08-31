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
#ifndef _REDSHIFT_LINEMODEL_FITTING_
#define _REDSHIFT_LINEMODEL_FITTING_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/line/catalogsTplRatio.h"
#include "RedshiftLibrary/line/regulament.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

#include "RedshiftLibrary/linemodel/element.h"
#include "RedshiftLibrary/operator/linemodelresult.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"

#include "RedshiftLibrary/operator/pdfz.h"

#include <boost/shared_ptr.hpp>

#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/linemodel/abstractfitter.h"
#include "RedshiftLibrary/linemodel/continuummanager.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/linemodel/spectrummodel.h"

#include <memory>
namespace NSEpic {

// to avoid circular dependency
class CLineProfileASYM;
class CLineProfileASYMFIXED;

class CLineModelFitting {

public:
  CLineModelFitting();
  CLineModelFitting(
      const std::shared_ptr<const CSpectrum> &template_,
      const TLambdaRange
          &lambdaRange); // only used for template orthogonalization, TODO use
                         // only one of the future subclasses ? at least inherit
                         // from clinemodelfitting

  void initParameters();
  void initMembers();
  void LoadCatalog(const CLineCatalog::TLineVector &restLineList);
  void LoadCatalogOneMultiline(const CLineCatalog::TLineVector &restLineList);
  void
  LoadCatalogTwoMultilinesAE(const CLineCatalog::TLineVector &restLineList);

  void LogCatalogInfos();

  void setRedshift(Float64 redshift, bool reinterpolatedContinuum);
  void SetContinuumComponent(std::string component);

  std::shared_ptr<const CContinuumManager> getContinuumManager() const {
    return m_continuumManager;
  }
  std::shared_ptr<CContinuumManager> getContinuumManager() {
    return m_continuumManager;
  }

  bool initDtd();
  Float64 EstimateDTransposeD(const std::string &spcComponent) const;
  Float64 EstimateMTransposeM() const;
  Float64 EstimateLikelihoodCstLog() const;
  Float64 getDTransposeD();
  Float64 getLikelihood_cstLog();

  const std::string &getTplratio_bestTplName() const;
  Float64 getTplratio_bestTplIsmCoeff() const;
  Float64 getTplratio_bestAmplitude() const;
  Float64 getTplratio_bestDtm() const;
  Float64 getTplratio_bestMtm() const;
  Int32 getTplratio_count() const;
  const TFloat64List &getTplratio_priors();
  const TFloat64List &GetChisquareTplratio() const;
  TFloat64List GetPriorLinesTplratio() const;
  const TFloat64List &GetScaleMargTplratio() const;
  const TBoolList &GetStrongELPresentTplratio() const;
  const TBoolList &getHaELPresentTplratio() const;
  const TInt32List &GetNLinesAboveSNRTplratio() const;
  void SetTplratio_PriorHelper();

  Int32 GetNElements() const;

  void SetVelocityEmission(Float64 vel);
  void SetVelocityAbsorption(Float64 vel);
  void SetVelocityEmissionOneElement(Float64 vel, Int32 idxElt);
  void SetVelocityAbsorptionOneElement(Float64 vel, Int32 idxElt);

  void setVelocity(Float64 vel, Int32 lineType);
  void setVelocity(Float64 vel, Int32 idxElt, Int32 lineType);

  Float64 GetVelocityEmission() const;
  Float64 GetVelocityAbsorption() const;
  Int32 ApplyVelocityBound(Float64 inf, Float64 sup);

  bool initModelAtZ(Float64 redshift,
                    const CSpectrumSpectralAxis &spectralAxis);

  Float64 fit(Float64 redshift, CLineModelSolution &modelSolution,
              CContinuumModelSolution &continuumModelSolution,
              Int32 contreest_iterations = 0, bool enableLogging = 0);
  TFloat64Range &getLambdaRange() { return m_dTransposeDLambdaRange; };
  void initTplratioCatalogs(Int32 opt_tplratio_ismFit);

  bool setTplratioModel(Int32 itplratio, bool enableSetVelocity = false);
  bool setTplratioAmplitude(const TFloat64List &ampsElts,
                            const TFloat64List &errorsElts);

  void SetFittingMethod(const std::string &fitMethod);
  void SetSecondpassContinuumFitPrms();

  void SetAbsLinesLimit(Float64 limit);
  void SetLeastSquareFastEstimationEnabled(Int32 enabled);

  Float64 GetRedshift() const;

  CSpectrumFluxAxis getModel(Int32 lineTypeFilter = -1) const;

  CMask getOutsideLinesMask() const;
  Float64 getOutsideLinesSTD(Int32 which) const;

  Int32 getSpcNSamples() const;
  Float64 getLeastSquareMeritFast(Int32 idxLine = -1) const;
  Float64 getLeastSquareContinuumMeritFast() const;
  Float64 getLeastSquareMerit() const;
  Float64 getLeastSquareContinuumMerit() const;
  Float64 getLeastSquareMeritUnderElements() const;
  Float64 getScaleMargCorrection(Int32 idxLine = -1) const;

  Float64 getStrongerMultipleELAmpCoeff() const;
  TStringList getLinesAboveSNR(Float64 snrcut = 3.5) const;
  Float64 getCumulSNRStrongEL() const;
  Float64 getCumulSNROnRange(TInt32Range idxRange) const;
  bool GetModelStrongEmissionLinePresent() const;
  bool GetModelHaStrongest() const;

  void LoadModelSolution(const CLineModelSolution &modelSolution);
  CLineModelSolution GetModelSolution(Int32 opt_level = 0);

  void getFluxDirectIntegration(const TInt32List &eIdx_list,
                                const TInt32List &subeIdx_list,
                                bool substract_abslinesmodel, Float64 &fluxdi,
                                Float64 &snrdi) const;

  Float64 getModelFluxVal(Int32 idx) const;
  void logParameters();
  CLineModelElementList m_Elements;
  std::shared_ptr<CAbstractFitter> m_fitter;

  std::shared_ptr<const CSpectrum> m_inputSpc;
  const CLineCatalog::TLineVector m_RestLineList;

  const TStringList &GetModelRulesLog() const;

  Int32 setPassMode(Int32 iPass);
  Int32 GetPassNumber() const;

  void SetForcedisableTplratioISMfit(bool opt);
  void prepareAndLoadContinuum(Int32 icontfitting, Float64 redshift);
  void computeSpectrumFluxWithoutContinuum();
  bool isContinuumComponentTplfitxx() const {
    return m_ContinuumComponent == "tplfit" ||
           m_ContinuumComponent == "tplfitauto";
  }
  void duplicateTplratioResult(Int32 ifitting, TFloat64List &bestTplratioMerit,
                               TFloat64List &bestTplratioMeritPrior);
  void updateTplratioResults(Int32 ifitting, Float64 _merit,
                             Float64 _meritprior,
                             TFloat64List &bestTplratioMerit,
                             TFloat64List &bestTplratioMeritPrior);
  Float64 computelogLinePriorMerit(
      Int32 itratio,
      const std::vector<CPriorHelper::SPriorTZE> &logPriorDataTplRatio);
  std::shared_ptr<CSpectrumModel> getSpectrumModel() {
    return m_model;
  } // not const because of tplortho
  CLineCatalogsTplRatio m_CatalogTplRatio;
  TFloat64List m_ChisquareTplratio;
  std::vector<TFloat64List> m_FittedAmpTplratio;
  std::vector<TFloat64List> m_FittedErrorTplratio;
  std::vector<TFloat64List> m_MtmTplratio;
  std::vector<TFloat64List> m_DtmTplratio;
  std::vector<TFloat64List> m_LyaAsymCoeffTplratio;
  std::vector<TFloat64List> m_LyaWidthCoeffTplratio;
  std::vector<TFloat64List> m_LyaDeltaCoeffTplratio;
  std::vector<TInt32List> m_LyaIgmIdxTplratio;
  std::vector<TFloat64List> m_LinesLogPriorTplratio;

  Int32 m_pass = 1;
  bool m_enableAmplitudeOffsets;

  Float64 m_LambdaOffsetMin = -400.0;
  Float64 m_LambdaOffsetMax = 400.0;
  Float64 m_LambdaOffsetStep = 25.0;
  bool m_enableLambdaOffsetsFit;

  bool m_opt_firstpass_forcedisableTplratioISMfit = true;

private:
  void fitAmplitudesSimplex();

  bool m_forceDisableLyaFitting = false;
  bool m_forceLyaFitting = false;

  bool SetMultilineNominalAmplitudesFast(Int32 iCatalog);
  void setLyaProfile(Float64 redshift,
                     const CLineCatalog::TLineVector &lineList,
                     bool tplratio = false);
  void setAsymProfile(Int32 idxLyaE, Int32 idxLineLyaE, Float64 redshift,
                      const CLineCatalog::TLineVector &lineList,
                      bool tplratio = false);
  void setSymIgmProfile(Int32 iElts, const TInt32List &idxLineIGM,
                        Float64 redshift);
  TAsymParams fitAsymParameters(Float64 redshift, Int32 idxLyaE,
                                const Int32 &idxLineLyaE);
  Int32 fitAsymIGMCorrection(Float64 redshift, Int32 idxLyaE,
                             const TInt32List &idxLine);
  Int32 getLineIndexInCatalog(Int32 iElts, Int32 idxLine,
                              const CLineCatalog::TLineVector &catalog,
                              bool tplratio) const;
  void SetLSF();

  TInt32List getOverlappingElementsBySupport(Int32 ind,
                                             Float64 overlapThres = 0.1) const;

  TInt32List findLineIdxInCatalog(const CLineCatalog::TLineVector &restLineList,
                                  const std::string &strTag, Int32 type) const;
  TPolynomCoeffs getPolynomCoeffs(Int32 eIdx) const;
  void applyPolynomCoeffs(Int32 eIdx, const TPolynomCoeffs &polynom_coeffs);
  void addDoubleLine(const CLine &r1, const CLine &r2, Int32 index1,
                     Int32 index2, Float64 nominalWidth, Float64 a1,
                     Float64 a2);

  void applyRules(bool enableLogs = false);
  std::vector<CRange<Int32>>
  getlambdaIndexesUnderLines(const TInt32List &eIdx_list,
                             const TInt32List &subeIdx_list,
                             Float64 sigma_support) const;
  void
  integrateFluxes_usingTrapez(const CSpectrumFluxAxis &continuumFlux,
                              const std::vector<CRange<Int32>> &indexRangeList,
                              Float64 &sumFlux, Float64 &sumErr) const;

  CRegulament m_Regulament;
  std::shared_ptr<CSpectrumModel> m_model;
  std::shared_ptr<CContinuumManager> m_continuumManager;

  TFloat64List m_ScaleMargCorrTplratio;
  TBoolList m_StrongELPresentTplratio;
  TBoolList m_StrongHalphaELPresentTplratio;
  TInt32List m_NLinesAboveSNRTplratio;

  Float64 m_Redshift;

  Float64 m_dTransposeD; // the cached dtd (maximum chisquare value)
  TFloat64Range m_dTransposeDLambdaRange; // the lambdaRange used to computed
                                          // cached dTransposeD values
  Float64 m_likelihood_cstLog; // constant term for the Likelihood calculation

  // Float64* m_unscaleContinuumFluxAxisDerivZ;

  std::string m_ContinuumComponent;
  std::string m_LineWidthType;
  Float64 m_NSigmaSupport;
  Float64 m_velocityEmission;
  Float64 m_velocityAbsorption;
  Float64 m_velocityEmissionInit;
  Float64 m_velocityAbsorptionInit;

  Float64 m_nominalWidthDefaultEmission;
  Float64 m_nominalWidthDefaultAbsorption;

  std::string m_fittingmethod;

  std::string m_rulesoption;
  std::string m_rigidity;
  bool m_forcedisableTplratioISMfit = false;

  std::string m_tplratioBestTplName = "undefined";
  Float64 m_tplratioBestTplIsmCoeff = NAN;
  Float64 m_tplratioBestTplAmplitude = NAN;
  Float64 m_tplratioBestTplDtm = NAN;
  Float64 m_tplratioBestTplMtm = NAN;
  Int32 m_tplratioLeastSquareFast =
      0; // for rigidity=tplratio: switch to use fast least square estimation
  std::shared_ptr<CPriorHelper> m_tplratio_priorhelper;

  Int32 m_secondpass_fitContinuum_dustfit;
  Int32 m_secondpass_fitContinuum_igm;

  bool m_forcedisableMultipleContinuumfit = false;

  bool m_lmfit_noContinuumTemplate;
  bool m_lmfit_bestTemplate;
  bool m_lmfit_fitContinuum;
  bool m_lmfit_fitEmissionVelocity;
  bool m_lmfit_fitAbsorptionVelocity;

  std::shared_ptr<const TFloat64Range> m_lambdaRange;

  linetags ltags;

  bool m_opt_lya_forcefit = false;
  bool m_opt_lya_forcedisablefit = false;
  Float64 m_opt_lya_fit_asym_min = 0.0;
  Float64 m_opt_lya_fit_asym_max = 4.0;
  Float64 m_opt_lya_fit_asym_step = 1.0;
  Float64 m_opt_lya_fit_width_min = 1.;
  Float64 m_opt_lya_fit_width_max = 4.;
  Float64 m_opt_lya_fit_width_step = 1.;
  Float64 m_opt_lya_fit_delta_min = 0.;
  Float64 m_opt_lya_fit_delta_max = 0.;
  Float64 m_opt_lya_fit_delta_step = 1.;

  bool m_opt_firstpass_forcedisableMultipleContinuumfit = true;

  std::string m_opt_firstpass_fittingmethod = "hybrid";
  std::string m_opt_secondpass_fittingmethod = "hybrid";
  bool m_ignoreLinesSupport = false;

  //  bool m_opt_enable_improveBalmerFit = false;
  Float64 m_opt_haprior = -1.;
  bool m_useloglambdasampling = false;
};

} // namespace NSEpic

#endif // _REDSHIFT_LINEMODEL_ELEMENTLIST_
