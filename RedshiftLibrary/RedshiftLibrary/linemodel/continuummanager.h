#ifndef _REDSHIFT_CONTINUUM_MANAGER_H
#define _REDSHIFT_CONTINUUM_MANAGER_H

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"

#include <cmath>

namespace NSEpic

{
class CTemplate;
class CMask;
class CPriorHelper;
class CTemplatesFitStore;
class COperatorTemplateFitting;
class CSpectrum;
class CTemplateCatalog;
class CSpectrumFluxCorrectionCalzetti;
class CSpectrumModel;
class CContinuumModelSolution;
class CContinuumManager {
public:
  CContinuumManager(const std::shared_ptr<const CSpectrum> &spc,
                    const std::shared_ptr<CSpectrumModel> &model,
                    const TLambdaRange &lambdaRange);

  bool SolveContinuum(const std::shared_ptr<const CTemplate> &tpl,
                      const TFloat64List &redshifts, Float64 overlapThreshold,
                      std::vector<CMask> maskList, std::string opt_interp,
                      Int32 opt_extinction, Int32 opt_dustFit, Float64 &merit,
                      Float64 &fitAmplitude, Float64 &fitAmplitudeError,
                      Float64 &fitAmplitudeSigma, Float64 &fitEbmvCoeff,
                      Int32 &fitMeiksinIdx, Float64 &fitDtM, Float64 &fitMtM,
                      Float64 &fitLogprior);

  const std::string &getFitContinuum_tplName() const;
  Float64 getFitContinuum_tplAmplitude() const;
  Float64 getFitContinuum_tplAmplitudeError() const;
  Float64 getFitContinuum_snr() const;
  Float64 getFitContinuum_tplMerit() const;
  Float64 getFitContinuum_tplMeritPhot() const;
  void setFitContinuum_tplAmplitude(Float64 tplAmp, Float64 tplAmpErr,
                                    const TFloat64List &polyCoeffs);
  Float64 getFitContinuum_tplIsmEbmvCoeff() const;
  Float64 getFitContinuum_tplIgmMeiksinIdx() const;
  Int32 SetFitContinuum_FitStore(
      const std::shared_ptr<const CTemplatesFitStore> &fitStore);
  void SetFitContinuum_SNRMax(Float64 snr_max);
  void SetFitContinuum_Option(Int32 opt);
  Int32 GetFitContinuum_Option() const;

  const std::shared_ptr<const CTemplatesFitStore> &
  GetFitContinuum_FitStore() const;
  std::shared_ptr<CPriorHelper> SetFitContinuum_PriorHelper();
  void LoadFitContinuum(Int32 icontinuum, Int32 autoSelect, Float64 redshift);

  void SetFitContinuum_FitValues(std::string tplfit_name, Float64 tplfit_amp,
                                 Float64 tplfit_amperr, Float64 tplfit_chi2,
                                 Float64 tplfit_chi2_phot, Float64 tplfit_ebmv,
                                 Int32 tplfit_meiksinidx,
                                 Float64 tplfit_continuumredshift,
                                 Float64 tplfit_dtm, Float64 tplfit_mtm,
                                 Float64 tplfit_logprior,
                                 const TFloat64List &polyCoeffs);
  Int32 ApplyContinuumOnGrid(const std::shared_ptr<const CTemplate> &tpl,
                             Float64 zcontinuum);
  Float64 getContinuumScaleMargCorrection() const;
  bool isContinuumComponentTplfitxx() const {
    return m_ContinuumComponent == "tplfit" ||
           m_ContinuumComponent == "tplfitauto";
  }
  CContinuumModelSolution GetContinuumModelSolution() const;
  void setContinuumComponent(std::string component);

  // new methods
  void logParameters();
  std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
  getIsmCorrectionFromTpl();
  void reinterpolateContinuum(Float64 redshift);
  void reinterpolateContinuum();
  void initObserveGridContinuumFlux(Int32 size);
  bool isContFittedToNull();
  Int32 getFittedMeiksinIndex() { return m_fitContinuum_tplFitMeiksinIdx; }
  Float64 getFitSum() {
    if (!isContinuumComponentTplfitxx())
      return 0.0;
    return m_fitContinuum_tplFitMerit_phot +
           m_fitContinuum_tplFitLogprior; // unconditionnal sum (if photometry
                                          // disabled, will sum 0.0)
  }
  Float64 getTerm1() {
    return m_fitContinuum_tplFitAmplitude * m_fitContinuum_tplFitAmplitude *
           m_fitContinuum_tplFitMtM;
  }
  Float64 getTerm2() {
    return -2. * m_fitContinuum_tplFitAmplitude * m_fitContinuum_tplFitDtM;
  }
  Float64 getFittedLogPrior() { return m_fitContinuum_tplFitLogprior; }
  void setFitContinuumFromFittedAmps(TFloat64List &ampsfitted,
                                     TInt32List &validEltsIdx);

  Int32 m_opt_fitcontinuum_maxCount = 2;

private:
  std::shared_ptr<const CTemplateCatalog> m_tplCatalog;
  TStringList m_tplCategoryList;

  std::shared_ptr<CPriorHelper> m_fitContinuum_priorhelper;
  std::shared_ptr<COperatorTemplateFitting> m_templateFittingOperator;

  std::shared_ptr<const CTemplatesFitStore> m_fitContinuum_tplfitStore;

  std::shared_ptr<CSpectrumModel> m_model;

  std::string m_ContinuumComponent;
  Int32 m_fitContinuum_option;
  Float64 m_opt_fitcontinuum_neg_threshold = -INFINITY;
  Float64 m_opt_fitcontinuum_null_amp_threshold = 0.;

  std::string m_fitContinuum_tplName;
  Float64 m_fitContinuum_tplFitAmplitude = NAN;
  Float64 m_fitContinuum_tplFitAmplitudeError = NAN;
  Float64 m_fitContinuum_tplFitAmplitudeSigmaMAX = NAN;
  Float64 m_fitContinuum_tplFitMerit = NAN;
  Float64 m_fitContinuum_tplFitMerit_phot = NAN;
  Float64 m_fitContinuum_tplFitEbmvCoeff = NAN;
  Int32 m_fitContinuum_tplFitMeiksinIdx = -1;
  Float64 m_fitContinuum_tplFitRedshift =
      NAN; // only used with m_fitContinuum_option==2 for now
  Float64 m_fitContinuum_tplFitDtM = NAN;
  Float64 m_fitContinuum_tplFitMtM = NAN;
  Float64 m_fitContinuum_tplFitLogprior = NAN;
  Float64 m_fitContinuum_tplFitSNRMax = NAN;
  TFloat64List m_fitContinuum_tplFitPolyCoeffs; // only used with
  // m_fitContinuum_option==2 for now
  Float64 m_fitContinuum_tplFitAlpha = 0.;

  TAxisSampleList
      m_observeGridContinuumFlux; // the continuum spectre without the
                                  // amplitude coeff; m_ContinuumFLux = amp *
                                  // m_observeGridContinuumFlux
};

} // namespace NSEpic
#endif
