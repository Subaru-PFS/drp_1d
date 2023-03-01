#ifndef _REDSHIFT_CONTINUUM_MANAGER_H
#define _REDSHIFT_CONTINUUM_MANAGER_H

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/linemodel/templatesfitstore.h"

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
  CContinuumManager(const std::shared_ptr<CSpectrumModel> &model,
                    const TLambdaRange &lambdaRange,
                    std::shared_ptr<CContinuumModelSolution>);

  Int32 SetFitContinuum_FitStore(
      const std::shared_ptr<const CTemplatesFitStore> &fitStore);
  void SetFitContinuum_SNRMax(Float64 snr_max);
  void SetFitContinuum_Option(Int32 opt);
  Int32 GetFitContinuum_Option() const;

  const std::shared_ptr<const CTemplatesFitStore> &
  GetFitContinuum_FitStore() const;
  std::shared_ptr<CPriorHelper> SetFitContinuum_PriorHelper();
  void LoadFitContinuum(Int32 icontinuum, Int32 autoSelect, Float64 redshift);

  void SetFitContinuum_FitValues(const CContinuumModelSolution &cms) {
    *m_fitContinuum = cms;
  }

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
  Float64 getFitContinuum_snr() const;
  CContinuumModelSolution GetContinuumModelSolutionCopy() const;
  std::shared_ptr<const CContinuumModelSolution>
  GetContinuumModelSolution() const {
    return m_fitContinuum;
  }
  void setContinuumComponent(std::string component);
  const std::string &getContinuumComponent() const {
    return m_ContinuumComponent;
  };

  // new methods
  void logParameters();
  std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
  getIsmCorrectionFromTpl();
  void reinterpolateContinuum(Float64 redshift);
  void reinterpolateContinuumResetAmp();
  void initObserveGridContinuumFlux(Int32 size);
  bool isContFittedToNull();
  Int32 getFittedMeiksinIndex() { return m_fitContinuum->tplMeiksinIdx; }
  Float64 getFitSum() {
    if (!isContinuumComponentTplfitxx())
      return 0.0;
    return m_fitContinuum->tplMeritPhot +
           m_fitContinuum->tplLogPrior; // unconditionnal sum (if photometry
                                        // disabled, will sum 0.0)
  }
  Float64 getTerm1() {
    return m_fitContinuum->tplAmplitude * m_fitContinuum->tplAmplitude *
           m_fitContinuum->tplMtM;
  }
  Float64 getTerm2() {
    return -2. * m_fitContinuum->tplAmplitude * m_fitContinuum->tplDtM;
  }
  Float64 getFittedLogPrior() { return m_fitContinuum->tplLogPrior; }
  void setFitContinuumFromFittedAmps(TFloat64List &ampsfitted,
                                     TInt32List &validEltsIdx);

private:
  std::shared_ptr<const CTemplateCatalog> m_tplCatalog;
  TStringList m_tplCategoryList;

  Int32 m_spcIndex = 0;
  std::shared_ptr<CPriorHelper> m_fitContinuum_priorhelper;
  std::shared_ptr<COperatorTemplateFitting> m_templateFittingOperator;

  std::shared_ptr<const CTemplatesFitStore> m_fitContinuum_tplfitStore;

  std::shared_ptr<CSpectrumModel> m_model;

  std::string m_ContinuumComponent;
  Int32 m_fitContinuum_option;
  Float64 m_opt_fitcontinuum_neg_threshold = -INFINITY;
  Float64 m_opt_fitcontinuum_null_amp_threshold = 0.;

  std::shared_ptr<CContinuumModelSolution> m_fitContinuum;
  std::shared_ptr<fitMaxValues> m_fitContinuumMaxValues;

  // m_fitContinuum_option==2 for now
  Float64 m_fitContinuum_tplFitAlpha = 0.;

  TAxisSampleList
      m_observeGridContinuumFlux; // the continuum spectre without the
                                  // amplitude coeff; m_ContinuumFLux = amp *
                                  // m_observeGridContinuumFlux

  void setFitContinuum_tplAmplitude(Float64 tplAmp, Float64 tplAmpErr,
                                    const TFloat64List &polyCoeffs);
};

} // namespace NSEpic
#endif
