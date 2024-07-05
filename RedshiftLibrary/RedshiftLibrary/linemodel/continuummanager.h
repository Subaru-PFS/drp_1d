#ifndef _REDSHIFT_CONTINUUM_MANAGER_H
#define _REDSHIFT_CONTINUUM_MANAGER_H

#include <cmath>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/linemodel/obsiterator.h"
#include "RedshiftLibrary/linemodel/spectrummodel.h"
#include "RedshiftLibrary/linemodel/templatesfitstore.h"

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

class CContinuumModelSolution;

class CContinuumManager {
public:
  CContinuumManager(const CSpcModelVectorPtr &models,
                    std::shared_ptr<CContinuumModelSolution>,
                    const CSpectraGlobalIndex &spcGlobIndex);

  const CSpectrumModel &getModel() const {
    return m_models->getSpectrumModel();
  }
  CSpectrumModel &getModel() { return m_models->getSpectrumModel(); }

  Int32 SetFitContinuum_FitStore(
      const std::shared_ptr<const CContinuumFitStore> &fitStore);
  void SetFitContinuum_SNRMax(Float64 snr_max);
  void SetFitContinuum_Option(Int32 opt);
  Int32 GetFitContinuum_Option() const;

  const std::shared_ptr<const CContinuumFitStore> &
  GetFitContinuum_FitStore() const;
  std::shared_ptr<CPriorHelper> SetFitContinuum_PriorHelper();
  void LoadFitContinuum(Int32 icontinuum, Float64 redshift);

  void SetFitContinuum_FitValues(const CContinuumModelSolution &cms) {
    *m_fitContinuum = cms;
  }

  Float64 getContinuumScaleMargCorrection() const;
  bool isContinuumComponentTplfitxx() const {
    return m_ContinuumComponent == "tplFit" ||
           m_ContinuumComponent == "tplFitAuto";
  }
  // TODO use wherever possible
  bool isContinuumComponentPowerLaw() const {
    return m_ContinuumComponent == "powerLaw";
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

  void logParameters();
  std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
  getIsmCorrectionFromTpl();
  void reinterpolateContinuum(Float64 redshift);
  void reinterpolateContinuumResetAmp();

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
  std::string m_tplCategory;

  std::shared_ptr<CPriorHelper> m_fitContinuum_priorhelper;

  std::shared_ptr<const CContinuumFitStore> m_fitContinuum_tplfitStore;

  CSpcModelVectorPtr m_models;

  CSpectraGlobalIndex m_spectraIndex;

  std::string m_ContinuumComponent;
  Int32 m_fitContinuum_option;
  Float64 m_opt_fitcontinuum_neg_threshold = -INFINITY;
  Float64 m_opt_fitcontinuum_null_amp_threshold = 0.;

  std::shared_ptr<CContinuumModelSolution> m_fitContinuum;
  std::shared_ptr<fitMaxValues> m_fitContinuumMaxValues;

  // m_fitContinuum_option==2 for now
  Float64 m_fitContinuum_tplFitAlpha = 0.;

  void setFitContinuum_tplAmplitude(Float64 tplAmp, Float64 tplAmpErr,
                                    const TFloat64List &polyCoeffs);
};

} // namespace NSEpic
#endif
