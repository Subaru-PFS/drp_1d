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
  enum class EFitType {
    interactiveFitting = 0,
    precomputedFitStore = 1,
    fixedValues = 2
  };

  CContinuumManager(const CSpcModelVectorPtr &models,
                    std::shared_ptr<CContinuumModelSolution>,
                    const CSpectraGlobalIndex &spcGlobIndex);

  const CSpectrumModel &getModel() const {
    return m_models->getSpectrumModel();
  }
  CSpectrumModel &getModel() { return m_models->getSpectrumModel(); }

  void SetFitContinuum_FitStore(
      const std::shared_ptr<const CContinuumFitStore> &fitStore);
  void SetFitContinuum_Option(EFitType opt);
  EFitType GetFitContinuum_Option() const;

  const std::shared_ptr<const CContinuumFitStore> &
  GetFitContinuum_FitStore() const;
  std::shared_ptr<CPriorHelper> SetFitContinuum_PriorHelper();
  void LoadFitContinuum(Int32 icontinuum, Float64 redshift);

  void SetFitContinuum_FitValues(const CContinuumModelSolution &cms) {
    *m_fitContinuum = cms;
  }

  Float64 getContinuumScaleMargCorrection() const;
  bool isContinuumComponentTplFitXXX() const {
    return m_ContinuumComponent.isTplFitXXX();
  }
  bool isContinuumComponentPowerLawXXX() const {
    return m_ContinuumComponent.isPowerLawXXX();
  }
  bool isContinuumComponentFitter() const {
    return m_ContinuumComponent.isContinuumFit();
  }

  bool isContinuumComponentNoContinuum() const {
    return m_ContinuumComponent.isNoContinuum();
  }

  Float64 getFitContinuum_snr() const;
  CContinuumModelSolution GetContinuumModelSolutionCopy() const;
  std::shared_ptr<const CContinuumModelSolution>
  GetContinuumModelSolution() const {
    return m_fitContinuum;
  }
  void setContinuumComponent(TContinuumComponent component);
  const TContinuumComponent &getContinuumComponent() const {
    return m_ContinuumComponent;
  };

  void logParameters();
  std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
  getIsmCorrectionFromTpl();
  void reinterpolateContinuum(Float64 redshift);
  void reinterpolateContinuumResetAmp();

  bool isContFittedToNull();
  Int32 getFittedMeiksinIndex() { return m_fitContinuum->meiksinIdx; }
  Float64 getFitSum() {
    if (!isContinuumComponentTplFitXXX())
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

  TContinuumComponent m_ContinuumComponent;
  EFitType m_fitContinuum_option;
  Float64 m_opt_fitcontinuum_neg_threshold = -INFINITY;
  Float64 m_opt_fitcontinuum_null_amp_threshold = 0.;

  std::shared_ptr<CContinuumModelSolution> m_fitContinuum;

  Float64 m_fitContinuum_tplFitAlpha = 0.;

  void setFitContinuum_tplAmplitude(Float64 tplAmp, Float64 tplAmpErr,
                                    const TFloat64List &polyCoeffs);
};

} // namespace NSEpic
#endif
