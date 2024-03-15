#include "RedshiftLibrary/linemodel/continuummanager.h"
#include "RedshiftLibrary/linemodel/spectrummodel.h"

#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"

#include "RedshiftLibrary/operator/templatefitting.h"
#include "RedshiftLibrary/statistics/priorhelper.h"
using namespace NSEpic;
using namespace std;

CContinuumManager::CContinuumManager(const CSpcModelVectorPtr &models,
                                     std::shared_ptr<CTplModelSolution> tfv,
                                     std::shared_ptr<Int32> curObs)
    : m_tplCatalog(Context.GetTemplateCatalog()),
      m_tplCategory(Context.GetCurrentCategory()), m_models(models),
      m_fitContinuum(tfv), m_curObs(curObs) {

  // NB: fitContinuum_option: this is the initialization (default value),
  // eventually overriden in SetFitContinuum_FitStore() when a fitStore gets
  // available
  m_fitContinuum_option =
      0; // 0=interactive fitting, 1=use precomputed fit store, 2=use fixed
  // values (typical use for second pass recompute)

  CAutoScope autoscope(Context.m_ScopeStack, "lineModel");

  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();
  setContinuumComponent(ps->GetScoped<std::string>("continuumComponent"));

  if (m_ContinuumComponent != "noContinuum") {
    m_opt_fitcontinuum_neg_threshold =
        ps->GetScoped<Float64>("continuumFit.negativeThreshold");
    m_opt_fitcontinuum_null_amp_threshold =
        ps->GetScoped<Float64>("continuumFit.nullThreshold");
  }
}

std::shared_ptr<CPriorHelper> CContinuumManager::SetFitContinuum_PriorHelper() {
  CAutoScope a = CAutoScope(Context.m_ScopeStack, "lineModel");
  CAutoScope b = CAutoScope(Context.m_ScopeStack, "continuumFit");
  CAutoScope c = CAutoScope(Context.m_ScopeStack, "priors");
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();

  m_fitContinuum_priorhelper = std::make_shared<CPriorHelper>();
  m_fitContinuum_priorhelper->Init(ps->GetScoped<std::string>("catalogDirPath"),
                                   0);
  m_fitContinuum_priorhelper->SetBetaA(ps->GetScoped<Float64>("betaA"));
  m_fitContinuum_priorhelper->SetBetaTE(ps->GetScoped<Float64>("betaTE"));
  m_fitContinuum_priorhelper->SetBetaZ(ps->GetScoped<Float64>("betaZ"));
  return m_fitContinuum_priorhelper;
}

/**
 * \brief Generates a continuum from the fitting with a set of templates :
 * uses the templatefitting operator
 * TODO: LoadFitContinuum should be limited to reading continuum values from
 * the variable class, especially that we want that continuum fitting results
 * are saved in tplfitStore container outside CElementList and these stores
 * will be injected in the class whenever required !
 * TODO: study this possibility before doing the change
 */
void CContinuumManager::LoadFitContinuum(Int32 icontinuum, Int32 autoSelect,
                                         Float64 redshift) {
  Log.LogDebug(Formatter() << "Elementlist, m_fitContinuum_option="
                           << m_fitContinuum_option);
  /*if (m_observeGridContinuumFlux.empty())
    THROWG(INTERNAL_ERROR,
           "Cannot loadfitcontinuum without precomputedGridTplFlux");
  */
  if (m_fitContinuum_option ==
      1) { // using precomputed fit store, i.e., fitValues
    CTplModelSolution fitValues =
        m_fitContinuum_tplfitStore->GetFitValues(redshift, icontinuum);
    if (fitValues.tplName.empty()) {
      THROWG(INTERNAL_ERROR, "Empty template name");
    }

    *m_fitContinuum = fitValues;

    m_fitContinuumMaxValues = m_fitContinuum_tplfitStore->getFitMaxValues();

  } else if (m_fitContinuum_option == 2) {
    // values unmodified nothing to do
  } else {
    THROWG(INTERNAL_ERROR, "Cannot parse fitContinuum_option");
  }

  // check that continuum is well loaded
  if (m_fitContinuum->tplName.empty())
    THROWG(INTERNAL_ERROR, Formatter()
                               << "Failed to load-fit continuum for cfitopt="
                               << m_fitContinuum_option);

  // Retrieve the best template, otherwise Getter throws an error
  std::shared_ptr<const CTemplate> tpl =
      m_tplCatalog->GetTemplateByName({m_tplCategory}, m_fitContinuum->tplName);

  if (autoSelect) {
    // Float64 contsnr = getFitContinuum_snr();
    /*Float64 contsnr = m_fitContinuum->tplSNRMax;
    m_fitContinuum->tplAlpha = 1.0;
    if(contsnr>50.)
    {
        m_fitContinuum->tplAlpha=0.0;
    }*/
    m_fitContinuum_tplFitAlpha = 0.0;
    if (m_fitContinuumMaxValues->fitAmplitudeSigmaMAX <
        m_opt_fitcontinuum_neg_threshold)
      m_fitContinuum_tplFitAlpha = 1.0; // switch to spectrum continuum
  }
  //  for (; *m_curObs < m_models->size(); (*m_curObs)++) {
  getModel().ApplyContinuumOnGrid(tpl, m_fitContinuum->tplRedshift);
  setFitContinuum_tplAmplitude(m_fitContinuum->tplAmplitude,
                               m_fitContinuum->tplAmplitudeError,
                               m_fitContinuum->pCoeffs);

  //}
  //*m_curObs = 0;

  Log.LogDebug(Formatter() << "    model : LoadFitContinuum, loaded: "
                           << m_fitContinuum->tplName);
  Log.LogDebug(Formatter() << "    model : LoadFitContinuum, loaded with A="
                           << m_fitContinuum->tplAmplitude << ", with A_error="
                           << m_fitContinuum->tplAmplitudeError);
  Log.LogDebug(
      Formatter() << "    model : LoadFitContinuum, loaded with DustCoeff="
                  << m_fitContinuum->tplEbmvCoeff
                  << ", with MeiksinIdx=" << m_fitContinuum->tplMeiksinIdx);
  Log.LogDebug(Formatter() << "    model : LoadFitContinuum, loaded with dtm="
                           << m_fitContinuum->tplDtM
                           << ", with mtm=" << m_fitContinuum->tplMtM
                           << "with logprior=" << m_fitContinuum->tplLogPrior);
  Log.LogDebug(Formatter() << "    model : LoadFitContinuum, loaded with snr="
                           << m_fitContinuum->tplSNR);
}

void CContinuumManager::setFitContinuum_tplAmplitude(
    Float64 tplAmp, Float64 tplAmpErr, const TFloat64List &polyCoeffs) {

  Float64 alpha =
      m_fitContinuum_tplFitAlpha; // alpha blend = 1: only
                                  // m_inputSpc->GetContinuumFluxAxis(),
                                  // alpha=0: only tplfit

  m_fitContinuum->tplAmplitude = tplAmp;
  m_fitContinuum->tplAmplitudeError = tplAmpErr;
  m_fitContinuum->pCoeffs = polyCoeffs;
  for (; *m_curObs < m_models->size(); (*m_curObs)++) {
    getModel().setContinuumFromTplFit(alpha, tplAmp, polyCoeffs);
  }
}

Int32 CContinuumManager::SetFitContinuum_FitStore(
    const std::shared_ptr<const CTemplatesFitStore> &fitStore) {
  if (fitStore) {
    m_fitContinuum_option = 1; // enable use of the fit store
    Log.LogDetail("Elementlist: enabling fitContinuum store.");
  }
  m_fitContinuum_tplfitStore = fitStore;
  return 1;
}

const std::shared_ptr<const CTemplatesFitStore> &
CContinuumManager::GetFitContinuum_FitStore() const {
  return m_fitContinuum_tplfitStore;
}

void CContinuumManager::SetFitContinuum_Option(Int32 opt) {
  m_fitContinuum_option = opt;
}

void CContinuumManager::SetFitContinuum_SNRMax(Float64 snr_max) {
  m_fitContinuumMaxValues->tplFitSNRMax = snr_max;
}

Int32 CContinuumManager::GetFitContinuum_Option() const {
  return m_fitContinuum_option;
}

CTplModelSolution CContinuumManager::GetContinuumModelSolutionCopy() const {
  CTplModelSolution continuumModelSolution = *m_fitContinuum;

  return continuumModelSolution;
}

Float64 CContinuumManager::getContinuumScaleMargCorrection() const {
  Float64 corr = 0.0;

  // scale marg for continuum
  if (isContinuumComponentTplfitxx()) // the support has to be already
                                      // computed when LoadFitContinuum() is
                                      // called
    corr += log(m_fitContinuum->tplMtM);

  return corr;
}

std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
CContinuumManager::getIsmCorrectionFromTpl() {
  return m_tplCatalog->GetTemplate(m_tplCategory, 0)->m_ismCorrectionCalzetti;
}

void CContinuumManager::logParameters() {
  Log.LogInfo(Formatter() << "ContinuumComponent=" << m_ContinuumComponent);

  Log.LogInfo(Formatter() << "fitContinuum_option=" << m_fitContinuum_option);
  Log.LogInfo(Formatter() << "fitContinuum_tplName="
                          << m_fitContinuum->tplName);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitAmplitude="
                          << m_fitContinuum->tplAmplitude);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitAmplitudeError="
                          << m_fitContinuum->tplAmplitudeError);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitAmplitudeSigmaMAX="
                          << m_fitContinuumMaxValues->fitAmplitudeSigmaMAX);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitMerit="
                          << m_fitContinuum->tplMerit);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitMerit_phot="
                          << m_fitContinuum->tplMeritPhot);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitEbmvCoeff="
                          << m_fitContinuum->tplEbmvCoeff);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitMeiksinIdx="
                          << m_fitContinuum->tplMeiksinIdx);
  Log.LogInfo(
      Formatter()
      << "fitContinuum_tplFitRedshift="
      << m_fitContinuum->tplRedshift); // only used with
                                       // m_fitContinuum_option==2 for now
  Log.LogInfo(Formatter() << "fitContinuum_tplFitDtM="
                          << m_fitContinuum->tplDtM);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitMtM="
                          << m_fitContinuum->tplMtM);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitLogprior="
                          << m_fitContinuum->tplLogPrior);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitSNRMax="
                          << m_fitContinuumMaxValues->tplFitSNRMax);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitAlpha="
                          << m_fitContinuum_tplFitAlpha);
}

bool CContinuumManager::isContFittedToNull() {
  return m_fitContinuum->tplAmplitude < m_opt_fitcontinuum_null_amp_threshold *
                                            m_fitContinuum->tplAmplitudeError;
}

void CContinuumManager::setContinuumComponent(std::string component) {
  m_ContinuumComponent = std::move(component);
  for (*m_curObs = 0; *m_curObs < m_models->size(); (*m_curObs)++) {
    getModel().setContinuumComponent(m_ContinuumComponent);
  }
  *m_fitContinuum = {};

  if (m_ContinuumComponent == "noContinuum") {
    m_fitContinuum->tplName = "noContinuum"; // to keep track in resultstore
    // the continuum is set to zero and the observed spectrum is the spectrum
    // without continuum
  }
  if (m_ContinuumComponent == "fromSpectrum") {
    m_fitContinuum->tplName = "fromSpectrum"; // to keep track in resultstore
    // the continuum is set to the spectrum continuum and the observed
    // spectrum is the raw spectrum
  }
}

void CContinuumManager::reinterpolateContinuum(const Float64 redshift) {
  for (*m_curObs = 0; *m_curObs < m_models->size(); (*m_curObs)++) {

    std::shared_ptr<const CTemplate> tpl = m_tplCatalog->GetTemplateByName(
        {m_tplCategory}, m_fitContinuum->tplName);
    getModel().ApplyContinuumOnGrid(tpl, redshift);
  }
}

void CContinuumManager::reinterpolateContinuumResetAmp() {
  reinterpolateContinuum(m_fitContinuum->tplRedshift);
  m_fitContinuum->tplAmplitude = 1.0;
  m_fitContinuum->tplAmplitudeError = 1.0;
  TFloat64List polyCoeffs_unused;
  setFitContinuum_tplAmplitude(m_fitContinuum->tplAmplitude,
                               m_fitContinuum->tplAmplitudeError,
                               polyCoeffs_unused);
}

void CContinuumManager::setFitContinuumFromFittedAmps(
    TFloat64List &ampsfitted, TInt32List &validEltsIdx) {
  m_fitContinuum->tplAmplitude = ampsfitted[validEltsIdx.size()];
  TFloat64List polyCoeffs;
  for (Int32 kpoly = validEltsIdx.size() + 1; kpoly < ampsfitted.size();
       kpoly++) {
    polyCoeffs.push_back(ampsfitted[kpoly]);
  }
  setFitContinuum_tplAmplitude(m_fitContinuum->tplAmplitude,
                               m_fitContinuum->tplAmplitudeError, polyCoeffs);
}
