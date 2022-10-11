#include "RedshiftLibrary/linemodel/continuummanager.h"
#include "RedshiftLibrary/linemodel/spectrummodel.h"
#include "RedshiftLibrary/linemodel/templatesfitstore.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"

#include "RedshiftLibrary/operator/templatefitting.h"
#include "RedshiftLibrary/statistics/priorhelper.h"
using namespace NSEpic;
using namespace std;

CContinuumManager::CContinuumManager(
    const std::shared_ptr<const CSpectrum> &spc,
    const std::shared_ptr<CSpectrumModel> &model,
    const TLambdaRange &lambdaRange)
    : m_tplCatalog(Context.GetTemplateCatalog()),
      m_tplCategoryList({Context.GetCurrentCategory()}), m_model(model) {
  m_templateFittingOperator =
      std::make_shared<COperatorTemplateFitting>(*(spc), lambdaRange);

  // NB: fitContinuum_option: this is the initialization (default value),
  // eventually overriden in SetFitContinuum_FitStore() when a fitStore gets
  // available
  m_fitContinuum_option =
      0; // 0=interactive fitting, 1=use precomputed fit store, 2=use fixed
  // values (typical use for second pass recompute)

  CAutoScope autoscope(Context.m_ScopeStack, "linemodel");

  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();
  setContinuumComponent(ps->GetScoped<std::string>("continuumcomponent"));

  if (m_ContinuumComponent != "nocontinuum") {
    m_opt_fitcontinuum_neg_threshold =
        ps->GetScoped<Float64>("continuumfit.negativethreshold");
    m_opt_fitcontinuum_null_amp_threshold =
        ps->GetScoped<Float64>("continuumfit.nullthreshold");
  }
}

/**
 * Apply the template continuum by interpolating the grid as define in Init
 * Continuum
 */
Int32 CContinuumManager::ApplyContinuumOnGrid(
    const std::shared_ptr<const CTemplate> &tpl, Float64 zcontinuum) {
  m_fitContinuum_tplName = tpl->GetName();
  Int32 n = tpl->GetSampleCount();

  Int32 idxDust = -1;
  if (m_fitContinuum_tplFitEbmvCoeff > 0.) {
    if (tpl->CalzettiInitFailed()) {
      THROWG(INTERNAL_ERROR, "  no calzetti calib. file in template");
    }
    idxDust = tpl->m_ismCorrectionCalzetti->GetEbmvIndex(
        m_fitContinuum_tplFitEbmvCoeff);
  }
  const CSpectrumSpectralAxis &tplSpectralAxis = tpl->GetSpectralAxis();
  TFloat64Range range(tplSpectralAxis[0], tplSpectralAxis[n - 1]);

  std::string inter_opt = "spline";
  tpl->setRebinInterpMethod(inter_opt);
  Float64 overlapThreshold = 1., amplitude = 1.;
  std::shared_ptr<CModelSpectrumResult> spcmodel =
      m_templateFittingOperator->ComputeSpectrumModel(
          tpl, zcontinuum, m_fitContinuum_tplFitEbmvCoeff,
          m_fitContinuum_tplFitMeiksinIdx, amplitude, overlapThreshold);
  if (spcmodel == nullptr)
    THROWG(INTERNAL_ERROR, "Couldnt compute spectrum model");

  // m_observeGridContinuumFlux should be a CSpectrumFluxAxis not
  // AxisSampleList
  m_observeGridContinuumFlux = std::move((*spcmodel).ModelFlux);

  return 0;
}

std::shared_ptr<CPriorHelper> CContinuumManager::SetFitContinuum_PriorHelper() {
  CAutoScope a = CAutoScope(Context.m_ScopeStack, "linemodel");
  CAutoScope b = CAutoScope(Context.m_ScopeStack, "continuumfit");
  CAutoScope c = CAutoScope(Context.m_ScopeStack, "priors");
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();

  m_fitContinuum_priorhelper = std::make_shared<CPriorHelper>();
  m_fitContinuum_priorhelper->Init(
      ps->GetScoped<std::string>("catalog_dirpath"), 0);
  m_fitContinuum_priorhelper->SetBetaA(ps->GetScoped<Float64>("betaA"));
  m_fitContinuum_priorhelper->SetBetaTE(ps->GetScoped<Float64>("betaTE"));
  m_fitContinuum_priorhelper->SetBetaZ(ps->GetScoped<Float64>("betaZ"));
  return m_fitContinuum_priorhelper;
}

void CContinuumManager::SolveContinuum(
    const std::shared_ptr<const CTemplate> &tpl, const TFloat64List &redshifts,
    Float64 overlapThreshold, std::vector<CMask> maskList,
    std::string opt_interp, Int32 opt_extinction, Int32 opt_dustFit,
    Float64 &merit, Float64 &fitAmplitude, Float64 &fitAmplitudeError,
    Float64 &fitAmplitudeSigma, Float64 &FitEbmvCoeff, Int32 &fitMeiksinIdx,
    Float64 &fitDtM, Float64 &fitMtM, Float64 &fitLogprior) {
  CPriorHelper::TPriorZEList zePriorData;

  m_fitContinuum_priorhelper->GetTplPriorData(tpl->GetName(), redshifts,
                                              zePriorData);
  bool keepigmism = false;
  if (FitEbmvCoeff + fitMeiksinIdx != -2) {
    keepigmism = true;
  }

  // Compute merit function
  // Log.LogInfo("Solving continuum for %s at z=%.4e", tpl.GetName().c_str(),
  // redshifts[0]); CRef<CChisquareResult>  chisquareResult =
  // (CChisquareResult*)chiSquare.ExportChi2versusAZ( _spc, _tpl, lambdaRange,
  // redshifts, overlapThreshold );
  m_templateFittingOperator->SetRedshifts(redshifts);
  tpl->setRebinInterpMethod(opt_interp);
  auto templateFittingResult =
      std::dynamic_pointer_cast<CTemplateFittingResult>(
          m_templateFittingOperator->Compute(
              tpl, overlapThreshold, maskList, opt_interp, opt_extinction,
              opt_dustFit, zePriorData, keepigmism, FitEbmvCoeff,
              fitMeiksinIdx));

  if (!templateFittingResult)
    THROWG(INTERNAL_ERROR, "Failed to compute chi square value");

  // Store results
  merit = templateFittingResult->ChiSquare[0];
  fitAmplitude = templateFittingResult->FitAmplitude[0];
  fitAmplitudeError = templateFittingResult->FitAmplitudeError[0];
  fitAmplitudeSigma = templateFittingResult->FitAmplitudeSigma[0];
  FitEbmvCoeff = templateFittingResult->FitEbmvCoeff[0];
  fitMeiksinIdx = templateFittingResult->FitMeiksinIdx[0];
  fitDtM = templateFittingResult->FitDtM[0];
  fitMtM = templateFittingResult->FitMtM[0];
  fitLogprior = templateFittingResult->LogPrior[0];
}

const std::string &CContinuumManager::getFitContinuum_tplName() const {
  return m_fitContinuum_tplName;
}

Float64 CContinuumManager::getFitContinuum_tplAmplitude() const {
  return m_fitContinuum_tplFitAmplitude;
}

Float64 CContinuumManager::getFitContinuum_tplAmplitudeError() const {
  return m_fitContinuum_tplFitAmplitudeError;
}

// This SNR estimate maybe needs to use observed spectrum with lines removed ?
Float64 CContinuumManager::getFitContinuum_snr() const {
  Float64 snr = -1.; // shoudnt be NAN here?
  if (m_fitContinuum_tplFitMtM > 0.) {
    snr = m_fitContinuum_tplFitDtM / std::sqrt(m_fitContinuum_tplFitMtM);
  }
  return snr;
}

Float64 CContinuumManager::getFitContinuum_tplMerit() const {
  return m_fitContinuum_tplFitMerit;
}

Float64 CContinuumManager::getFitContinuum_tplMeritPhot() const {
  return m_fitContinuum_tplFitMerit_phot;
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
  Log.LogDebug("Elementlist, m_fitContinuum_option=%d", m_fitContinuum_option);
  if (m_observeGridContinuumFlux.empty())
    THROWG(INTERNAL_ERROR,
           "Cannot loadfitcontinuum without precomputedGridTplFlux");

  if (m_fitContinuum_option ==
      1) { // using precomputed fit store, i.e., fitValues
    CTemplatesFitStore::TemplateFitValues fitValues =
        m_fitContinuum_tplfitStore->GetFitValues(redshift, icontinuum);
    if (fitValues.tplName.empty()) {
      THROWG(INTERNAL_ERROR, "Empty template name");
    }

    m_fitContinuum_tplName = fitValues.tplName;
    m_fitContinuum_tplFitAmplitude = fitValues.fitAmplitude;
    m_fitContinuum_tplFitAmplitudeError = fitValues.fitAmplitudeError;
    m_fitContinuum_tplFitMerit = fitValues.merit;
    m_fitContinuum_tplFitMerit_phot = fitValues.chiSquare_phot;
    m_fitContinuum_tplFitEbmvCoeff = fitValues.ismEbmvCoeff;
    m_fitContinuum_tplFitMeiksinIdx = fitValues.igmMeiksinIdx;
    m_fitContinuum_tplFitRedshift = redshift;
    m_fitContinuum_tplFitDtM = fitValues.fitDtM;
    m_fitContinuum_tplFitMtM = fitValues.fitMtM;
    m_fitContinuum_tplFitLogprior = fitValues.logprior;
    m_fitContinuum_tplFitPolyCoeffs = {};

    m_fitContinuum_tplFitAmplitudeSigmaMAX =
        m_fitContinuum_tplfitStore->m_fitContinuum_fitAmplitudeSigmaMAX;
    m_fitContinuum_tplFitSNRMax =
        m_fitContinuum_tplfitStore->m_fitContinuum_tplFitSNRMax;
    m_opt_fitcontinuum_maxCount =
        m_fitContinuum_tplfitStore->m_opt_fitcontinuum_maxCount;

  } else if (m_fitContinuum_option == 2) {
    // values unmodified nothing to do
  } else {
    THROWG(INTERNAL_ERROR, "Cannot parse fitContinuum_option");
  }

  // check that continuum is well loaded
  if (m_fitContinuum_tplName.empty())
    THROWG(INTERNAL_ERROR, Formatter()
                               << "Failed to load-fit continuum for cfitopt="
                               << m_fitContinuum_option);

  // Retrieve the best template, otherwise Getter throws an error
  std::shared_ptr<const CTemplate> tpl = m_tplCatalog->GetTemplateByName(
      m_tplCategoryList, m_fitContinuum_tplName);

  if (autoSelect) {
    // Float64 contsnr = getFitContinuum_snr();
    /*Float64 contsnr = m_fitContinuum_tplFitSNRMax;
    m_fitContinuum_tplFitAlpha = 1.0;
    if(contsnr>50.)
    {
        m_fitContinuum_tplFitAlpha=0.0;
    }*/
    m_fitContinuum_tplFitAlpha = 0.0;
    if (m_fitContinuum_tplFitAmplitudeSigmaMAX <
        m_opt_fitcontinuum_neg_threshold)
      m_fitContinuum_tplFitAlpha = 1.0; // switch to spectrum continuum
  }

  ApplyContinuumOnGrid(tpl, m_fitContinuum_tplFitRedshift);

  setFitContinuum_tplAmplitude(m_fitContinuum_tplFitAmplitude,
                               m_fitContinuum_tplFitAmplitudeError,
                               m_fitContinuum_tplFitPolyCoeffs);

  Log.LogDebug("    model : LoadFitContinuum, loaded: %s",
               m_fitContinuum_tplName.c_str());
  Log.LogDebug(
      "    model : LoadFitContinuum, loaded with A=%e, with A_error=%e",
      m_fitContinuum_tplFitAmplitude, m_fitContinuum_tplFitAmplitudeError);
  Log.LogDebug("    model : LoadFitContinuum, loaded with DustCoeff=%e, with "
               "MeiksinIdx=%d",
               m_fitContinuum_tplFitEbmvCoeff, m_fitContinuum_tplFitMeiksinIdx);
  Log.LogDebug("    model : LoadFitContinuum, loaded with dtm=%e, with mtm=%e, "
               "with logprior=%e",
               m_fitContinuum_tplFitDtM, m_fitContinuum_tplFitMtM,
               m_fitContinuum_tplFitLogprior);
  Float64 tplfitsnr = NAN;
  if (m_fitContinuum_tplFitMtM > 0.0) {
    tplfitsnr = m_fitContinuum_tplFitDtM / std::sqrt(m_fitContinuum_tplFitMtM);
  }
  Log.LogDebug("    model : LoadFitContinuum, loaded with snr=%e", tplfitsnr);
}

void CContinuumManager::setFitContinuum_tplAmplitude(
    Float64 tplAmp, Float64 tplAmpErr, const TFloat64List &polyCoeffs) {

  Float64 alpha =
      m_fitContinuum_tplFitAlpha; // alpha blend = 1: only
                                  // m_inputSpc->GetContinuumFluxAxis(),
                                  // alpha=0: only tplfit

  m_fitContinuum_tplFitAmplitude = tplAmp;
  m_fitContinuum_tplFitAmplitudeError = tplAmpErr;
  m_fitContinuum_tplFitPolyCoeffs = polyCoeffs;
  m_model->setContinuumFromTplFit(alpha, tplAmp, polyCoeffs,
                                  m_observeGridContinuumFlux);
}

Float64 CContinuumManager::getFitContinuum_tplIsmEbmvCoeff() const {
  return m_fitContinuum_tplFitEbmvCoeff;
}

Float64 CContinuumManager::getFitContinuum_tplIgmMeiksinIdx() const {
  return m_fitContinuum_tplFitMeiksinIdx;
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
  m_fitContinuum_tplFitSNRMax = snr_max;
}

Int32 CContinuumManager::GetFitContinuum_Option() const {
  return m_fitContinuum_option;
}

void CContinuumManager::SetFitContinuum_FitValues(
    std::string tplfit_name, Float64 tplfit_amp, Float64 tplfit_amperr,
    Float64 tplfit_chi2, Float64 tplfit_chi2_phot, Float64 tplfit_ebmv,
    Int32 tplfit_meiksinidx, Float64 tplfit_continuumredshift,
    Float64 tplfit_dtm, Float64 tplfit_mtm, Float64 tplfit_logprior,
    const TFloat64List &polyCoeffs) {
  m_fitContinuum_tplName = tplfit_name;
  m_fitContinuum_tplFitAmplitude = tplfit_amp;
  m_fitContinuum_tplFitAmplitudeError = tplfit_amperr;
  m_fitContinuum_tplFitMerit = tplfit_chi2;
  m_fitContinuum_tplFitMerit_phot = tplfit_chi2_phot;
  m_fitContinuum_tplFitEbmvCoeff = tplfit_ebmv;
  m_fitContinuum_tplFitMeiksinIdx = tplfit_meiksinidx;
  m_fitContinuum_tplFitRedshift = tplfit_continuumredshift;

  m_fitContinuum_tplFitDtM = tplfit_dtm;
  m_fitContinuum_tplFitMtM = tplfit_mtm;
  m_fitContinuum_tplFitLogprior = tplfit_logprior;
  m_fitContinuum_tplFitPolyCoeffs = polyCoeffs;
}

CContinuumModelSolution CContinuumManager::GetContinuumModelSolution() const {
  CContinuumModelSolution continuumModelSolution;

  continuumModelSolution.tplName = m_fitContinuum_tplName;
  continuumModelSolution.tplEbmvCoeff = m_fitContinuum_tplFitEbmvCoeff;
  continuumModelSolution.tplMeiksinIdx = m_fitContinuum_tplFitMeiksinIdx;
  continuumModelSolution.tplAmplitude = m_fitContinuum_tplFitAmplitude;
  continuumModelSolution.tplAmplitudeError =
      m_fitContinuum_tplFitAmplitudeError;
  continuumModelSolution.tplMerit = m_fitContinuum_tplFitMerit;
  continuumModelSolution.tplMeritPhot = m_fitContinuum_tplFitMerit_phot;
  continuumModelSolution.tplDtm = m_fitContinuum_tplFitDtM;
  continuumModelSolution.tplMtm = m_fitContinuum_tplFitMtM;
  continuumModelSolution.tplLogPrior = m_fitContinuum_tplFitLogprior;
  continuumModelSolution.tplRedshift = m_fitContinuum_tplFitRedshift;
  continuumModelSolution.pCoeffs = m_fitContinuum_tplFitPolyCoeffs;

  return continuumModelSolution;
}

Float64 CContinuumManager::getContinuumScaleMargCorrection() const {
  Float64 corr = 0.0;

  // scale marg for continuum
  if (isContinuumComponentTplfitxx()) // the support has to be already
                                      // computed when LoadFitContinuum() is
                                      // called
    corr += log(m_fitContinuum_tplFitMtM);

  return corr;
}

std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
CContinuumManager::getIsmCorrectionFromTpl() {
  return m_tplCatalog->GetTemplate(m_tplCategoryList[0], 0)
      ->m_ismCorrectionCalzetti;
}

void CContinuumManager::logParameters() {
  Log.LogInfo(Formatter() << "ContinuumComponent=" << m_ContinuumComponent);

  Log.LogInfo(Formatter() << "fitContinuum_option=" << m_fitContinuum_option);
  Log.LogInfo(Formatter() << "fitContinuum_tplName=" << m_fitContinuum_tplName);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitAmplitude="
                          << m_fitContinuum_tplFitAmplitude);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitAmplitudeError="
                          << m_fitContinuum_tplFitAmplitudeError);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitAmplitudeSigmaMAX="
                          << m_fitContinuum_tplFitAmplitudeSigmaMAX);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitMerit="
                          << m_fitContinuum_tplFitMerit);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitMerit_phot="
                          << m_fitContinuum_tplFitMerit_phot);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitEbmvCoeff="
                          << m_fitContinuum_tplFitEbmvCoeff);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitMeiksinIdx="
                          << m_fitContinuum_tplFitMeiksinIdx);
  Log.LogInfo(
      Formatter()
      << "fitContinuum_tplFitRedshift="
      << m_fitContinuum_tplFitRedshift); // only used with
                                         // m_fitContinuum_option==2 for now
  Log.LogInfo(Formatter() << "fitContinuum_tplFitDtM="
                          << m_fitContinuum_tplFitDtM);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitMtM="
                          << m_fitContinuum_tplFitMtM);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitLogprior="
                          << m_fitContinuum_tplFitLogprior);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitSNRMax="
                          << m_fitContinuum_tplFitSNRMax);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitAlpha="
                          << m_fitContinuum_tplFitAlpha);
}

void CContinuumManager::initObserveGridContinuumFlux(Int32 size) {
  m_observeGridContinuumFlux.resize(size);
}

bool CContinuumManager::isContFittedToNull() {
  return m_fitContinuum_tplFitAmplitude >
         m_opt_fitcontinuum_null_amp_threshold *
             m_fitContinuum_tplFitAmplitudeError;
}

void CContinuumManager::setContinuumComponent(std::string component) {
  m_ContinuumComponent = std::move(component);
  m_model->setContinuumComponent(m_ContinuumComponent);

  if (m_ContinuumComponent == "nocontinuum") {
    m_fitContinuum_tplName = "nocontinuum"; // to keep track in resultstore
    // the continuum is set to zero and the observed spectrum is the spectrum
    // without continuum
  }
  if (m_ContinuumComponent == "fromspectrum") {
    m_fitContinuum_tplName = "fromspectrum"; // to keep track in resultstore
    // the continuum is set to the spectrum continuum and the observed
    // spectrum is the raw spectrum
  }
}

void CContinuumManager::reinterpolateContinuum(const Float64 redshift) {
  std::shared_ptr<const CTemplate> tpl = m_tplCatalog->GetTemplateByName(
      m_tplCategoryList, m_fitContinuum_tplName);
  ApplyContinuumOnGrid(tpl, redshift);
}

void CContinuumManager::reinterpolateContinuumResetAmp() {
  reinterpolateContinuum(m_fitContinuum_tplFitRedshift);
  m_fitContinuum_tplFitAmplitude = 1.0;
  m_fitContinuum_tplFitAmplitudeError = 1.0;
  TFloat64List polyCoeffs_unused;
  setFitContinuum_tplAmplitude(m_fitContinuum_tplFitAmplitude,
                               m_fitContinuum_tplFitAmplitudeError,
                               polyCoeffs_unused);
}

void CContinuumManager::setFitContinuumFromFittedAmps(
    TFloat64List &ampsfitted, TInt32List &validEltsIdx) {
  m_fitContinuum_tplFitAmplitude = ampsfitted[validEltsIdx.size()];
  TFloat64List polyCoeffs;
  for (Int32 kpoly = validEltsIdx.size() + 1; kpoly < ampsfitted.size();
       kpoly++) {
    polyCoeffs.push_back(ampsfitted[kpoly]);
  }
  setFitContinuum_tplAmplitude(m_fitContinuum_tplFitAmplitude,
                               m_fitContinuum_tplFitAmplitudeError, polyCoeffs);
}
