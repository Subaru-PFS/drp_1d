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

#include "RedshiftLibrary/linemodel/continuummanager.h"
#include "RedshiftLibrary/linemodel/spectrummodel.h"

#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"

#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/operator/powerlaw.h"
#include "RedshiftLibrary/operator/templatefitting.h"
#include "RedshiftLibrary/statistics/priorhelper.h"

using namespace NSEpic;
using namespace std;

CContinuumManager::CContinuumManager(
    const CSpcModelVectorPtr &models,
    std::shared_ptr<CContinuumModelSolution> continuumModelSolution,
    const CSpectraGlobalIndex &spcGlobIndex)
    : m_tplCatalog(Context.GetTemplateCatalog()),
      m_tplCategory(Context.GetCurrentCategory()), m_models(models),
      m_spectraIndex(spcGlobIndex), m_fitContinuum(continuumModelSolution) {

  // NB: fitContinuum_option: this is the initialization (default value),
  // eventually overriden in SetFitContinuum_FitStore() when a fitStore gets
  // available
  m_fitContinuum_option = EFitType::interactiveFitting;

  CAutoScope autoscope(Context.m_ScopeStack, "lineModel");

  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();
  setContinuumComponent(
      TContinuumComponent(ps->GetScoped<std::string>("continuumComponent")));

  if (!m_ContinuumComponent.isNoContinuum()) {
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
void CContinuumManager::LoadFitContinuum(Int32 icontinuum, Float64 redshift) {
  Log.LogDebug(Formatter() << "Elementlist, m_fitContinuum_option="
                           << static_cast<Int32>(m_fitContinuum_option));
  if (m_fitContinuum_option ==
      EFitType::precomputedFitStore) { // using precomputed fit store, i.e.,
                                       // fitValues
    CContinuumModelSolution fitValues =
        m_fitContinuum_tplfitStore->GetFitValues(redshift, icontinuum);
    if (fitValues.name.empty()) {
      THROWG(ErrorCode::INTERNAL_ERROR, "Empty template name");
    }

    *m_fitContinuum = fitValues;
  } else if (m_fitContinuum_option == EFitType::fixedValues) {
    // values unmodified nothing to do
  } else {
    THROWG(ErrorCode::INTERNAL_ERROR, "Cannot parse fitContinuum_option");
  }
  if (m_fitContinuum->name.empty())
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Failed to load-fit continuum for cfitopt="
                       << static_cast<Int32>(m_fitContinuum_option));
  // Retrieve the best template, otherwise Getter throws an error
  if (isContinuumComponentPowerLawXXX()) {
    getModel().ApplyContinuumPowerLawOnGrid(m_fitContinuum);
  } else {
    std::shared_ptr<const CTemplate> tpl =
        m_tplCatalog->GetTemplateByName({m_tplCategory}, m_fitContinuum->name);

    getModel().ApplyContinuumTplOnGrid(tpl, m_fitContinuum->redshift);

    setFitContinuum_tplAmplitude(m_fitContinuum->tplAmplitude,
                                 m_fitContinuum->tplAmplitudeError,
                                 m_fitContinuum->pCoeffs);

    Log.LogDebug(Formatter() << "    model : LoadFitContinuum, loaded: "
                             << m_fitContinuum->name);
    Log.LogDebug(Formatter()
                 << "    model : LoadFitContinuum, loaded with A="
                 << m_fitContinuum->tplAmplitude
                 << ", with A_error=" << m_fitContinuum->tplAmplitudeError);
    Log.LogDebug(Formatter()
                 << "    model : LoadFitContinuum, loaded with DustCoeff="
                 << m_fitContinuum->ebmvCoef
                 << ", with meiksinIdx=" << m_fitContinuum->meiksinIdx);
    Log.LogDebug(Formatter()
                 << "    model : LoadFitContinuum, loaded with dtm="
                 << m_fitContinuum->tplDtM
                 << ", with mtm=" << m_fitContinuum->tplMtM
                 << "with logprior=" << m_fitContinuum->tplLogPrior);
    Log.LogDebug(Formatter() << "    model : LoadFitContinuum, loaded with snr="
                             << m_fitContinuum->SNR);
  }
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
  getModel().setContinuumFromTplFit(alpha, tplAmp, polyCoeffs);
}

void CContinuumManager::SetFitContinuum_FitStore(
    const std::shared_ptr<const CContinuumFitStore> &fitStore) {
  if (fitStore) {
    m_fitContinuum_option =
        EFitType::precomputedFitStore; // enable use of the fit store
    Log.LogDetail("Elementlist: enabling fitContinuum store.");
  }
  m_fitContinuum_tplfitStore = fitStore;
}

const std::shared_ptr<const CContinuumFitStore> &
CContinuumManager::GetFitContinuum_FitStore() const {
  return m_fitContinuum_tplfitStore;
}

void CContinuumManager::SetFitContinuum_Option(EFitType opt) {
  m_fitContinuum_option = opt;
}

CContinuumManager::EFitType CContinuumManager::GetFitContinuum_Option() const {
  return m_fitContinuum_option;
}

CContinuumModelSolution
CContinuumManager::GetContinuumModelSolutionCopy() const {
  CContinuumModelSolution continuumModelSolution = *m_fitContinuum;

  return continuumModelSolution;
}

Float64 CContinuumManager::getContinuumScaleMargCorrection() const {
  Float64 corr = 0.0;

  // scale marg for continuum
  if (isContinuumComponentTplFitXXX()) // the support has to be already
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
  Log.LogInfo(Formatter() << "ContinuumComponent="
                          << std::string(m_ContinuumComponent));

  Log.LogInfo(Formatter() << "fitContinuum_option="
                          << static_cast<Int32>(m_fitContinuum_option));
  Log.LogInfo(Formatter() << "fitContinuum_tplName=" << m_fitContinuum->name);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitAmplitude="
                          << m_fitContinuum->tplAmplitude);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitAmplitudeError="
                          << m_fitContinuum->tplAmplitudeError);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitMerit="
                          << m_fitContinuum->merit);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitReducedChi2="
                          << m_fitContinuum->reducedChi2);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitMerit_phot="
                          << m_fitContinuum->tplMeritPhot);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitEbmvCoeff="
                          << m_fitContinuum->ebmvCoef);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitMeiksinIdx="
                          << m_fitContinuum->meiksinIdx);
  Log.LogInfo(
      Formatter()
      << "fitContinuum_tplFitRedshift="
      << m_fitContinuum->redshift); // only used with
                                    // m_fitContinuum_option==fixedValues
                                    // for now
  Log.LogInfo(Formatter() << "fitContinuum_tplFitDtM="
                          << m_fitContinuum->tplDtM);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitMtM="
                          << m_fitContinuum->tplMtM);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitLogprior="
                          << m_fitContinuum->tplLogPrior);
  Log.LogInfo(Formatter() << "fitContinuum_tplFitAlpha="
                          << m_fitContinuum_tplFitAlpha);
}

bool CContinuumManager::isContFittedToNull() {
  return m_fitContinuum->tplAmplitude < m_opt_fitcontinuum_null_amp_threshold *
                                            m_fitContinuum->tplAmplitudeError;
}

void CContinuumManager::setContinuumComponent(TContinuumComponent component) {
  m_ContinuumComponent = std::move(component);
  for (auto &spcIndex : m_spectraIndex) {
    getModel().setContinuumComponent(m_ContinuumComponent);
  }
  *m_fitContinuum = {};

  if (m_ContinuumComponent.isNoContinuum()) {
    m_fitContinuum->name = "noContinuum"; // to keep track in resultstore
    // the continuum is set to zero and the observed spectrum is the spectrum
    // without continuum
  }
  if (m_ContinuumComponent.isFromSpectrum()) {
    m_fitContinuum->name = "fromSpectrum"; // to keep track in resultstore
    // the continuum is set to the spectrum continuum and the observed
    // spectrum is the raw spectrum
  }
  if (m_ContinuumComponent.isPowerLawXXX()) {
    m_fitContinuum->name = "powerLaw";
  }
}

void CContinuumManager::reinterpolateContinuum(const Float64 redshift) {
  for (auto &spcIndex : m_spectraIndex) {
    std::shared_ptr<const CTemplate> tpl =
        m_tplCatalog->GetTemplateByName({m_tplCategory}, m_fitContinuum->name);
    getModel().ApplyContinuumTplOnGrid(tpl, redshift);
  }
}

void CContinuumManager::reinterpolateContinuumResetAmp() {
  reinterpolateContinuum(m_fitContinuum->redshift);
  m_fitContinuum->tplAmplitude = 1.0;
  m_fitContinuum->tplAmplitudeError = 1.0;
  TFloat64List polyCoeffs_unused;
  for (auto &spcIndex : m_spectraIndex) {
    setFitContinuum_tplAmplitude(m_fitContinuum->tplAmplitude,
                                 m_fitContinuum->tplAmplitudeError,
                                 polyCoeffs_unused);
  }
}

void CContinuumManager::setFitContinuumFromFittedAmps(
    TFloat64List &ampsfitted, TInt32List &validEltsIdx) {
  m_fitContinuum->tplAmplitude = ampsfitted[ssize(validEltsIdx)];
  TFloat64List polyCoeffs;
  for (Int32 kpoly = ssize(validEltsIdx) + 1; kpoly < ssize(ampsfitted);
       kpoly++) {
    polyCoeffs.push_back(ampsfitted[kpoly]);
  }
  setFitContinuum_tplAmplitude(m_fitContinuum->tplAmplitude,
                               m_fitContinuum->tplAmplitudeError, polyCoeffs);
}
