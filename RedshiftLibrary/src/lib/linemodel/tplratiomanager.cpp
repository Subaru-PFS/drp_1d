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

#include "RedshiftLibrary/linemodel/tplratiomanager.h"
#include "RedshiftLibrary/line/catalogsTplRatio.h"
#include "RedshiftLibrary/linemodel/continuummanager.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/linemodel/spectrummodel.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;

CTplratioManager::CTplratioManager(
    CLineModelElementList &elements, std::shared_ptr<CSpectrumModel> model,
    std::shared_ptr<const CSpectrum> inputSpc,
    std::shared_ptr<const TFloat64Range> lambdaRange,
    std::shared_ptr<CContinuumManager> continuumManager,
    const CLineCatalog::TLineVector &restLineList)
    : CLineRatioManager(elements, model, inputSpc, lambdaRange,
                        continuumManager, restLineList) {
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();

  CAutoScope autoscope(Context.m_ScopeStack, "linemodel");

  m_CatalogTplRatio = Context.GetTplRatioCatalog();
  initTplratioCatalogs(
      Context.GetParameterStore()->GetScoped<bool>("tplratio_ismfit"));
  SetTplratio_PriorHelper();

  m_opt_firstpass_forcedisableTplratioISMfit =
      !ps->GetScoped<bool>("firstpass.tplratio_ismfit");
  m_NSigmaSupport = ps->GetScoped<Float64>("nsigmasupport");
  m_opt_haprior = ps->GetScoped<Float64>("haprior");
}

void CTplratioManager::initMerit(Int32 ntplratio) {
  m_bestTplratioMerit.clear();
  m_bestTplratioMeritPrior.clear();
  m_bestTplratioMerit.resize(ntplratio, INFINITY);
  m_bestTplratioMeritPrior.resize(ntplratio, 0.0);
}

void CTplratioManager::SetTplratio_PriorHelper() {
  CAutoScope b = CAutoScope(Context.m_ScopeStack, "tplratio");
  CAutoScope c = CAutoScope(Context.m_ScopeStack, "priors");
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();

  m_tplratio_priorhelper = std::make_shared<CPriorHelper>();
  m_tplratio_priorhelper->Init(ps->GetScoped<std::string>("catalog_dirpath"),
                               1);
  m_tplratio_priorhelper->SetBetaA(ps->GetScoped<Float64>("betaA"));
  m_tplratio_priorhelper->SetBetaTE(ps->GetScoped<Float64>("betaTE"));
  m_tplratio_priorhelper->SetBetaZ(ps->GetScoped<Float64>("betaZ"));
}

Int32 CTplratioManager::getLineIndexInCatalog(
    Int32 iElts, Int32 idxLine,
    const CLineCatalog::TLineVector &catalog) const {
  Int32 lineIndex = undefIdx;

  // get index of line inside tplratio catalog
  const std::string &strID = m_Elements[iElts]->m_Lines[idxLine].GetStrID();
  lineIndex = std::find_if(catalog.begin(), catalog.end(),
                           [strID](const CLine &line) {
                             return line.GetStrID() == strID;
                           }) -
              catalog.begin();
  if (lineIndex >= catalog.size())
    lineIndex = undefIdx;

  return lineIndex;
}

Int32 CTplratioManager::prepareFit(Float64 redshift) {
  m_logPriorDataTplRatio.clear();
  m_savedIdxFitted = -1;
  Int32 ntplratio = m_CatalogTplRatio->GetCatalogsCount();
  initMerit(ntplratio);
  if (!m_tplratio_priorhelper->mInitFailed) {
    // prior initilization for tplratio EL only
    if (m_Elements.size() > 1)
      THROWG(INTERNAL_ERROR, "model: Unable to use tplratio line priors "
                             "with nElts>1 for now");
    // NB: this could be done if the EL element idx in searched (see later
    // in the itratio loop, UV Abs lines would be not affected by priors
    // then)

    for (Int32 itratio = 0; itratio < ntplratio; itratio++) {
      // prepare the lines prior data
      Int32 ebvfilter = m_CatalogTplRatio->GetIsmIndex(itratio);
      CPriorHelper::SPriorTZE logPriorData;
      std::string tplrationame = m_CatalogTplRatio->GetCatalogName(itratio);
      bool retGetPrior = m_tplratio_priorhelper->GetTZEPriorData(
          tplrationame, ebvfilter, redshift, logPriorData);
      if (retGetPrior == false)
        THROWG(INTERNAL_ERROR,
               "model: Failed to get prior for chi2 solvecontinuum.");
      else
        m_logPriorDataTplRatio.push_back(logPriorData);
    }
  }
  return ntplratio;
}

bool CTplratioManager::init(Float64 redshift, Int32 itratio) {
  if (m_forcedisableTplratioISMfit && itratio > 0 &&
      m_CatalogTplRatio->GetIsmIndex(itratio) > 0) {
    duplicateTplratioResult(itratio, m_bestTplratioMerit,
                            m_bestTplratioMeritPrior);
    return true;
  }
  setTplratioModel(itratio, false);
  // prepare the Lya width and asym coefficients if the asymfit profile
  // option is met INFO: tpl-shape are often ASYMFIXED in the tplratio
  // catalog files, for the lyaE profile, as of 2016-01-11 INFO:
  // tplratio can override the lyafitting, see m_opt_lya_forcefit
  setLyaProfile(redshift, m_CatalogTplRatio->GetCatalog(itratio).GetList());
  return false;
}

/**
 * @brief :copy the values for ebmv=ebmv_fixed (=0)
 *
 * @param idx
 */
void CTplratioManager::duplicateTplratioResult(
    Int32 idx, TFloat64List &bestTplratioMerit,
    TFloat64List &bestTplratioMeritPrior) {
  bestTplratioMerit[idx] = bestTplratioMerit[idx - 1];
  bestTplratioMeritPrior[idx] = bestTplratioMeritPrior[idx - 1];
  m_ChisquareTplratio[idx] = m_ChisquareTplratio[idx - 1];
  m_ScaleMargCorrTplratio[idx] = m_ScaleMargCorrTplratio[idx - 1];
  m_StrongELPresentTplratio[idx] = m_StrongELPresentTplratio[idx - 1];
  m_StrongHalphaELPresentTplratio[idx] =
      m_StrongHalphaELPresentTplratio[idx - 1];
  m_NLinesAboveSNRTplratio[idx] = m_NLinesAboveSNRTplratio[idx - 1];

  for (Int32 iElt = 0; iElt < m_Elements.size(); iElt++) {
    m_FittedAmpTplratio[idx][iElt] = m_FittedAmpTplratio[idx - 1][iElt];
    m_FittedErrorTplratio[idx][iElt] = m_FittedErrorTplratio[idx - 1][iElt];
    m_DtmTplratio[idx][iElt] = m_DtmTplratio[idx - 1][iElt];
    m_MtmTplratio[idx][iElt] = m_MtmTplratio[idx - 1][iElt];
    m_LyaAsymCoeffTplratio[idx][iElt] = m_LyaAsymCoeffTplratio[idx - 1][iElt];
    m_LyaWidthCoeffTplratio[idx][iElt] = m_LyaWidthCoeffTplratio[idx - 1][iElt];
    m_LyaDeltaCoeffTplratio[idx][iElt] = m_LyaDeltaCoeffTplratio[idx - 1][iElt];
    m_LyaIgmIdxTplratio[idx][iElt] = m_LyaIgmIdxTplratio[idx - 1][iElt];
    m_LinesLogPriorTplratio[idx][iElt] = m_LinesLogPriorTplratio[idx - 1][iElt];
  }
  return;
}

void CTplratioManager::initTplratioCatalogs(Int32 opt_tplratio_ismFit) {
  // TODO: use the passed tplRatioCatalog
  // TODO: check if m_CatalogTplRatio changes between iterations

  m_LineCatalogLinesCorrespondingNominalAmp =
      m_CatalogTplRatio->InitLineCorrespondingAmplitudes(
          m_Elements, opt_tplratio_ismFit,
          m_continuumManager->getIsmCorrectionFromTpl(), m_NSigmaSupport);
  m_opt_dust_calzetti = opt_tplratio_ismFit;
  SetMultilineNominalAmplitudesFast(0);
  // m_RestLineList = m_CatalogTplRatio->GetRestLinesList(0);
  // LoadCatalog(m_RestLineList);
  // LogCatalogInfos();

  // Resize tplratio buffers
  m_ChisquareTplratio.resize(m_CatalogTplRatio->GetCatalogsCount());
  m_ScaleMargCorrTplratio.resize(m_CatalogTplRatio->GetCatalogsCount());
  m_StrongELPresentTplratio.resize(m_CatalogTplRatio->GetCatalogsCount());
  m_StrongHalphaELPresentTplratio.resize(m_CatalogTplRatio->GetCatalogsCount());
  m_NLinesAboveSNRTplratio.resize(m_CatalogTplRatio->GetCatalogsCount());
  m_FittedAmpTplratio.resize(m_CatalogTplRatio->GetCatalogsCount());
  m_LyaAsymCoeffTplratio.resize(m_CatalogTplRatio->GetCatalogsCount());
  m_LyaWidthCoeffTplratio.resize(m_CatalogTplRatio->GetCatalogsCount());
  m_LyaDeltaCoeffTplratio.resize(m_CatalogTplRatio->GetCatalogsCount());
  m_LyaIgmIdxTplratio.resize(m_CatalogTplRatio->GetCatalogsCount());
  m_FittedErrorTplratio.resize(m_CatalogTplRatio->GetCatalogsCount());
  m_MtmTplratio.resize(m_CatalogTplRatio->GetCatalogsCount());
  m_DtmTplratio.resize(m_CatalogTplRatio->GetCatalogsCount());
  m_LinesLogPriorTplratio.resize(m_CatalogTplRatio->GetCatalogsCount());
  for (Int32 ktplratio = 0; ktplratio < m_CatalogTplRatio->GetCatalogsCount();
       ktplratio++) {
    m_FittedAmpTplratio[ktplratio].resize(m_Elements.size());
    m_FittedErrorTplratio[ktplratio].resize(m_Elements.size());
    m_MtmTplratio[ktplratio].resize(m_Elements.size());
    m_DtmTplratio[ktplratio].resize(m_Elements.size());
    m_LyaAsymCoeffTplratio[ktplratio].resize(m_Elements.size());
    m_LyaWidthCoeffTplratio[ktplratio].resize(m_Elements.size());
    m_LyaDeltaCoeffTplratio[ktplratio].resize(m_Elements.size());
    m_LyaIgmIdxTplratio[ktplratio].resize(m_Elements.size());
    m_LinesLogPriorTplratio[ktplratio].resize(m_Elements.size());
  }

  m_tplratioLeastSquareFast = false;
}

/**
 * @brief CLineModelFitting::SetMultilineNominalAmplitudesFast
 * This method sets the linemodel unique elt nominal amplitudes to the
 * corresponding value of the iCatalog st catalog. INFO: fast method,
 * InitLineCorrespondence() should have been called previously with the same
 * LineModelElementList arg.
 * @param iCatalog
 * @return
 */
bool CTplratioManager::SetMultilineNominalAmplitudesFast(Int32 iCatalog) {
  if (iCatalog < 0) {
    return false;
  }
  Float64 nominalAmp = 0.0;
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    Int32 nLines = m_Elements[iElts]->GetSize();
    for (Int32 j = 0; j < nLines; j++) {
      nominalAmp =
          m_LineCatalogLinesCorrespondingNominalAmp[iElts][iCatalog][j];
      m_Elements[iElts]->SetNominalAmplitude(j, nominalAmp);
    }
  }
  return true;
}

const std::string &CTplratioManager::getTplratio_bestTplName() const {
  return m_tplratioBestTplName;
}
Float64 CTplratioManager::getTplratio_bestTplIsmCoeff() const {
  return m_tplratioBestTplIsmCoeff;
}

Float64 CTplratioManager::getTplratio_bestAmplitude() const {
  return m_tplratioBestTplAmplitude;
}

Float64 CTplratioManager::getTplratio_bestDtm() const {
  return m_tplratioBestTplDtm;
}

Float64 CTplratioManager::getTplratio_bestMtm() const {
  return m_tplratioBestTplMtm;
}

void CTplratioManager::logParameters() {
  CLineRatioManager::logParameters();
  Log.LogDetail(Formatter() << " m_opt_haprior" << m_opt_haprior);
  Log.LogDetail(Formatter() << "NSigmaSupport=" << m_NSigmaSupport);
  Log.LogDetail(Formatter() << " m_opt_firstpass_forcedisableTplratioISMfit "
                            << m_opt_firstpass_forcedisableTplratioISMfit);
  Log.LogDetail(Formatter() << "forcedisableTplratioISMfit="
                            << m_forcedisableTplratioISMfit);
  Log.LogDetail(Formatter() << "tplratioBestTplName=" << m_tplratioBestTplName);
  Log.LogDetail(Formatter()
                << "tplratioBestTplIsmCoeff=" << m_tplratioBestTplIsmCoeff);
  Log.LogDetail(Formatter()
                << "tplratioBestTplAmplitude=" << m_tplratioBestTplAmplitude);
  Log.LogDetail(Formatter() << "tplratioBestTplDtm=" << m_tplratioBestTplDtm);
  Log.LogDetail(Formatter() << "tplratioBestTplMtm=" << m_tplratioBestTplMtm);
  Log.LogDetail(
      Formatter()
      << "tplratioLeastSquareFast="
      << m_tplratioLeastSquareFast); // for rigidity=tplratio: switch to
                                     // use fast least square estimation
}

void CTplratioManager::setPassMode(Int32 iPass) {
  CLineRatioManager::setPassMode(iPass);
  if (iPass == 1) {
    m_forcedisableTplratioISMfit = m_opt_firstpass_forcedisableTplratioISMfit;
  }
  if (iPass == 2) {

    m_forcedisableTplratioISMfit = false;
  }
}

void CTplratioManager::SetForcedisableTplratioISMfit(bool opt) {
  m_forcedisableTplratioISMfit = opt;
}
Int32 CTplratioManager::getTplratio_count() const {
  return m_CatalogTplRatio->GetCatalogsCount();
}

TFloat64List CTplratioManager::getTplratio_priors() {
  return m_CatalogTplRatio->getCatalogsPriors();
}

const TFloat64List &CTplratioManager::GetChisquareTplratio() const {
  return m_ChisquareTplratio;
}

/**
 * @brief CTplratioManager::GetPriorLinesTplratio
 * WARNING: as stated in fit(), the prior is valid in this code structure only
 * for tplratio with only 1 element containing the EL component, hence using
 * idx=0 here
 * @return
 */
TFloat64List CTplratioManager::GetPriorLinesTplratio() const {
  TFloat64List plinestplratio;
  Int32 eltIdx = 0;
  for (Int32 ktpl = 0; ktpl < m_LinesLogPriorTplratio.size(); ktpl++) {
    plinestplratio.push_back(m_LinesLogPriorTplratio[ktpl][eltIdx]);
  }
  return plinestplratio;
}

const TFloat64List &CTplratioManager::GetScaleMargTplratio() const {
  return m_ScaleMargCorrTplratio;
}

const TBoolList &CTplratioManager::GetStrongELPresentTplratio() const {
  return m_StrongELPresentTplratio;
}

const TBoolList &CTplratioManager::getHaELPresentTplratio() const {
  return m_StrongHalphaELPresentTplratio;
}

const TInt32List &CTplratioManager::GetNLinesAboveSNRTplratio() const {
  return m_NLinesAboveSNRTplratio;
}

bool CTplratioManager::setTplratioAmplitude(const TFloat64List &ampsElts,
                                            const TFloat64List &errorsElts) {
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    m_Elements[iElts]->SetFittedAmplitude(ampsElts[iElts], errorsElts[iElts]);
  }
  return true;
}

void CTplratioManager::SetLeastSquareFastEstimationEnabled(Int32 enabled) {
  m_tplratioLeastSquareFast = enabled;
}

/**
 * @brief
 *
 * @param idx
 * @param _merit
 * @param _meritprior
 */
void CTplratioManager::updateTplratioResults(
    Int32 idx, Float64 _merit, Float64 _meritprior,
    TFloat64List &bestTplratioMerit, TFloat64List &bestTplratioMeritPrior) {
  bestTplratioMerit[idx] = _merit;
  bestTplratioMeritPrior[idx] = _meritprior;
  m_ChisquareTplratio[idx] = _merit;
  m_ScaleMargCorrTplratio[idx] = m_Elements.getScaleMargCorrection(idx);
  m_StrongELPresentTplratio[idx] =
      m_Elements.GetModelStrongEmissionLinePresent();
  // given that Ha is a strong emission line,
  if (m_opt_haprior > 0.) // check first that haprior is activated
    m_StrongHalphaELPresentTplratio[idx] =
        m_StrongELPresentTplratio[idx] ? m_Elements.GetModelHaStrongest()
                                       : false; // result per tplratio

  TStringList strongELSNRAboveCut; // = getLinesAboveSNR(3.5); //this
                                   // is costing a lot of processing
                                   // time, so deactivated for now.
  m_NLinesAboveSNRTplratio[idx] = strongELSNRAboveCut.size();

  // Saving the model A, errorA, and dtm, mtm, ... (for all tplratios,
  // needed ?) NB: this is only needed for the index=savedIdxFitted
  // ultimately
  for (Int32 iElt = 0; iElt < m_Elements.size(); iElt++) {
    bool savedAmp = false;
    bool allampzero = true;
    m_FittedAmpTplratio[idx][iElt] = NAN;
    m_FittedErrorTplratio[idx][iElt] = NAN;
    m_DtmTplratio[idx][iElt] = NAN;
    m_MtmTplratio[idx][iElt] = NAN;
    m_LyaAsymCoeffTplratio[idx][iElt] = NAN;
    m_LyaWidthCoeffTplratio[idx][iElt] = NAN;
    m_LyaDeltaCoeffTplratio[idx][iElt] = NAN;
    m_LyaIgmIdxTplratio[idx][iElt] = undefIdx;
    m_LinesLogPriorTplratio[idx][iElt] = _meritprior;
    Int32 nLines = m_Elements[iElt]->GetSize();

    for (Int32 j = 0; j < nLines; j++) {
      if (savedAmp) {
        break;
      }
      Float64 amp = m_Elements[iElt]->GetFittedAmplitude(j);
      if (amp > 0) {
        allampzero = false;
      }
      if (amp > 0 && !m_Elements[iElt]->IsOutsideLambdaRange(j)) {

        Float64 amp_error = m_Elements[iElt]->GetFittedAmplitudeErrorSigma(j);
        Float64 nominal_amp = m_Elements[iElt]->GetNominalAmplitude(j);
        m_FittedAmpTplratio[idx][iElt] = amp / nominal_amp;
        Log.LogDebug("    model : fit tplratio mode, tplratio_fittedamp: %e",
                     m_FittedAmpTplratio[idx][iElt]);

        m_FittedErrorTplratio[idx][iElt] = amp_error / nominal_amp;
        m_DtmTplratio[idx][iElt] = m_Elements[iElt]->GetSumCross();
        m_MtmTplratio[idx][iElt] = m_Elements[iElt]->GetSumGauss();

        TAsymParams params = m_Elements[iElt]->GetAsymfitParams(0);
        m_LyaAsymCoeffTplratio[idx][iElt] = params.alpha;
        m_LyaWidthCoeffTplratio[idx][iElt] = params.sigma;
        m_LyaDeltaCoeffTplratio[idx][iElt] = params.delta;

        TSymIgmParams params_igm = m_Elements[iElt]->GetSymIgmParams(0);
        m_LyaIgmIdxTplratio[idx][iElt] = params_igm.m_igmidx;

        savedAmp = true;
        break;
      }
    }
    // TODO: this case should be treated more
    // carefully, save dtm, mtm, and more...
    if (allampzero && !savedAmp)
      m_FittedAmpTplratio[idx][iElt] = 0.0;
  }
  return;
}

Float64 CTplratioManager::computelogLinePriorMerit(
    Int32 itratio,
    const std::vector<CPriorHelper::SPriorTZE> &logPriorDataTplRatio) {

  if (!logPriorDataTplRatio.size())
    return 0.;

  // lines prior
  Float64 _meritprior = -2. * logPriorDataTplRatio[itratio].betaTE *
                        logPriorDataTplRatio[itratio].logprior_precompTE;
  _meritprior += -2. * logPriorDataTplRatio[itratio].betaA *
                 logPriorDataTplRatio[itratio].logprior_precompA;
  _meritprior += -2. * logPriorDataTplRatio[itratio].betaZ *
                 logPriorDataTplRatio[itratio].logprior_precompZ;
  if (logPriorDataTplRatio[itratio].A_sigma <= 0.0)
    return _meritprior;

  Float64 ampl = 0.0;
  for (Int32 iElt = 0; iElt < m_Elements.size(); iElt++) {
    bool foundAmp = false;
    Int32 nLines = m_Elements[iElt]->GetSize();
    for (Int32 j = 0; j < nLines; j++) {
      Float64 amp = m_Elements[iElt]->GetFittedAmplitude(j);
      if (amp > 0 && !m_Elements[iElt]->IsOutsideLambdaRange(j)) {
        Float64 nominal_amp = m_Elements[iElt]->GetNominalAmplitude(j);
        ampl = amp / nominal_amp;
        foundAmp = true;
        break;
      }
    }
    /*if (foundAmp)
      break;
      //Didier: probably this is missing here??? I suppose we are
      // looking for the first non-null and valid amplitude?
     */
  }
  _meritprior += logPriorDataTplRatio[itratio].betaA *
                 (ampl - logPriorDataTplRatio[itratio].A_mean) *
                 (ampl - logPriorDataTplRatio[itratio].A_mean) /
                 (logPriorDataTplRatio[itratio].A_sigma *
                  logPriorDataTplRatio[itratio].A_sigma);

  return _meritprior;
}

Float64 CTplratioManager::computeMerit(Int32 itratio) {
  // TODO enableLogging should be renamed and be a CRigidityManager member,
  // initialized with parameter store
  bool enableLogging = true;
  Float64 _merit;
  if (!enableLogging && m_tplratioLeastSquareFast)
    _merit = getLeastSquareMeritFast();
  else {
    m_model->refreshModel();
    _merit = getLeastSquareMerit();
  }

  Float64 _meritprior =
      computelogLinePriorMerit(itratio, m_logPriorDataTplRatio);

  if (_merit + _meritprior <
      m_bestTplratioMerit[itratio] + m_bestTplratioMeritPrior[itratio]) {
    // update result variables
    updateTplratioResults(itratio, _merit, _meritprior, m_bestTplratioMerit,
                          m_bestTplratioMeritPrior);
  }
  return _merit;
}

void CTplratioManager::finish(Float64 redshift) {
  bool retSetMultiAmplFast =
      SetMultilineNominalAmplitudesFast(m_savedIdxFitted);
  if (!retSetMultiAmplFast) {
    Log.LogError("Linemodel: tplratio, Unable to set Multiline "
                 "NominalAmplitudes from Tplratio !");
  }

  // Set the velocities from templates: todo auto switch when velfit is ON
  // m_CatalogTplRatio->GetCatalogVelocities(savedIdxFitted,
  // m_velocityEmission, m_velocityAbsorption);
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    Log.LogDetail("    model - Linemodel: tplratio = %d (%s, with "
                  "ebmv=%.3f), and A=%e",
                  m_savedIdxFitted, m_tplratioBestTplName.c_str(),
                  m_tplratioBestTplIsmCoeff,
                  m_FittedAmpTplratio[m_savedIdxFitted][iElts]);
    m_Elements[iElts]->SetFittedAmplitude(
        m_FittedAmpTplratio[m_savedIdxFitted][iElts],
        m_FittedErrorTplratio[m_savedIdxFitted][iElts]);
    m_Elements[iElts]->SetSumCross(m_DtmTplratio[m_savedIdxFitted][iElts]);
    m_Elements[iElts]->SetSumGauss(m_MtmTplratio[m_savedIdxFitted][iElts]);
  }

  // Lya
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++)
    m_Elements[iElts]->SetAsymfitParams(
        {m_LyaWidthCoeffTplratio[m_savedIdxFitted][iElts],
         m_LyaAsymCoeffTplratio[m_savedIdxFitted][iElts],
         m_LyaDeltaCoeffTplratio[m_savedIdxFitted][iElts]});

  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++)
    m_Elements[iElts]->SetSymIgmParams(
        TSymIgmParams(m_LyaIgmIdxTplratio[m_savedIdxFitted][iElts], redshift));

  m_model->refreshModel();
}

void CTplratioManager::saveResults(Int32 itratio) {
  m_savedIdxFitted = itratio;
  m_tplratioBestTplName = m_CatalogTplRatio->GetCatalogName(m_savedIdxFitted);
  m_tplratioBestTplIsmCoeff = GetIsmCoeff(m_savedIdxFitted);
  m_tplratioBestTplAmplitude =
      m_FittedAmpTplratio[m_savedIdxFitted][0]; // Should be only 1 elt
  // in tpl ratio mode...
  m_tplratioBestTplDtm =
      m_DtmTplratio[m_savedIdxFitted]
                   [0]; // Should be only 1 elt in tpl ratio mode...
  m_tplratioBestTplMtm =
      m_MtmTplratio[m_savedIdxFitted]
                   [0]; // Should be only 1 elt in tpl ratio mode...
}

Float64 CTplratioManager::GetIsmCoeff(Int32 idx) const {
  if (m_continuumManager->getIsmCorrectionFromTpl() == nullptr &&
      m_opt_dust_calzetti)
    THROWG(
        INTERNAL_ERROR,
        "ismCorrectionCalzetti is not loaded while tplratio_ism is activated");
  if (!m_opt_dust_calzetti)
    return NAN;
  return m_continuumManager->getIsmCorrectionFromTpl()->GetEbmvValue(
      m_CatalogTplRatio->GetIsmIndex(idx));
}

bool CTplratioManager::setTplratioModel(Int32 itplratio,
                                        bool enableSetVelocity) {
  SetMultilineNominalAmplitudesFast(itplratio);

  /* TODO reactivate this if once called with enableSetVelocity=true . ->
  velocities must be imported from linemodelfitting if (enableSetVelocity) {
    // Set the velocities from templates: todo auto switch when velfit is ON
    m_CatalogTplRatio.GetCatalogVelocities(itplratio, m_velocityEmission,
                                           m_velocityAbsorption);
  }
  */

  Log.LogDebug("    model : setTplratioModel, loaded: %d = %s", itplratio,
               m_CatalogTplRatio->GetCatalogName(itplratio).c_str());
  return true;
}
