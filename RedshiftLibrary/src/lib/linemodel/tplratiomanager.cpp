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
    const CLMEltListVectorPtr &elementsVector, const CSpcModelVectorPtr &models,
    const CCSpectrumVectorPtr &inputSpcs,
    const CTLambdaRangePtrVector &lambdaRanges,
    std::shared_ptr<CContinuumManager> continuumManager,
    const CLineMap &restLineList, const std::shared_ptr<Int32> &curObs)
    : CLineRatioManager(elementsVector, models, inputSpcs, lambdaRanges,
                        continuumManager, restLineList, curObs) {
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
  m_MeritTplratio.assign(ntplratio, INFINITY);
  m_PriorMeritTplratio.assign(ntplratio, 0.0);
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

Int32 CTplratioManager::prepareFit(Float64 redshift) {
  m_logPriorDataTplRatio.clear();
  m_savedIdxFitted = undefIdx;
  Int32 ntplratio = m_CatalogTplRatio->GetCatalogsCount();
  initMerit(ntplratio);
  if (!m_tplratio_priorhelper->mInitFailed) {
    // prior initilization for tplratio EL only
    if (getElementList().size() > 1)
      THROWG(INTERNAL_ERROR, "model: Unable to use tplratio line priors "
                             "with nElts>1 for now");
    // NB: this could be done if the EL element idx in searched (see later
    // in the itratio loop, UV Abs lines would be not affected by priors
    // then)

    for (Int32 itratio = 0; itratio < ntplratio; itratio++) {
      // prepare the lines prior data
      Int32 ebvfilter = m_CatalogTplRatio->GetIsmIndex(itratio);

      std::string tplrationame = m_CatalogTplRatio->GetCatalogName(itratio);
      m_logPriorDataTplRatio.push_back(m_tplratio_priorhelper->GetTZEPriorData(
          tplrationame, ebvfilter, redshift));
    }
  }
  return ntplratio;
}

bool CTplratioManager::init(Float64 redshift, Int32 itratio) {
  if (m_forcedisableTplratioISMfit && itratio > 0 &&
      m_CatalogTplRatio->GetIsmIndex(itratio) > 0) {
    duplicateTplratioResult(itratio);
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
void CTplratioManager::duplicateTplratioResult(Int32 idx) {
  m_PriorMeritTplratio[idx] = m_PriorMeritTplratio[idx - 1];
  m_MeritTplratio[idx] = m_MeritTplratio[idx - 1];
  m_ScaleMargCorrTplratio[idx] = m_ScaleMargCorrTplratio[idx - 1];
  m_StrongELPresentTplratio[idx] = m_StrongELPresentTplratio[idx - 1];
  m_StrongHalphaELPresentTplratio[idx] =
      m_StrongHalphaELPresentTplratio[idx - 1];
  m_NLinesAboveSNRTplratio[idx] = m_NLinesAboveSNRTplratio[idx - 1];

  for (Int32 iElt = 0; iElt < getElementList().size(); iElt++) {
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

  m_LineCatalogCorrespondingNominalAmp =
      m_CatalogTplRatio->InitLineCorrespondingAmplitudes(
          getElementList(), opt_tplratio_ismFit,
          m_continuumManager->getIsmCorrectionFromTpl());
  m_opt_dust_calzetti = opt_tplratio_ismFit;
  m_tplratioBestTplName = undefStr;
  Int32 s = m_CatalogTplRatio->GetCatalogsCount();
  Int32 elCount = getElementList().size();
  // Resize tplratio buffers
  m_MeritTplratio.assign(s, NAN);
  m_ScaleMargCorrTplratio.assign(s, NAN);
  m_StrongELPresentTplratio.assign(s, false);
  m_StrongHalphaELPresentTplratio.assign(s, false);
  m_NLinesAboveSNRTplratio.assign(s, -1);
  m_FittedAmpTplratio.assign(s, TFloat64List(elCount, NAN));
  m_LyaAsymCoeffTplratio.assign(s, TFloat64List(elCount, NAN));
  m_LyaWidthCoeffTplratio.assign(s, TFloat64List(elCount, NAN));
  m_LyaDeltaCoeffTplratio.assign(s, TFloat64List(elCount, NAN));
  m_LyaIgmIdxTplratio.assign(s, TInt32List(elCount, -1));
  m_FittedErrorTplratio.assign(s, TFloat64List(elCount, NAN));
  m_MtmTplratio.assign(s, TFloat64List(elCount, NAN));
  m_DtmTplratio.assign(s, TFloat64List(elCount, NAN));
  m_LinesLogPriorTplratio.assign(s, TFloat64List(elCount, 0.));

  m_tplratioLeastSquareFast = false;
}

/**
 * @brief CTplratioManager::SetNominalAmplitudes
 * This method sets the linemodel unique elt nominal amplitudes to the
 * corresponding value of the iCatalog st catalog. INFO: fast method,
 * InitLineCorrespondence() should have been called previously with the same
 * LineModelElementList arg.
 * @param iCatalog
 * @return
 */
void CTplratioManager::SetNominalAmplitudes(Int32 iCatalog) {
  if (iCatalog < 0 || iCatalog >= m_CatalogTplRatio->GetCatalogsCount())
    THROWG(INTERNAL_ERROR, Formatter()
                               << "wrong line catalog index: " << iCatalog);
  for (Int32 elt_index = 0; elt_index != getElementList().size(); ++elt_index) {
    auto &elt_ptr = getElementList()[elt_index];
    for (Int32 line_index = 0; line_index != elt_ptr->GetSize(); ++line_index) {
      Float64 const nominalAmp =
          m_LineCatalogCorrespondingNominalAmp[iCatalog][elt_index][line_index];
      elt_ptr->SetNominalAmplitude(line_index, nominalAmp);
    }
  }
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
  Log.LogDetail(Formatter() << "tplratioBestTplAmplitudeEm="
                            << m_tplratioBestTplAmplitudeEm);
  Log.LogDetail(Formatter() << "tplratioBestTplAmplitudeAbs="
                            << m_tplratioBestTplAmplitudeAbs);
  Log.LogDetail(Formatter()
                << "tplratioBestTplDtmEm=" << m_tplratioBestTplDtmEm);
  Log.LogDetail(Formatter()
                << "tplratioBestTplDtmAbs=" << m_tplratioBestTplDtmAbs);

  Log.LogDetail(Formatter()
                << "tplratioBestTplMtmEm=" << m_tplratioBestTplMtmEm);
  Log.LogDetail(Formatter()
                << "tplratioBestTplMtmAbs=" << m_tplratioBestTplMtmAbs);
  Log.LogDetail(
      Formatter()
      << "tplratioLeastSquareFast="
      << m_tplratioLeastSquareFast); // for rigidity=tplratio: switch to
                                     // use fast least square estimation
}

void CTplratioManager::setPassMode(Int32 iPass) {
  CLineRatioManager::setPassMode(iPass);
  if (iPass == 1)
    m_forcedisableTplratioISMfit = m_opt_firstpass_forcedisableTplratioISMfit;

  if (iPass == 2)
    m_forcedisableTplratioISMfit = false;
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
  return m_MeritTplratio;
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

void CTplratioManager::SetLeastSquareFastEstimationEnabled(Int32 enabled) {
  m_tplratioLeastSquareFast = enabled;
}
// moved into a fct to reduce cognitive complexity of updateTplratioResults
void CTplratioManager::fillHalphaArray(Int32 idx) {
  // check first that haprior is activated
  if (m_opt_haprior <= 0.)
    return;
  m_StrongHalphaELPresentTplratio[idx] =
      m_StrongELPresentTplratio[idx] ? getElementList().GetModelHaStrongest()
                                     : false; // result per tplratio
}
/**
 * @brief
 *
 * @param idx
 * @param _merit
 * @param _meritprior
 */
void CTplratioManager::updateTplratioResults(Int32 idx, Float64 _merit,
                                             Float64 _meritprior) {
  m_PriorMeritTplratio[idx] = _meritprior;
  m_MeritTplratio[idx] = _merit;
  m_ScaleMargCorrTplratio[idx] = getElementList().getScaleMargCorrection();
  m_StrongELPresentTplratio[idx] =
      getElementList().GetModelStrongEmissionLinePresent();
  // given that Ha is a strong emission line,
  fillHalphaArray(idx);

  TStringList strongELSNRAboveCut; // = getLinesAboveSNR(3.5); //this
                                   // is costing a lot of processing
                                   // time, so deactivated for now.
  m_NLinesAboveSNRTplratio[idx] = strongELSNRAboveCut.size();

  Int32 s = getElementList().size();
  // reinit
  m_FittedAmpTplratio[idx].assign(s, NAN);
  m_FittedErrorTplratio[idx].assign(s, NAN);
  m_DtmTplratio[idx].assign(s, NAN);
  m_MtmTplratio[idx].assign(s, NAN);
  m_LyaAsymCoeffTplratio[idx].assign(s, NAN);
  m_LyaWidthCoeffTplratio[idx].assign(s, NAN);
  m_LyaDeltaCoeffTplratio[idx].assign(s, NAN);
  m_LyaIgmIdxTplratio[idx].assign(s, undefIdx);
  m_LinesLogPriorTplratio[idx].assign(s, _meritprior);
  // Saving the model A, errorA, and dtm, mtm, ... (for all tplratios,
  // needed ?) NB: this is only needed for the index=savedIdxFitted
  // ultimately
  for (Int32 iElt = 0; iElt < s; iElt++) {
    bool savedAmp = false;
    bool allampzero = true;
    auto const &elt = getElementList()[iElt];
    for (Int32 line_idx = 0; line_idx != elt->GetSize(); ++line_idx) {
      Float64 amp = elt->GetFittedAmplitude(line_idx);
      // we need to take into case where lineratioAmplitude is NAN
      if (isnan(amp) || amp <= 0. || elt->IsOutsideLambdaRange(line_idx))
        continue;
      allampzero = false;

      Float64 amp_error = elt->GetFittedAmplitudeErrorSigma(line_idx);
      Float64 nominal_amp = elt->GetNominalAmplitude(line_idx);
      m_FittedAmpTplratio[idx][iElt] = amp / nominal_amp;
      Log.LogDebug("    model : fit tplratio mode, tplratio_fittedamp: %e",
                   m_FittedAmpTplratio[idx][iElt]);

      m_FittedErrorTplratio[idx][iElt] = amp_error / nominal_amp;
      m_DtmTplratio[idx][iElt] = elt->GetSumCross();
      m_MtmTplratio[idx][iElt] = elt->GetSumGauss();

      TAsymParams params = elt->GetAsymfitParams(0);
      m_LyaAsymCoeffTplratio[idx][iElt] = params.alpha;
      m_LyaWidthCoeffTplratio[idx][iElt] = params.sigma;
      m_LyaDeltaCoeffTplratio[idx][iElt] = params.delta;

      TSymIgmParams params_igm = elt->GetSymIgmParams(0);
      m_LyaIgmIdxTplratio[idx][iElt] = params_igm.m_igmidx;

      savedAmp = true;
      break;
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
  for (const auto &elt : getElementList()) {
    bool foundAmp = false;
    for (Int32 line_idx = 0; line_idx != elt->GetSize(); ++line_idx) {
      Float64 amp = elt->GetFittedAmplitude(line_idx);
      if (amp <= 0. || elt->IsOutsideLambdaRange(line_idx))
        continue;
      Float64 nominal_amp = elt->GetNominalAmplitude(line_idx);
      ampl = amp / nominal_amp;
      foundAmp = true;
      break;
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
  Float64 _merit = 0;
  Float64 _meritprior = 0;
  /*if (!enableLogging && m_tplratioLeastSquareFast)
    _merit = getLeastSquareMeritFast();
  else */

  for (*m_curObs = 0; *m_curObs < m_models->size(); (*m_curObs)++) {
    getModel().refreshModel();
    _merit += getLeastSquareMerit();

    _meritprior += computelogLinePriorMerit(itratio, m_logPriorDataTplRatio);
  }
  //  *m_curObs = 0;
  if (_merit + _meritprior <
      m_MeritTplratio[itratio] + m_PriorMeritTplratio[itratio]) {
    // update result variables
    updateTplratioResults(itratio, _merit, _meritprior);
  }
  return _merit;
}

void CTplratioManager::resetToBestRatio(Float64 redshift) {

  // first reinit all the elements:
  init(redshift, m_savedIdxFitted);

  for (Int32 iElts = 0; iElts < getElementList().size(); iElts++) {
    Log.LogDetail("    model - Linemodel: tplratio = %d (%s, with "
                  "ebmv=%.3f), and A=%e",
                  m_savedIdxFitted, m_tplratioBestTplName.c_str(),
                  m_tplratioBestTplIsmCoeff,
                  m_FittedAmpTplratio[m_savedIdxFitted][iElts]);
    getElementList()[iElts]->SetElementAmplitude(
        m_FittedAmpTplratio[m_savedIdxFitted][iElts],
        m_FittedErrorTplratio[m_savedIdxFitted][iElts]);
    getElementList()[iElts]->SetSumCross(
        m_DtmTplratio[m_savedIdxFitted][iElts]);
    getElementList()[iElts]->SetSumGauss(
        m_MtmTplratio[m_savedIdxFitted][iElts]);
  }

  // Lya
  for (Int32 iElts = 0; iElts < getElementList().size(); iElts++)
    getElementList()[iElts]->SetAsymfitParams(
        {m_LyaWidthCoeffTplratio[m_savedIdxFitted][iElts],
         m_LyaAsymCoeffTplratio[m_savedIdxFitted][iElts],
         m_LyaDeltaCoeffTplratio[m_savedIdxFitted][iElts]});

  for (Int32 iElts = 0; iElts < getElementList().size(); iElts++)
    getElementList()[iElts]->SetSymIgmParams(
        TSymIgmParams(m_LyaIgmIdxTplratio[m_savedIdxFitted][iElts], redshift));

  getModel().refreshModel();
}

void CTplratioManager::saveResults(Int32 itratio) {
  m_savedIdxFitted = itratio;
  m_tplratioBestTplName = m_CatalogTplRatio->GetCatalogName(m_savedIdxFitted);
  m_tplratioBestTplIsmCoeff = GetIsmCoeff(m_savedIdxFitted);
  TInt32List idx_em =
      getElementList().findElementTypeIndices(CLine::EType::nType_Emission);
  if (!idx_em.empty()) {
    m_tplratioBestTplAmplitudeEm =
        m_FittedAmpTplratio[m_savedIdxFitted]
                           [idx_em.front()]; // Should be only 1 elt in tpl
                                             // ratio mode...
    m_tplratioBestTplDtmEm =
        m_DtmTplratio[m_savedIdxFitted][idx_em.front()]; // Should be only 1 elt
                                                         // in tpl ratio mode...
    m_tplratioBestTplMtmEm =
        m_MtmTplratio[m_savedIdxFitted][idx_em.front()]; // Should be only 1 elt
                                                         // in tpl ratio mode...
  }
  TInt32List idx_abs =
      getElementList().findElementTypeIndices(CLine::EType::nType_Absorption);
  if (!idx_abs.empty()) {
    m_tplratioBestTplAmplitudeAbs =
        m_FittedAmpTplratio[m_savedIdxFitted]
                           [idx_abs.front()]; // Should be only 1 elt in tpl
                                              // ratio mode...
    m_tplratioBestTplDtmAbs =
        m_DtmTplratio[m_savedIdxFitted]
                     [idx_abs.front()]; // Should be only 1 elt
                                        // in tpl ratio mode...
    m_tplratioBestTplMtmAbs =
        m_MtmTplratio[m_savedIdxFitted]
                     [idx_abs.front()]; // Should be only 1 elt
                                        // in tpl ratio mode...
  }
}

Float64 CTplratioManager::GetIsmCoeff(Int32 idx) const {
  if (m_continuumManager->getIsmCorrectionFromTpl() == nullptr &&
      m_opt_dust_calzetti)
    THROWG(INTERNAL_ERROR, "ismCorrectionCalzetti is not loaded while "
                           "tplratio_ism is activated");
  if (!m_opt_dust_calzetti)
    return NAN;
  return m_continuumManager->getIsmCorrectionFromTpl()->GetEbmvValue(
      m_CatalogTplRatio->GetIsmIndex(idx));
}

bool CTplratioManager::setTplratioModel(Int32 itplratio,
                                        bool enableSetVelocity) {
  SetNominalAmplitudes(itplratio);

  /* TODO reactivate this if once called with enableSetVelocity=true . ->
  velocities must be imported from linemodelfitting if (enableSetVelocity) {
    // Set the velocities from templates: todo auto switch when velfit is ON
    m_CatalogTplRatio.GetCatalogVelocities(itplratio, m_velocityEmission,
                                           m_velocityAbsorption);
  }


  Log.LogDebug("    model : setTplratioModel, loaded: %d = %s", itplratio,
               m_CatalogTplRatio->GetCatalogName(itplratio).c_str());*/
  return true;
}
