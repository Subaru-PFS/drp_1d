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

#include "RedshiftLibrary/linemodel/lineratiomanager.h"
#include "RedshiftLibrary/linemodel/abstractfitter.h"
#include "RedshiftLibrary/linemodel/continuummanager.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/linemodel/rulesmanager.h"
#include "RedshiftLibrary/linemodel/spectrummodel.h"
#include "RedshiftLibrary/linemodel/tplcorrmanager.h"
#include "RedshiftLibrary/linemodel/tplratiomanager.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

using namespace NSEpic;

CLineRatioManager::CLineRatioManager(
    const std::shared_ptr<CLMEltListVector> &elementsVector,
    const CSpcModelVectorPtr &models, const CCSpectrumVectorPtr &inputSpcs,
    const CTLambdaRangePtrVector &lambdaRanges,
    std::shared_ptr<CContinuumManager> continuumManager,
    const CLineMap &restLineList, const CSpectraGlobalIndex &spcIndex)
    : m_elementsVector(elementsVector), m_inputSpcs(inputSpcs),
      m_lambdaRanges(lambdaRanges), m_models(models), m_spectraIndex(spcIndex),
      m_continuumManager(continuumManager), m_RestLineList(restLineList) {

  CAutoScope autoscope(Context.m_ScopeStack, "lineModel");
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();

  bool useAsymProfile = ps->GetScoped<std::string>("lya.profile") == "asym";
  if (Context.GetCurrentMethod() == "lineModelSolve" && useAsymProfile) {
    m_opt_lya_forcefit =
        ps->GetScoped<bool>("lya.asymProfile.switchFixedToFit");
    m_opt_lya_forcedisablefit =
        ps->GetScoped<bool>("lya.asymProfile.switchFitToFixed");
  }
}

/**
 * @brief CLineModelFitting::setLyaProfile
 * If a Lya line is present with SYMIGM profile, fit igmIdx
 * If a Lya line is present with ASYMFIT profile, fit the width and asymmetry
 * parameters If a Lya line is present with ASYMFIXED profile, set the width
 * and asymmetry parameters according to profile parameters given in the
 * string
 * @param redshift, catalog (could be tplratio or linecatalog)
 * @return 1 if successfully fitted, 0 if error, 2 if Lya not present, 3 if Lya
 * not configured to be fitted in the catalog
 * note: this applies to all asym or symigm profiles (for symigm: all lines
 * below Lya)
 */

void CLineRatioManager::setLyaProfile(Float64 redshift,
                                      const CLineMap &catalog) {
  auto const indices_Igm = m_elementsVector->getIgmLinesIndices();
  if (indices_Igm.empty())
    return;

  // ASymProfile
  {
    // assuming only one asymfit/fixed profile
    auto const &[elt_idx_LyaE, line_indices_LyaE] = indices_Igm.front();
    Int32 line_idx_LyaE = line_indices_LyaE.front();

    auto const &param_LyaE = m_elementsVector->getElementParam()[elt_idx_LyaE];
    auto const &profile = param_LyaE->getLineProfile(line_idx_LyaE);

    if (profile->isAsym())
      setAsymProfile(elt_idx_LyaE, line_idx_LyaE, redshift, catalog);
  }
  // SymIgm Profile
  {
    for (auto const &[elt_idx_igm, line_indices_igm] : indices_Igm) {
      auto const &param_EltIgm =
          m_elementsVector->getElementParam()[elt_idx_igm];

      // const auto &elt = getElementList()[elt_idx_igm];
      auto line_indices_filtered = line_indices_igm;
      auto end = std::remove_if(
          line_indices_filtered.begin(), line_indices_filtered.end(),
          [&param_EltIgm](Int32 idx) {
            return !param_EltIgm->getLineProfile(idx)->isSymIgm();
          });
      line_indices_filtered.erase(end, line_indices_filtered.end());
      if (!line_indices_filtered.empty())
        setSymIgmProfile(elt_idx_igm, line_indices_filtered, redshift);
    }
  }
}

void CLineRatioManager::setAsymProfile(Int32 idxLyaE, Int32 idxLineLyaE,
                                       Float64 redshift,
                                       const CLineMap &catalog) {
  auto const &param_LyaE = m_elementsVector->getElementParam()[idxLyaE];
  Int32 lineId = param_LyaE->m_Lines[idxLineLyaE].GetID();
  auto const &ref_line = catalog.at(lineId);

  // finding or setting the correct profile
  auto const &ref_profile = ref_line.GetProfile();
  CLineProfile_ptr profile;
  if (m_forceDisableLyaFitting && ref_profile->isAsymFit())
    // convert asymfit to asymfixed profile
    profile = dynamic_cast<const CLineProfileASYMFIT *>(ref_profile.get())
                  ->cloneToASYM();
  else if (m_forceLyaFitting && ref_profile->isAsymFixed())
    // convert asymfixed to asymfit
    profile = dynamic_cast<const CLineProfileASYM *>(ref_profile.get())
                  ->cloneToASYMFIT();
  else
    profile = ref_profile->Clone();

  param_LyaE->SetLineProfile(idxLineLyaE, std::move(profile));
}

void CLineRatioManager::setSymIgmProfile(Int32 iElts,
                                         const TInt32List &idxLineIGM,
                                         Float64 redshift) {

  auto const &param_EltIgm = m_elementsVector->getElementParam()[iElts];

  bool fixedIGM =
      (m_continuumManager->isContinuumComponentFitter()) &&
      !m_continuumManager->isContFittedToNull(); // set to false when continuum
                                                 // is fitted to null

  if (fixedIGM) {
    Int32 igmidx = m_continuumManager->getFittedMeiksinIndex();
    param_EltIgm->SetSymIgmFit(false);
    param_EltIgm->SetSymIgmParams(TSymIgmParams(igmidx, redshift));
  } else {
    param_EltIgm->SetSymIgmFit(true);
  }
}

bool CLineRatioManager::init(Float64 redshift, Int32 itratio) {
  // prepare the Lya width and asym coefficients if the asymfit profile
  // option is met
  setLyaProfile(redshift, m_RestLineList);

  return false;
}

/**
 * @brief setPassMode
 * @param iPass
 * set the fitting parameters according the the iPass argument.
 * @return
 */
void CLineRatioManager::setPassMode(Int32 iPass) {

  if (iPass == 1) {
    m_forceDisableLyaFitting = true;
    m_forceLyaFitting = false;
  }
  if (iPass == 2) {
    m_forceDisableLyaFitting = m_opt_lya_forcedisablefit;
    m_forceLyaFitting = m_opt_lya_forcefit;
    Log.LogDetail(
        Formatter()
        << "    model: set forceLyaFitting ASYMFIT for Tpl-ratio mode : "
        << m_forceLyaFitting);
  }
  if (iPass == 3) {
    m_forceDisableLyaFitting = false;
  }
}

/**
 * \brief Accumulates the squared differences between model and spectrum in
 *the argument lambdaRange and returns the sum.
 **/
Float64 CLineRatioManager::getLeastSquareMerit() const {
  Float64 fit = 0.0;

  for (auto &spcIndex : m_spectraIndex) {

    const CSpectrumSpectralAxis &spcSpectralAxis =
        getSpectrum().GetSpectralAxis();
    const auto &ErrorNoContinuum = getSpectrum().GetErrorAxis();

    const CSpectrumFluxAxis &Yspc = getModel().getSpcFluxAxis();
    const CSpectrumFluxAxis &Ymodel =
        getModel().GetModelSpectrum().GetFluxAxis();

    Float64 diff = 0.0;

    Int32 imin =
        spcSpectralAxis.GetIndexAtWaveLength(getLambdaRange().GetBegin());
    Int32 imax =
        spcSpectralAxis.GetIndexAtWaveLength(getLambdaRange().GetEnd());
    for (Int32 j = imin; j < imax; j++) {
      diff = (Yspc[j] - Ymodel[j]);
      fit += (diff * diff) / (ErrorNoContinuum[j] * ErrorNoContinuum[j]);
    }
    if (std::isnan(fit)) {
      Log.LogDetail(
          Formatter()
          << "CLineModelFitting::getLeastSquareMerit: NaN value found on "
             "the lambdarange = ("
          << getLambdaRange().GetBegin() << ", " << getLambdaRange().GetEnd()
          << ")");
      Log.LogDetail(
          Formatter()
          << "CLineModelFitting::getLeastSquareMerit: NaN value found on "
             "the true observed spectral axis lambdarange = ("
          << spcSpectralAxis[imin] << ", " << spcSpectralAxis[imax] << ")");
      for (Int32 j = imin; j < imax; j++) {
        if (std::isnan(Yspc[j])) {
          Log.LogDetail(
              Formatter()
              << "CLineModelFitting::getLeastSquareMerit: NaN value found "
                 "for the observed spectrum at lambda="
              << spcSpectralAxis[j]);
          break;
        }
        if (std::isnan(Ymodel[j])) {
          Log.LogDetail(
              Formatter()
              << "CLineModelFitting::getLeastSquareMerit: NaN value found "
                 "for the model at lambda="
              << spcSpectralAxis[j]);
          break;
        }

        if (std::isnan(ErrorNoContinuum[j])) {
          Log.LogDetail(
              Formatter()
              << "CLineModelFitting::getLeastSquareMerit: NaN value found "
                 "for the sqrt(variance) at lambda="
              << spcSpectralAxis[j]);
          break;
        }
        if (ErrorNoContinuum[j] == 0.0) {
          Log.LogDetail(
              Formatter()
              << "CLineModelFitting::getLeastSquareMerit: 0 value found "
                 "for the sqrt(variance) at lambda="
              << spcSpectralAxis[j]);
          break;
        }
      }
      THROWG(ErrorCode::INTERNAL_ERROR, "computed fit is NaN");
    }
  }
  fit += m_continuumManager->getFitSum();

  return fit;
}

void CLineRatioManager::logParameters() {

  Log.LogDetail(Formatter() << " m_opt_lya_forcefit" << m_opt_lya_forcefit);
  Log.LogDetail(Formatter()
                << " m_opt_lya_forcedisablefit" << m_opt_lya_forcedisablefit);
}

std::shared_ptr<CLineRatioManager> CLineRatioManager::makeLineRatioManager(
    const std::string &lineRatioType,
    const std::shared_ptr<CLMEltListVector> &elementsVector,
    const CSpcModelVectorPtr &models, const CCSpectrumVectorPtr &inputSpcs,
    const CTLambdaRangePtrVector &lambdaRanges,
    std::shared_ptr<CContinuumManager> continuumManager,
    const CLineMap &restLineList, std::shared_ptr<CAbstractFitter> fitter,
    const CSpectraGlobalIndex &spcIndex) {
  std::shared_ptr<CLineRatioManager> ret;
  if (lineRatioType == "tplRatio")
    ret = std::make_shared<CTplratioManager>(
        CTplratioManager(elementsVector, models, inputSpcs, lambdaRanges,
                         continuumManager, restLineList, spcIndex));
  else if (lineRatioType == "tplCorr")
    ret = std::make_shared<CTplCorrManager>(
        CTplCorrManager(elementsVector, models, inputSpcs, lambdaRanges,
                        continuumManager, restLineList, spcIndex));
  else if (lineRatioType == "rules")
    ret = std::make_shared<CRulesManager>(
        CRulesManager(elementsVector, models, inputSpcs, lambdaRanges,
                      continuumManager, restLineList, spcIndex));
  else
    THROWG(ErrorCode::INVALID_PARAMETER,
           "Only {tplratio, rules, tpcorr} values are "
           "supported for linemodel.lineRatioType");
  ret->setFitter(fitter);

  return ret;
}

void CLineRatioManager::refreshAllModels() {
  for (auto &spcIndex : m_spectraIndex) {
    getModel().refreshModel();
  }
}
