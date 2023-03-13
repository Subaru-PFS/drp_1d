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
    CLineModelElementList &elements, std::shared_ptr<CSpectrumModel> model,
    std::shared_ptr<const CSpectrum> inputSpc,
    std::shared_ptr<const TFloat64Range> lambdaRange,
    std::shared_ptr<CContinuumManager> continuumManager,
    const CLineCatalog::TLineVector &restLineList)
    : m_Elements(elements), m_model(model), m_inputSpc(inputSpc),
      m_lambdaRange(lambdaRange), m_continuumManager(continuumManager),
      m_RestLineList(restLineList) {

  CAutoScope autoscope(Context.m_ScopeStack, "linemodel");
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();

  if (Context.GetCurrentMethod() == "LineModelSolve") {
    m_opt_lya_forcefit = ps->GetScoped<bool>("lyaforcefit");
    m_opt_lya_forcedisablefit = ps->GetScoped<bool>("lyaforcedisablefit");
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

void CLineRatioManager::setLyaProfile(
    Float64 redshift, const CLineCatalog::TLineVector &catalog) {
  auto idxLineIGM_ = m_Elements.getIgmLinesIndices();
  auto const idxEltIGM = std::move(idxLineIGM_.front());
  std::vector<TInt32List> idxLineIGM(
      std::make_move_iterator(idxLineIGM_.begin() + 1),
      std::make_move_iterator(idxLineIGM_.end()));

  if (idxEltIGM.empty())
    return;

  // assuming only one asymfit/fixed profile
  Int32 idxLyaE = idxEltIGM.front();
  Int32 idxLineLyaE = idxLineIGM.front().front();

  if (!m_Elements[idxLyaE]->IsOutsideLambdaRange(idxLineLyaE)) {
    const auto &profile =
        m_Elements[idxLyaE]->GetLines()[idxLineLyaE].GetProfile();
    if (profile.isAsym())
      setAsymProfile(idxLyaE, idxLineLyaE, redshift, catalog);
  }

  for (Int32 i = 0; i < idxEltIGM.size(); ++i) {
    const auto &Elt = m_Elements[idxEltIGM[i]];
    if (!Elt->IsOutsideLambdaRange()) {
      TInt32List &idxLine = idxLineIGM[i];
      auto end =
          std::remove_if(idxLine.begin(), idxLine.end(), [Elt](Int32 idx) {
            return !Elt->GetLines()[idx].GetProfile().isSymIgm();
          });
      idxLine.erase(end, idxLine.end());
      if (!idxLine.empty())
        setSymIgmProfile(idxEltIGM[i], idxLine, redshift);
    }
  }
}

void CLineRatioManager::setAsymProfile(
    Int32 idxLyaE, Int32 idxLineLyaE, Float64 redshift,
    const CLineCatalog::TLineVector &catalog) {
  Int32 lineIndex = getLineIndexInCatalog(idxLyaE, idxLineLyaE, catalog);
  if (lineIndex == undefIdx)
    return;

  // finding or setting the correct profile
  CLineProfile_ptr profile;
  if (m_forceDisableLyaFitting && catalog[lineIndex].GetProfile().isAsymFit())
    // convert asymfit to asymfixed profile
    profile = dynamic_cast<const CLineProfileASYMFIT *>(
                  &catalog[lineIndex].GetProfile())
                  ->cloneToASYM();
  else if (m_forceLyaFitting && catalog[lineIndex].GetProfile().isAsymFixed())
    // convert asymfixed to asymfit
    profile =
        dynamic_cast<const CLineProfileASYM *>(&catalog[lineIndex].GetProfile())
            ->cloneToASYMFIT();
  else
    profile = catalog[lineIndex].GetProfile().Clone();

  bool doasymfit = profile->isAsymFit();

  m_Elements[idxLyaE]->SetLineProfile(idxLineLyaE, std::move(profile));

  if (!doasymfit)
    return;

  // find the best width and asym coeff. parameters
  TAsymParams bestfitParams =
      m_fitter->fitAsymParameters(redshift, idxLyaE, idxLineLyaE);

  // set the associated Lya members in the element definition
  m_Elements[idxLyaE]->SetAsymfitParams(bestfitParams);
}

Int32 CLineRatioManager::getLineIndexInCatalog(
    Int32 iElts, Int32 idxLine,
    const CLineCatalog::TLineVector &catalog) const {
  return m_Elements[iElts]->getLineIndexInCatalog(idxLine, catalog);
}

void CLineRatioManager::setSymIgmProfile(Int32 iElts,
                                         const TInt32List &idxLineIGM,
                                         Float64 redshift) {

  bool fixedIGM = m_continuumManager->isContinuumComponentTplfitxx();

  // set to false when continuum is fitted to null
  fixedIGM &= m_continuumManager->isContFittedToNull();

  Int32 bestigmidx =
      fixedIGM ? m_continuumManager->getFittedMeiksinIndex()
               : m_fitter->fitAsymIGMCorrection(redshift, iElts, idxLineIGM);

  m_Elements[iElts]->SetSymIgmParams(TSymIgmParams(bestigmidx, redshift));
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
        "    model: set forceLyaFitting ASYMFIT for Tpl-ratio mode : %d",
        m_forceLyaFitting);
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
  const CSpectrumSpectralAxis &spcSpectralAxis = m_inputSpc->GetSpectralAxis();
  const auto &ErrorNoContinuum = m_inputSpc->GetErrorAxis();

  const CSpectrumFluxAxis &Yspc = m_model->getSpcFluxAxis();
  const CSpectrumFluxAxis &Ymodel = m_model->GetModelSpectrum().GetFluxAxis();

  Float64 fit = 0.0;
  Float64 diff = 0.0;

  Int32 imin = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange->GetBegin());
  Int32 imax = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange->GetEnd());
  for (Int32 j = imin; j < imax; j++) {
    diff = (Yspc[j] - Ymodel[j]);
    fit += (diff * diff) / (ErrorNoContinuum[j] * ErrorNoContinuum[j]);
    //        if ( 1E6 * diff < ErrorNoContinuum[j] )
    //        {
    //            Log.LogDebug( "Warning: noise is at least 6 orders greater
    //            than the residue!" ); Log.LogDebug(
    //            "CLineModelFitting::getLeastSquareMerit diff = %f", diff );
    //            Log.LogDebug( "CLineModelFitting::getLeastSquareMerit
    //            ErrorNoContinuum[%d] = %f", j, ErrorNoContinuum[j] );
    //        }
  }

  fit += m_continuumManager->getFitSum();

  if (std::isnan(fit)) {
    Log.LogDetail("CLineModelFitting::getLeastSquareMerit: NaN value found on "
                  "the lambdarange = (%f, %f)",
                  m_lambdaRange->GetBegin(), m_lambdaRange->GetEnd());
    Log.LogDetail("CLineModelFitting::getLeastSquareMerit: NaN value found on "
                  "the true observed spectral axis lambdarange = (%f, %f)",
                  spcSpectralAxis[imin], spcSpectralAxis[imax]);
    for (Int32 j = imin; j < imax; j++) {
      if (std::isnan(Yspc[j])) {
        Log.LogDetail("CLineModelFitting::getLeastSquareMerit: NaN value found "
                      "for the observed spectrum at lambda=%f",
                      spcSpectralAxis[j]);
        break;
      }
      if (std::isnan(Ymodel[j])) {
        Log.LogDetail("CLineModelFitting::getLeastSquareMerit: NaN value found "
                      "for the model at lambda=%f",
                      spcSpectralAxis[j]);
        break;
      }

      if (std::isnan(ErrorNoContinuum[j])) {
        Log.LogDetail("CLineModelFitting::getLeastSquareMerit: NaN value found "
                      "for the sqrt(variance) at lambda=%f",
                      spcSpectralAxis[j]);
        break;
      }
      if (ErrorNoContinuum[j] == 0.0) {
        Log.LogDetail("CLineModelFitting::getLeastSquareMerit: 0 value found "
                      "for the sqrt(variance) at lambda=%f",
                      spcSpectralAxis[j]);
        break;
      }
    }
    THROWG(INTERNAL_ERROR, "computed fit is NaN");
  }
  return fit;
}

void CLineRatioManager::logParameters() {

  Log.LogDetail(Formatter() << " m_opt_lya_forcefit" << m_opt_lya_forcefit);
  Log.LogDetail(Formatter()
                << " m_opt_lya_forcedisablefit" << m_opt_lya_forcedisablefit);
}

std::shared_ptr<CLineRatioManager> CLineRatioManager::makeLineRatioManager(
    const std::string &lineRatioType, CLineModelElementList &elements,
    std::shared_ptr<CSpectrumModel> model,
    std::shared_ptr<const CSpectrum> inputSpc,
    std::shared_ptr<const TFloat64Range> lambdaRange,
    std::shared_ptr<CContinuumManager> continuumManager,
    const CLineCatalog::TLineVector &restLineList,
    std::shared_ptr<CAbstractFitter> fitter) {
  std::shared_ptr<CLineRatioManager> ret;
  if (lineRatioType == "tplratio")
    ret = std::make_shared<CTplratioManager>(
        CTplratioManager(elements, model, inputSpc, lambdaRange,
                         continuumManager, restLineList));
  else if (lineRatioType == "tplcorr")
    ret = std::make_shared<CTplCorrManager>(
        CTplCorrManager(elements, model, inputSpc, lambdaRange,
                        continuumManager, restLineList));
  else if (lineRatioType == "rules")
    ret = std::make_shared<CRulesManager>(
        CRulesManager(elements, model, inputSpc, lambdaRange, continuumManager,
                      restLineList));
  else
    THROWG(INVALID_PARAMETER, "Only {tplratio, rules, tpcorr} values are "
                              "supported for linemodel.lineRatioType");
  ret->setFitter(fitter);
  return ret;
}
