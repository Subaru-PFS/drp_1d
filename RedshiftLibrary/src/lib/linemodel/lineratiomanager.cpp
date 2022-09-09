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

  m_ContinuumComponent = ps->GetScoped<std::string>("continuumcomponent");

  if (Context.GetCurrentMethod() == "LineModelSolve") {

    m_opt_lya_forcefit = ps->GetScoped<bool>("lyaforcefit");
    m_opt_lya_forcedisablefit = ps->GetScoped<bool>("lyaforcedisablefit");

    m_opt_lya_fit_asym_min = ps->GetScoped<Float64>("lyafit.asymfitmin");
    m_opt_lya_fit_asym_max = ps->GetScoped<Float64>("lyafit.asymfitmax");
    m_opt_lya_fit_asym_step = ps->GetScoped<Float64>("lyafit.asymfitstep");
    m_opt_lya_fit_width_min = ps->GetScoped<Float64>("lyafit.widthfitmin");
    m_opt_lya_fit_width_max = ps->GetScoped<Float64>("lyafit.widthfitmax");
    m_opt_lya_fit_width_step = ps->GetScoped<Float64>("lyafit.widthfitstep");
    m_opt_lya_fit_delta_min = ps->GetScoped<Float64>("lyafit.deltafitmin");
    m_opt_lya_fit_delta_max = ps->GetScoped<Float64>("lyafit.deltafitmax");
    m_opt_lya_fit_delta_step = ps->GetScoped<Float64>("lyafit.deltafitstep");
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
        m_Elements[idxLyaE]->m_Lines[idxLineLyaE].GetProfile();
    if (profile.isAsym())
      setAsymProfile(idxLyaE, idxLineLyaE, redshift, catalog);
  }

  for (Int32 i = 0; i < idxEltIGM.size(); ++i) {
    const auto &Elt = m_Elements[idxEltIGM[i]];
    if (!Elt->IsOutsideLambdaRange()) {
      TInt32List &idxLine = idxLineIGM[i];
      auto end =
          std::remove_if(idxLine.begin(), idxLine.end(), [Elt](Int32 idx) {
            return !Elt->m_Lines[idx].GetProfile().isSymIgm();
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

  m_Elements[idxLyaE]->m_Lines[idxLineLyaE].SetProfile(std::move(profile));

  if (!doasymfit)
    return;

  // find the best width and asym coeff. parameters
  TAsymParams bestfitParams = fitAsymParameters(redshift, idxLyaE, idxLineLyaE);

  // set the associated Lya members in the element definition
  m_Elements[idxLyaE]->SetAsymfitParams(bestfitParams);
}

TAsymParams CLineRatioManager::fitAsymParameters(Float64 redshift,
                                                 Int32 idxLyaE,
                                                 const Int32 &idxLineLyaE) {

  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc->GetSpectralAxis();
  // 3. find the best width and asym coeff. parameters
  Float64 widthCoeffStep = m_opt_lya_fit_width_step;
  Float64 widthCoeffMin = m_opt_lya_fit_width_min;
  Float64 widthCoeffMax = m_opt_lya_fit_width_max;
  Int32 nWidthSteps =
      int((widthCoeffMax - widthCoeffMin) / widthCoeffStep + 1.5);
  Float64 asymCoeffStep = m_opt_lya_fit_asym_step;
  Float64 asymCoeffMin = m_opt_lya_fit_asym_min;
  Float64 asymCoeffMax = m_opt_lya_fit_asym_max;
  Int32 nAsymSteps = int((asymCoeffMax - asymCoeffMin) / asymCoeffStep + 1.5);
  Float64 deltaStep = m_opt_lya_fit_delta_step;
  Float64 deltaMin = m_opt_lya_fit_delta_min;
  Float64 deltaMax = m_opt_lya_fit_delta_max;
  Int32 nDeltaSteps = int((deltaMax - deltaMin) / deltaStep + 1.5);

  TAsymParams bestparams = {widthCoeffMin, asymCoeffMin, deltaMin};
  Float64 meritMin = DBL_MAX;

  TInt32List filterEltsIdxLya(1, idxLyaE);

  for (Int32 iDelta = 0; iDelta < nDeltaSteps; iDelta++) {
    Float64 delta = deltaMin + deltaStep * iDelta;
    for (Int32 iWidth = 0; iWidth < nWidthSteps; iWidth++) {
      Float64 asymWidthCoeff = widthCoeffMin + widthCoeffStep * iWidth;
      for (Int32 iAsym = 0; iAsym < nAsymSteps; iAsym++) {
        Float64 asymAlphaCoeff = asymCoeffMin + asymCoeffStep * iAsym;
        m_Elements[idxLyaE]->SetAsymfitParams(
            {asymWidthCoeff, asymAlphaCoeff, delta});

        // idxLineLyaE = -1;
        m_Elements[idxLyaE]->fitAmplitude(
            spectralAxis, m_model->getSpcFluxAxisNoContinuum(),
            m_model->getContinuumFluxAxis(), redshift, idxLineLyaE);

        Float64 m = 0; // TODO DV why initializing to m_dTransposeD ?;
        if (1) {

          m_model->refreshModelUnderElements(filterEltsIdxLya, idxLineLyaE);
          m = m_model->getModelErrorUnderElement(idxLyaE,
                                                 m_model->getSpcFluxAxis());
        } else {
          m = getLeastSquareMeritFast(idxLyaE);
        }
        if (m < meritMin) {
          meritMin = m;
          bestparams = m_Elements[idxLyaE]->GetAsymfitParams(0);
        }

        Log.LogDebug("Fitting Lya Profile: width=%f, asym=%f, delta=%f",
                     asymWidthCoeff, asymAlphaCoeff, delta);
        Log.LogDebug("Fitting Lya Profile: merit=%e", m);
        Log.LogDebug("Fitting Lya Profile: idxLyaE=%d, idxLineLyaE=%d", idxLyaE,
                     idxLineLyaE);
      }
    }
  }
  Log.LogDebug("Lya Profile found: width=%f, asym=%f, delta=%f",
               bestparams.sigma, bestparams.alpha, bestparams.delta);
  return bestparams;
}

Int32 CLineRatioManager::getLineIndexInCatalog(
    Int32 iElts, Int32 idxLine,
    const CLineCatalog::TLineVector &catalog) const {
  Int32 lineIndex = undefIdx;
  lineIndex = m_Elements[iElts]->m_LineCatalogIndexes[idxLine];
  if (lineIndex < 0 || lineIndex >= catalog.size())
    THROWG(INTERNAL_ERROR, "Lya idx out-of-bound");

  return lineIndex;
}

void CLineRatioManager::setSymIgmProfile(Int32 iElts,
                                         const TInt32List &idxLineIGM,
                                         Float64 redshift) {

  bool fixedIGM = isContinuumComponentTplfitxx();

  // set to false when continuum is fitted to null
  fixedIGM &= m_continuumManager->isContFittedToNull();

  Int32 bestigmidx = fixedIGM
                         ? m_continuumManager->getFittedMeiksinIndex()
                         : fitAsymIGMCorrection(redshift, iElts, idxLineIGM);

  m_Elements[iElts]->SetSymIgmParams(TSymIgmParams(bestigmidx, redshift));
}

Int32 CLineRatioManager::fitAsymIGMCorrection(Float64 redshift, Int32 iElts,
                                              const TInt32List &idxLine) {

  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc->GetSpectralAxis();
  if (spectralAxis[0] / (1 + redshift) > RESTLAMBDA_LYA)
    return -1;

  Float64 meritMin = DBL_MAX;
  Int32 bestIgmIdx = -1;

  Int32 igmCount =
      m_Elements[iElts]->getLineProfile(idxLine.front()).getIGMIdxCount();
  for (Int32 igmIdx = 0; igmIdx < igmCount; igmIdx++) {
    m_Elements[iElts]->SetSymIgmParams(TSymIgmParams(igmIdx, redshift));
    m_Elements[iElts]->fitAmplitude(spectralAxis,
                                    m_model->getSpcFluxAxisNoContinuum(),
                                    m_model->getContinuumFluxAxis(), redshift);

    m_model->refreshModelUnderElements(TInt32List(1, iElts));
    Float64 m =
        m_model->getModelErrorUnderElement(iElts, m_model->getSpcFluxAxis());

    if (m < meritMin) {
      meritMin = m;
      bestIgmIdx = igmIdx;
    }
  }
  return bestIgmIdx;
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
  const CSpectrumFluxAxis &spcFluxAxis = m_model->getSpcFluxAxis();
  const CSpectrumFluxAxis &modelFluxAxis =
      m_model->GetModelSpectrum().GetFluxAxis();

  // Int32 numDevs = 0;
  Float64 fit = 0.0;
  const Float64 *Ymodel = modelFluxAxis.GetSamples();
  const Float64 *Yspc = spcFluxAxis.GetSamples();
  const auto &ErrorNoContinuum = m_inputSpc->GetFluxAxis().GetError();
  Float64 diff = 0.0;

  Int32 imin = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange->GetBegin());
  Int32 imax = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange->GetEnd());
  for (Int32 j = imin; j < imax; j++) {
    // numDevs++;
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

  if (isContinuumComponentTplfitxx()) {
    fit += m_continuumManager->getFitSum(); // unconditionnal sum (if photometry
                                            // disabled, will sum 0.0)
  }

  //  Log.LogDebug("CLineModelFitting::getLeastSquareMerit fit = %f", fit);
  if (std::isnan(fit)) {
    Log.LogError("CLineModelFitting::getLeastSquareMerit: NaN value found on "
                 "the lambdarange = (%f, %f)",
                 m_lambdaRange->GetBegin(), m_lambdaRange->GetEnd());
    Log.LogError("CLineModelFitting::getLeastSquareMerit: NaN value found on "
                 "the true observed spectral axis lambdarange = (%f, %f)",
                 spcSpectralAxis[imin], spcSpectralAxis[imax]);
    for (Int32 j = imin; j < imax; j++) {
      if (std::isnan(Yspc[j])) {
        Log.LogError("CLineModelFitting::getLeastSquareMerit: NaN value found "
                     "for the observed spectrum at lambda=%f",
                     spcSpectralAxis[j]);
        break;
      }
      if (std::isnan(Ymodel[j])) {
        Log.LogError("CLineModelFitting::getLeastSquareMerit: NaN value found "
                     "for the model at lambda=%f",
                     spcSpectralAxis[j]);
        break;
      }

      if (std::isnan(ErrorNoContinuum[j])) {
        Log.LogError("CLineModelFitting::getLeastSquareMerit: NaN value found "
                     "for the sqrt(variance) at lambda=%f",
                     spcSpectralAxis[j]);
        break;
      }
      if (ErrorNoContinuum[j] == 0.0) {
        Log.LogError("CLineModelFitting::getLeastSquareMerit: 0 value found "
                     "for the sqrt(variance) at lambda=%f",
                     spcSpectralAxis[j]);
        break;
      }
    }
    THROWG(INTERNAL_ERROR, "NaN value found");
  }
  return fit;
}

/**
 * \brief Get the squared difference by fast method proposed by D. Vibert
 **/
Float64 CLineRatioManager::getLeastSquareMeritFast(Int32 idxLine) const {
  Float64 fit = -1; // TODO restore getLeastSquareContinuumMeritFast();

  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    if (idxLine != undefIdx && idxLine != iElts) {
      continue;
    }
    Float64 dtm = m_Elements[iElts]->GetSumCross();
    Float64 mtm = m_Elements[iElts]->GetSumGauss();
    Float64 a = m_Elements[iElts]->GetFitAmplitude();
    Float64 term1 = a * a * mtm;
    Float64 term2 = -2. * a * dtm;
    fit += term1 + term2;
  }

  Log.LogDebug("CLineModelFitting::getLeastSquareMerit fit fast = %f", fit);
  return fit;
}

void CLineRatioManager::logParameters() {

  Log.LogDetail(Formatter() << " m_opt_lya_forcefit" << m_opt_lya_forcefit);
  Log.LogDetail(Formatter()
                << " m_opt_lya_forcedisablefit" << m_opt_lya_forcedisablefit);
  Log.LogDetail(Formatter()
                << " m_opt_lya_fit_asym_min" << m_opt_lya_fit_asym_min);
  Log.LogDetail(Formatter()
                << " m_opt_lya_fit_asym_max" << m_opt_lya_fit_asym_max);
  Log.LogDetail(Formatter()
                << " m_opt_lya_fit_asym_step" << m_opt_lya_fit_asym_step);
  Log.LogDetail(Formatter()
                << " m_opt_lya_fit_width_min" << m_opt_lya_fit_width_min);
  Log.LogDetail(Formatter()
                << " m_opt_lya_fit_width_max" << m_opt_lya_fit_width_max);
  Log.LogDetail(Formatter()
                << " m_opt_lya_fit_width_step" << m_opt_lya_fit_width_step);
  Log.LogDetail(Formatter()
                << " m_opt_lya_fit_delta_min" << m_opt_lya_fit_delta_min);
  Log.LogDetail(Formatter()
                << " m_opt_lya_fit_delta_max" << m_opt_lya_fit_delta_max);
  Log.LogDetail(Formatter()
                << " m_opt_lya_fit_delta_step" << m_opt_lya_fit_delta_step);
}

std::shared_ptr<CLineRatioManager> CLineRatioManager::makeLineRatioManager(
    const std::string &lineRatioType, CLineModelElementList &elements,
    std::shared_ptr<CSpectrumModel> model,
    std::shared_ptr<const CSpectrum> inputSpc,
    std::shared_ptr<const TFloat64Range> lambdaRange,
    std::shared_ptr<CContinuumManager> continuumManager,
    const CLineCatalog::TLineVector &restLineList)

{

  if (lineRatioType == "tplratio")
    return std::make_shared<CTplratioManager>(
        CTplratioManager(elements, model, inputSpc, lambdaRange,
                         continuumManager, restLineList));
  else if (lineRatioType == "tplcorr")
    return std::make_shared<CTplCorrManager>(
        CTplCorrManager(elements, model, inputSpc, lambdaRange,
                        continuumManager, restLineList));
  else if (lineRatioType == "rules")
    return std::make_shared<CRulesManager>(
        CRulesManager(elements, model, inputSpc, lambdaRange, continuumManager,
                      restLineList));
  else
    THROWG(INVALID_PARAMETER, "Only {tplratio, rules, tpcorr} values are "
                              "supported for linemodel.lineRatioType");
}
