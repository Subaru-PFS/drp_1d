#include "RedshiftLibrary/linemodel/abstractfitter.h"
#include "RedshiftLibrary/linemodel/hybridfitter.h"
#include "RedshiftLibrary/linemodel/individualfitter.h"
#include "RedshiftLibrary/linemodel/onesfitter.h"
#include "RedshiftLibrary/linemodel/randomfitter.h"
#include "RedshiftLibrary/linemodel/svdfitter.h"
#include "RedshiftLibrary/linemodel/svdlcfitter.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;

CAbstractFitter::CAbstractFitter(
    CLineModelElementList &elements,
    std::shared_ptr<const CSpectrum> inputSpectrum,
    std::shared_ptr<const TLambdaRange> lambdaRange,
    std::shared_ptr<CSpectrumModel> spectrumModel,
    const TLineVector &restLineList,
    const std::vector<TLineModelElementParam_ptr> &elementParam)
    : m_Elements(elements), m_inputSpc(*(inputSpectrum)),
      m_RestLineList(restLineList), m_lambdaRange(*(lambdaRange)),
      m_model(spectrumModel), m_ElementParam(elementParam) {

  CAutoScope autoscope(Context.m_ScopeStack, "linemodel");
  if (Context.GetCurrentMethod() == "LineModelSolve") {
    std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();
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

std::shared_ptr<CAbstractFitter> CAbstractFitter::makeFitter(
    std::string fittingMethod, CLineModelElementList &elements,
    std::shared_ptr<const CSpectrum> inputSpectrum,
    std::shared_ptr<const TLambdaRange> lambdaRange,
    std::shared_ptr<CSpectrumModel> spectrumModel,
    const TLineVector &restLineList,
    std::shared_ptr<CContinuumManager> continuumManager,
    const std::vector<TLineModelElementParam_ptr> &elementParam) {
  if (fittingMethod == "hybrid")
    return std::make_shared<CHybridFitter>(
        CHybridFitter(elements, inputSpectrum, lambdaRange, spectrumModel,
                      restLineList, elementParam));
  else if (fittingMethod == "svd")
    return std::make_shared<CSvdFitter>(CSvdFitter(elements, inputSpectrum,
                                                   lambdaRange, spectrumModel,
                                                   restLineList, elementParam));
  else if (fittingMethod == "svdlc")
    return std::make_shared<CSvdlcFitter>(
        CSvdlcFitter(elements, inputSpectrum, lambdaRange, spectrumModel,
                     restLineList, elementParam, continuumManager));
  else if (fittingMethod == "svdlcp2")
    return std::make_shared<CSvdlcFitter>(
        CSvdlcFitter(elements, inputSpectrum, lambdaRange, spectrumModel,
                     restLineList, elementParam, continuumManager, 2));

  else if (fittingMethod == "ones")
    return std::make_shared<COnesFitter>(
        COnesFitter(elements, inputSpectrum, lambdaRange, spectrumModel,
                    restLineList, elementParam));
  else if (fittingMethod == "random")
    return std::make_shared<CRandomFitter>(
        CRandomFitter(elements, inputSpectrum, lambdaRange, spectrumModel,
                      restLineList, elementParam));
  else if (fittingMethod == "individual")
    return std::make_shared<CIndividualFitter>(
        CIndividualFitter(elements, inputSpectrum, lambdaRange, spectrumModel,
                          restLineList, elementParam));
  else
    THROWG(INTERNAL_ERROR, Formatter()
                               << "Unknown fitting method " << fittingMethod);
}

void CAbstractFitter::logParameters() {
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

void CAbstractFitter::computeCrossProducts(
    CLineModelElement &elt, const CSpectrumSpectralAxis &spectralAxis,
    const CSpectrumFluxAxis &noContinuumfluxAxis,
    const CSpectrumFluxAxis &continuumfluxAxis, Float64 redshift,
    Int32 lineIdx) {

  Float64 sumCross = 0.0;
  Float64 sumGauss = 0.0;
  Float64 dtmFree = 0.0;

  const CSpectrumNoiseAxis &error = noContinuumfluxAxis.GetError();

  Float64 y = 0.0;
  Float64 x = 0.0;
  Float64 yg = 0.0;
  Float64 c = 1.0;

  Float64 err2 = 0.0;
  Int32 num = 0;

  Int32 nLines = elt.GetSize();
  for (Int32 k = 0; k < nLines; k++) { // loop for the intervals
    if (elt.IsOutsideLambdaRange(k)) {
      continue;
    }
    if (lineIdx != undefIdx && !elt.isLineActiveOnSupport(k, lineIdx)) {
      continue;
    }

    for (Int32 i = elt.getStartNoOverlap(k); i <= elt.getEndNoOverlap(k); i++) {
      c = continuumfluxAxis[i];
      y = noContinuumfluxAxis[i];
      x = spectralAxis[i];

      yg = 0.0;

      for (Int32 k2 = 0; k2 < nLines; k2++) { // loop for the signal synthesis
        if (elt.IsOutsideLambdaRange(k2) || !elt.isLineActiveOnSupport(k2, k)) {
          continue;
        }
        Int32 sf = elt.getSignFactor(k2);
        if (sf == -1) {
          yg += sf * c * elt.GetNominalAmplitude(k2) *
                elt.GetLineProfileAtRedshift(k2, redshift, x);
        } else {
          yg += sf * elt.GetNominalAmplitude(k2) *
                elt.GetLineProfileAtRedshift(k2, redshift, x);
        }
      }
      num++;
      err2 = 1.0 / (error[i] * error[i]);
      dtmFree += yg * y * err2;
      sumGauss += yg * yg * err2;
    }
  }

  if (num == 0 || sumGauss == 0) {
    Log.LogDebug("CLineModelElement::fitAmplitude: Could not fit amplitude:    "
                 " num=%d, mtm=%f",
                 num, sumGauss);
    sumGauss = NAN;
    dtmFree = NAN;
  }

  elt.SetSumGauss(sumGauss);
  elt.SetDtmFree(dtmFree);
}

void CAbstractFitter::fitAmplitude(Int32 eltIndex,
                                   const CSpectrumSpectralAxis &spectralAxis,
                                   const CSpectrumFluxAxis &noContinuumfluxAxis,
                                   const CSpectrumFluxAxis &continuumfluxAxis,
                                   Float64 redshift, Int32 lineIdx) {
  auto &elt = m_Elements[eltIndex];
  Int32 nLines = elt->GetSize();
  auto &elementParam = m_ElementParam[eltIndex];

  elementParam->m_FittedAmplitudes.assign(nLines, NAN);
  elementParam->m_FittedAmplitudeErrorSigmas.assign(nLines, NAN);

  if (elt->IsOutsideLambdaRange()) {
    elementParam->m_fitAmplitude = NAN;
    elt->SetSumCross(NAN);
    elt->SetSumGauss(NAN);
    elt->SetDtmFree(NAN);
    return;
  }

  computeCrossProducts(*elt, spectralAxis, noContinuumfluxAxis,
                       continuumfluxAxis, redshift, lineIdx);

  if (std::isnan(elt->GetSumGauss())) {
    elementParam->m_fitAmplitude = NAN;
    elt->SetSumCross(NAN);
    return;
  }

  elt->SetSumCross(std::max(0.0, elt->GetDtmFree()));
  Float64 A = elt->GetSumCross() / elt->GetSumGauss();
  elementParam->m_fitAmplitude = A;

  for (Int32 k = 0; k < nLines; k++) {
    if (elt->IsOutsideLambdaRange(k)) {
      continue;
    }
    elementParam->m_FittedAmplitudes[k] =
        A * elementParam->m_NominalAmplitudes[k];

    // limit the absorption to 0.0-1.0, so that it's never <0
    //*
    if (elt->getSignFactor(k) == -1 && m_absLinesLimit > 0.0 &&
        elementParam->m_FittedAmplitudes[k] > m_absLinesLimit) {
      elementParam->m_FittedAmplitudes[k] = m_absLinesLimit;
    }

    elementParam->m_FittedAmplitudeErrorSigmas[k] =
        elementParam->m_NominalAmplitudes[k] * 1.0 / sqrt(elt->GetSumGauss());
  }
  return;
}

/**
 * @brief CAbstractFitter::fitAmplitudeAndLambdaOffset
 * Fit the amplitudes and a unique Offset for all lines by using the
 * fitAmplitude() in a loop to tabulate the offset
 * @param spectralAxis
 * @param noContinuumfluxAxis
 * @param continuumfluxAxis
 * @param redshift
 * @param lineIdx
 */
void CAbstractFitter::fitAmplitudeAndLambdaOffset(
    Int32 eltIndex, const CSpectrumSpectralAxis &spectralAxis,
    const CSpectrumFluxAxis &noContinuumfluxAxis,
    const CSpectrumFluxAxis &continuumfluxAxis, Float64 redshift, Int32 lineIdx,
    bool enableOffsetFitting, Float64 step, Float64 min, Float64 max) {
  auto &elt = m_Elements[eltIndex];
  Int32 nLines = elt->GetSize();
  Int32 nSteps = int((max - min) / step + 0.5);

  bool atLeastOneOffsetToFit = false;
  if (enableOffsetFitting) {
    for (Int32 iR = 0; iR < nLines; iR++) {
      // check if the line is to be fitted
      if (elt->GetLines()[iR].GetOffsetFitEnabled()) {
        atLeastOneOffsetToFit = true;
        break;
      }
    }
  }

  if (!atLeastOneOffsetToFit) {
    Log.LogDebug("    multiline: no offsets to fit");
    nSteps = 1;
  } else {
    Log.LogDebug("    multiline: offsets to fit n=%d", nSteps);
  }

  Float64 bestMerit = DBL_MAX;
  Int32 idxBestMerit = -1;
  for (Int32 iO = 0; iO < nSteps; iO++) {
    // set offset value
    if (atLeastOneOffsetToFit)
      elt->SetOffset(min + step * iO);

    // fit for this offset
    fitAmplitude(eltIndex, spectralAxis, noContinuumfluxAxis, continuumfluxAxis,
                 redshift, lineIdx);

    // check fitting
    if (atLeastOneOffsetToFit) {
      Float64 fit = getLeastSquareMeritFast(eltIndex);
      if (fit < bestMerit) {
        bestMerit = fit;
        idxBestMerit = iO;
      }
    }
  }

  if (idxBestMerit >= 0 && atLeastOneOffsetToFit) {
    // set offset value
    if (atLeastOneOffsetToFit)
      elt->SetOffset(min + step * idxBestMerit);

    // fit again for this offset
    fitAmplitude(eltIndex, spectralAxis, noContinuumfluxAxis, continuumfluxAxis,
                 redshift, lineIdx);
  }
}

/**
 * \brief Get the squared difference by fast method proposed by D. Vibert
 **/
Float64 CAbstractFitter::getLeastSquareMeritFast(Int32 idxLine) const {
  Float64 fit = -1; // TODO restore getLeastSquareContinuumMeritFast();
  Int32 istart = 0;
  Int32 iend = m_Elements.size();
  if (idxLine != undefIdx) {
    istart = idxLine;
    iend = idxLine + 1;
  }
  for (Int32 iElts = istart; iElts < iend; iElts++) {
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

TAsymParams CAbstractFitter::fitAsymParameters(Float64 redshift, Int32 idxLyaE,
                                               const Int32 &idxLineLyaE) {

  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc.GetSpectralAxis();
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
        std::ostringstream oss;
        m_Elements.debug(oss);
        Log.LogInfo(oss.str());
        fitAmplitude(idxLyaE, spectralAxis,
                     m_model->getSpcFluxAxisNoContinuum(),
                     m_model->getContinuumFluxAxis(), redshift, idxLineLyaE);
        oss.clear();
        m_Elements.debug(oss);
        Log.LogInfo(oss.str());

        Float64 m = 0; // TODO DV why initializing to m_dTransposeD ?;
        if (1) {

          m_model->refreshModelUnderElements(filterEltsIdxLya, idxLineLyaE);
          m = m_model->getModelErrorUnderElement(idxLyaE,
                                                 m_model->getSpcFluxAxis());
          Log.LogInfo("Fitting Lya profi: merit=%e", m);
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

Int32 CAbstractFitter::fitAsymIGMCorrection(Float64 redshift, Int32 iElts,
                                            const TInt32List &idxLine) {

  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc.GetSpectralAxis();
  if (spectralAxis[0] / (1 + redshift) > RESTLAMBDA_LYA)
    return -1;

  Float64 meritMin = DBL_MAX;
  Int32 bestIgmIdx = -1;

  Int32 igmCount =
      m_Elements[iElts]->getLineProfile(idxLine.front()).getIGMIdxCount();
  for (Int32 igmIdx = 0; igmIdx < igmCount; igmIdx++) {
    m_Elements[iElts]->SetSymIgmParams(TSymIgmParams(igmIdx, redshift));
    fitAmplitude(iElts, spectralAxis, m_model->getSpcFluxAxisNoContinuum(),
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
