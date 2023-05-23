#include "RedshiftLibrary/linemodel/abstractfitter.h"
#include "RedshiftLibrary/linemodel/hybridfitter.h"
#ifdef LBFGSFITTER
#include "RedshiftLibrary/linemodel/gaussianfit/lbfgsfitter.h"
#endif
#include "RedshiftLibrary/linemodel/individualfitter.h"
#include "RedshiftLibrary/linemodel/onesfitter.h"
#include "RedshiftLibrary/linemodel/randomfitter.h"
#include "RedshiftLibrary/linemodel/svdfitter.h"
#include "RedshiftLibrary/linemodel/svdlcfitter.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;

CAbstractFitter::CAbstractFitter(
    const CLMEltListVectorPtr &elementsVector,
    const CCSpectrumVectorPtr &inputSpcs,
    const CTLambdaRangePtrVector &lambdaRanges,
    const CSpcModelVectorPtr &spectrumModels, const TLineVector &restLineList,
    const std::vector<TLineModelElementParam_ptr> &elementParam,
    bool enableAmplitudeOffsets, bool enableLambdaOffsetsFit)
    : m_ElementsVector(elementsVector), m_inputSpcs(inputSpcs),
      m_RestLineList(restLineList), m_lambdaRanges(lambdaRanges),
      m_models(spectrumModels), m_ElementParam(elementParam),
      m_enableAmplitudeOffsets(enableAmplitudeOffsets),
      m_enableLambdaOffsetsFit(enableLambdaOffsetsFit) {

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
    std::string fittingMethod, const CLMEltListVectorPtr &elementsVector,
    const CCSpectrumVectorPtr &inputSpcs,
    const CTLambdaRangePtrVector &lambdaRanges,
    const CSpcModelVectorPtr &spectrumModels, const TLineVector &restLineList,
    std::shared_ptr<CContinuumManager> continuumManager,
    const std::vector<TLineModelElementParam_ptr> &elementParam,
    bool enableAmplitudeOffsets, bool enableLambdaOffsetsFit) {
  if (fittingMethod == "hybrid")
    return std::make_shared<CHybridFitter>(
        elementsVector, inputSpcs, lambdaRanges, spectrumModels, restLineList,
        elementParam, enableAmplitudeOffsets, enableLambdaOffsetsFit);
#ifdef LBFGSFITTER
  else if (fittingMethod == "lbfgs")
    return std::make_shared<CLbfgsFitter>(
        elementsVector, inputSpcs, lambdaRanges, spectrumModels, restLineList,
        elementParam, enableAmplitudeOffsets, enableLambdaOffsetsFit);
#endif
  else if (fittingMethod == "svd")
    return std::make_shared<CSvdFitter>(
        elementsVector, inputSpcs, lambdaRanges, spectrumModels, restLineList,
        elementParam, enableAmplitudeOffsets, enableLambdaOffsetsFit);
  else if (fittingMethod == "svdlc")
    return std::make_shared<CSvdlcFitter>(
        elementsVector, inputSpcs, lambdaRanges, spectrumModels, restLineList,
        elementParam, continuumManager);
  else if (fittingMethod == "svdlcp2")
    return std::make_shared<CSvdlcFitter>(
        elementsVector, inputSpcs, lambdaRanges, spectrumModels, restLineList,
        elementParam, continuumManager, 2);

  else if (fittingMethod == "ones")
    return std::make_shared<COnesFitter>(elementsVector, inputSpcs,
                                         lambdaRanges, spectrumModels,
                                         restLineList, elementParam);
  else if (fittingMethod == "random")
    return std::make_shared<CRandomFitter>(elementsVector, inputSpcs,
                                           lambdaRanges, spectrumModels,
                                           restLineList, elementParam);
  else if (fittingMethod == "individual")
    return std::make_shared<CIndividualFitter>(elementsVector, inputSpcs,
                                               lambdaRanges, spectrumModels,
                                               restLineList, elementParam);
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

void CAbstractFitter::fit(Float64 redshift) {
  initFit(redshift);
  doFit(redshift);
};

void CAbstractFitter::initFit(Float64 redshift) {
  resetSupport(redshift);

  // prepare the Lya width and asym coefficients if the asymfit profile
  // option is met
  fitLyaProfile(redshift);

  resetElementsFittingParam();
}

void CAbstractFitter::resetSupport(Float64 redshift) {

  resetLambdaOffsets();

  // prepare the elements support
  const CSpectrumSpectralAxis &spectralAxis = getSpectrum().GetSpectralAxis();
  for (Int32 iElts = 0; iElts < getElementList().size(); iElts++) {
    getElementList()[iElts]->resetAsymfitParams();
    getElementList()[iElts]->prepareSupport(spectralAxis, redshift,
                                            getLambdaRange());
  }
}

void CAbstractFitter::resetElementsFittingParam() {
  for (Int32 iElts = 0; iElts < getElementList().size(); iElts++) {
    getElementList()[iElts]->reset();
  }
  if (m_enableAmplitudeOffsets)
    getElementList().resetAmplitudeOffset();
}

void CAbstractFitter::resetLambdaOffsets() {
  for (auto &elt : getElementList())
    for (size_t lineIdx = 0; lineIdx < elt->GetSize(); ++lineIdx)
      elt->SetOffset(lineIdx, elt->GetLines()[lineIdx].GetOffset());
}

void CAbstractFitter::fitLyaProfile(Float64 redshift) {
  TInt32List idxEltIGM;
  std::vector<TInt32List> idxLineIGM;
  std::tie(idxEltIGM, idxLineIGM) = getElementList().getIgmLinesIndices();

  if (idxEltIGM.empty())
    return;

  // assuming only one asymfit/fixed profile
  Int32 idxLyaE = idxEltIGM.front();
  Int32 idxLineLyaE = idxLineIGM.front().front();

  const auto &profile = getElementList()[idxLyaE]->getLineProfile(idxLineLyaE);

  if (profile->isAsymFit()) {
    // find the best width and asym coeff. parameters
    TAsymParams bestfitParams =
        fitAsymParameters(redshift, idxLyaE, idxLineLyaE);

    // set the associated Lya members in the element definition
    getElementList()[idxLyaE]->SetAsymfitParams(bestfitParams);
  }

  // deal with symIgm profiles
  for (Int32 i = 0; i < idxEltIGM.size(); ++i) {
    const auto &Elt = getElementList()[idxEltIGM[i]];
    if (!Elt->IsOutsideLambdaRange()) {
      TInt32List &idxLine = idxLineIGM[i];
      auto end =
          std::remove_if(idxLine.begin(), idxLine.end(), [Elt](Int32 idx) {
            return !Elt->GetLines()[idx].GetProfile()->isSymIgmFit();
          });
      idxLine.erase(end, idxLine.end());
      if (!idxLine.empty()) {
        // setSymIgmProfile(idxEltIGM[i], idxLine, redshift);
        auto bestigmidx = fitAsymIGMCorrection(redshift, idxEltIGM[i], idxLine);
        getElementList()[idxEltIGM[i]]->SetSymIgmParams(
            TSymIgmParams(bestigmidx, redshift));
      }
    }
  }
}

void CAbstractFitter::fitAmplitude(Int32 eltIndex, Float64 redshift,
                                   Int32 lineIdx) {
  const CSpectrumSpectralAxis &spectralAxis = getSpectrum().GetSpectralAxis();
  const CSpectrumFluxAxis &noContinuumfluxAxis =
      getModel().getSpcFluxAxisNoContinuum();
  const CSpectrumFluxAxis &continuumfluxAxis =
      getModel().getContinuumFluxAxis();

  getElementList()[eltIndex]->fitAmplitude(
      redshift, spectralAxis, noContinuumfluxAxis, continuumfluxAxis, lineIdx);
}

void CAbstractFitter::setLambdaOffset(const TInt32List &EltsIdx,
                                      Int32 offsetCount) const {

  Float64 offset = m_LambdaOffsetMin + m_LambdaOffsetStep * offsetCount;
  for (Int32 iE : EltsIdx)
    getElementList()[iE]->SetAllOffsetsEnabled(offset);

  return;
}

bool CAbstractFitter::HasLambdaOffsetFitting(TInt32List EltsIdx,
                                             bool enableOffsetFitting) const {
  bool atLeastOneOffsetToFit = false;
  if (enableOffsetFitting) {
    for (Int32 iE : EltsIdx)
      for (const auto &line : getElementList()[iE]->GetLines())
        // check if the line is to be fitted
        if (line.GetOffsetFitEnabled()) {
          atLeastOneOffsetToFit = true;
          break;
        }
  }
  return atLeastOneOffsetToFit;
}

Int32 CAbstractFitter::GetLambdaOffsetSteps(bool atLeastOneOffsetToFit) const {
  Int32 nSteps =
      atLeastOneOffsetToFit
          ? int((m_LambdaOffsetMax - m_LambdaOffsetMin) / m_LambdaOffsetStep +
                0.5)
          : 1;

  return nSteps;
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
void CAbstractFitter::fitAmplitudeAndLambdaOffset(Int32 eltIndex,
                                                  Float64 redshift,
                                                  Int32 lineIdx,
                                                  bool enableOffsetFitting) {

  bool atLeastOneOffsetToFit =
      HasLambdaOffsetFitting({eltIndex}, enableOffsetFitting);
  Int32 nSteps = GetLambdaOffsetSteps(atLeastOneOffsetToFit);

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
      setLambdaOffset({eltIndex}, iO);

    // fit for this offset
    fitAmplitude(eltIndex, redshift, lineIdx);

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
      setLambdaOffset({eltIndex}, idxBestMerit);

    // fit again for this offset
    fitAmplitude(eltIndex, redshift, lineIdx);
  }
}

/**
 * \brief Get the squared difference by fast method proposed by D. Vibert
 **/
Float64 CAbstractFitter::getLeastSquareMeritFast(Int32 idxLine) const {
  Float64 fit = 0.; // TODO restore getLeastSquareContinuumMeritFast();
  Int32 istart = 0;
  Int32 iend = getElementList().size();
  if (idxLine != undefIdx) {
    istart = idxLine;
    iend = idxLine + 1;
  }
  for (Int32 iElts = istart; iElts < iend; iElts++) {
    Float64 dtm = getElementList()[iElts]->GetSumCross();
    Float64 mtm = getElementList()[iElts]->GetSumGauss();
    Float64 a = getElementList().GetElementAmplitude(
        iElts); //[iElts]->GetFitAmplitude();
    Float64 term1 = a * a * mtm;
    Float64 term2 = -2. * a * dtm;
    fit += term1 + term2;
  }

  Log.LogDebug("CLineModelFitting::getLeastSquareMerit fit fast = %f", fit);
  return fit;
}

TAsymParams CAbstractFitter::fitAsymParameters(Float64 redshift, Int32 idxLyaE,
                                               const Int32 &idxLineLyaE) {

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
        getElementList()[idxLyaE]->SetAsymfitParams(
            {asymWidthCoeff, asymAlphaCoeff, delta});

        // idxLineLyaE = -1;
        fitAmplitude(idxLyaE, redshift, idxLineLyaE);

        Float64 m = 0; // TODO DV why initializing to m_dTransposeD ?;
        if (1) {

          getModel().refreshModelUnderElements(filterEltsIdxLya, idxLineLyaE);
          m = getModel().getModelErrorUnderElement(idxLyaE,
                                                   getModel().getSpcFluxAxis());
        } else {
          m = getLeastSquareMeritFast(idxLyaE);
        }
        if (m < meritMin) {
          meritMin = m;
          bestparams = getElementList()[idxLyaE]->GetAsymfitParams(0);
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

  const CSpectrumSpectralAxis &spectralAxis = getSpectrum().GetSpectralAxis();
  if (spectralAxis[0] / (1 + redshift) > RESTLAMBDA_LYA)
    return -1;

  Float64 meritMin = DBL_MAX;
  Int32 bestIgmIdx = -1;

  Int32 igmCount = getElementList()[iElts]
                       ->getLineProfile(idxLine.front())
                       ->getIGMIdxCount();
  for (Int32 igmIdx = 0; igmIdx < igmCount; igmIdx++) {
    getElementList()[iElts]->SetSymIgmParams(TSymIgmParams(igmIdx, redshift));
    fitAmplitude(iElts, redshift);

    getModel().refreshModelUnderElements(TInt32List(1, iElts));
    Float64 m = getModel().getModelErrorUnderElement(
        iElts, getModel().getSpcFluxAxis());

    if (m < meritMin) {
      meritMin = m;
      bestIgmIdx = igmIdx;
    }
  }
  return bestIgmIdx;
}
