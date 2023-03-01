#include "RedshiftLibrary/linemodel/abstractfitter.h"
#include "RedshiftLibrary/linemodel/hybridfitter.h"
#include "RedshiftLibrary/linemodel/individualfitter.h"
#include "RedshiftLibrary/linemodel/onesfitter.h"
#include "RedshiftLibrary/linemodel/randomfitter.h"
#include "RedshiftLibrary/linemodel/svdfitter.h"
#include "RedshiftLibrary/linemodel/svdlcfitter.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;

CAbstractFitter::CAbstractFitter(
    CLineModelElementList &elements,
    std::shared_ptr<const CSpectrum> inputSpectrum,
    std::shared_ptr<const TLambdaRange> lambdaRange,
    std::shared_ptr<CSpectrumModel> spectrumModel,
    const CLineCatalog::TLineVector &restLineList,
    const std::vector<std::shared_ptr<TFittedData>> &fittedData)
    : m_Elements(elements), m_inputSpc(*(inputSpectrum)),
      m_RestLineList(restLineList), m_lambdaRange(*(lambdaRange)),
      m_model(spectrumModel), m_fittedData(fittedData) {}

std::shared_ptr<CAbstractFitter> CAbstractFitter::makeFitter(
    std::string fittingMethod, CLineModelElementList &elements,
    std::shared_ptr<const CSpectrum> inputSpectrum,
    std::shared_ptr<const TLambdaRange> lambdaRange,
    std::shared_ptr<CSpectrumModel> spectrumModel,
    const CLineCatalog::TLineVector &restLineList,
    std::shared_ptr<CContinuumManager> continuumManager,
    const std::vector<std::shared_ptr<TFittedData>> &fittedData) {
  if (fittingMethod == "hybrid")
    return std::make_shared<CHybridFitter>(
        CHybridFitter(elements, inputSpectrum, lambdaRange, spectrumModel,
                      restLineList, fittedData));
  else if (fittingMethod == "svd")
    return std::make_shared<CSvdFitter>(CSvdFitter(elements, inputSpectrum,
                                                   lambdaRange, spectrumModel,
                                                   restLineList, fittedData));
  else if (fittingMethod == "svdlc")
    return std::make_shared<CSvdlcFitter>(
        CSvdlcFitter(elements, inputSpectrum, lambdaRange, spectrumModel,
                     restLineList, fittedData, continuumManager));
  else if (fittingMethod == "svdlcp2")
    return std::make_shared<CSvdlcFitter>(
        CSvdlcFitter(elements, inputSpectrum, lambdaRange, spectrumModel,
                     restLineList, fittedData, continuumManager, 2));

  else if (fittingMethod == "ones")
    return std::make_shared<COnesFitter>(COnesFitter(elements, inputSpectrum,
                                                     lambdaRange, spectrumModel,
                                                     restLineList, fittedData));
  else if (fittingMethod == "random")
    return std::make_shared<CRandomFitter>(
        CRandomFitter(elements, inputSpectrum, lambdaRange, spectrumModel,
                      restLineList, fittedData));
  else if (fittingMethod == "individual")
    return std::make_shared<CIndividualFitter>(
        CIndividualFitter(elements, inputSpectrum, lambdaRange, spectrumModel,
                          restLineList, fittedData));
  else
    THROWG(INTERNAL_ERROR, Formatter()
                               << "Unknown fitting method " << fittingMethod);
}

void CAbstractFitter::fitAmplitude(Int32 eltIndex,
                                   const CSpectrumSpectralAxis &spectralAxis,
                                   const CSpectrumFluxAxis &noContinuumfluxAxis,
                                   const CSpectrumFluxAxis &continuumfluxAxis,
                                   Float64 redshift, Int32 lineIdx) {
  Int32 nLines = m_Elements[eltIndex]->m_Lines.size();

  m_fittedData[eltIndex]->m_FittedAmplitudes.assign(nLines, NAN);
  m_fittedData[eltIndex]->m_FittedAmplitudeErrorSigmas.assign(nLines, NAN);

  if (m_Elements[eltIndex]->IsOutsideLambdaRange()) {
    m_fittedData[eltIndex]->m_fitAmplitude = NAN;
    m_sumCross = NAN;
    m_sumGauss = NAN;
    m_dtmFree = NAN;
    return;
  }

  m_sumCross = 0.0;
  m_sumGauss = 0.0;
  m_dtmFree = 0.0;

  const CSpectrumNoiseAxis &error = noContinuumfluxAxis.GetError();

  Float64 y = 0.0;
  Float64 x = 0.0;
  Float64 yg = 0.0;
  Float64 c = 1.0;

  Float64 err2 = 0.0;
  Int32 num = 0;
  // Log.LogDebug("    multiline: nLines=%d", nLines);

  for (Int32 k = 0; k < nLines; k++) { // loop for the intervals
    if (m_Elements[eltIndex]->IsOutsideLambdaRange(k)) {
      continue;
    }
    if (lineIdx != undefIdx &&
        !m_Elements[eltIndex]->isLineActiveOnSupport(k, lineIdx)) {
      continue;
    }

    // A estimation
    //#pragma omp parallel for
    for (Int32 i = m_Elements[eltIndex]->getStartNoOverlap(k);
         i <= m_Elements[eltIndex]->getEndNoOverlap(k); i++) {
      c = continuumfluxAxis[i];
      y = noContinuumfluxAxis[i];
      x = spectralAxis[i];

      yg = 0.0;

      for (Int32 k2 = 0; k2 < nLines; k2++) { // loop for the signal synthesis
        if (m_Elements[eltIndex]->IsOutsideLambdaRange(k2) ||
            !m_Elements[eltIndex]->isLineActiveOnSupport(k2, k)) {
          continue;
        }
        Int32 sf = m_Elements[eltIndex]->getSignFactor(k2);
        if (sf == -1) {
          yg += sf * c * m_fittedData[eltIndex]->m_NominalAmplitudes[k2] *
                m_Elements[eltIndex]->GetLineProfileAtRedshift(k2, redshift, x);
        } else {
          yg += sf * m_fittedData[eltIndex]->m_NominalAmplitudes[k2] *
                m_Elements[eltIndex]->GetLineProfileAtRedshift(k2, redshift, x);
        }
      }
      num++;
      err2 = 1.0 / (error[i] * error[i]);
      m_dtmFree += yg * y * err2;
      m_sumGauss += yg * yg * err2;
    }
  }

  if (num == 0 || m_sumGauss == 0) {
    Log.LogDebug("CLineModelElement::fitAmplitude: Could not fit amplitude:    "
                 " num=%d, mtm=%f",
                 num, m_sumGauss);
    /*for (Int32 k2 = 0; k2 < nLines; k2++) {
      Log.LogDebug("    multiline failed:     subE=%d, nominal_amp=%f", k2,
                   m_fittedData[eltIndex]->m_NominalAmplitudes[k2]);
    }*/
    m_fittedData[eltIndex]->m_fitAmplitude = NAN;
    m_sumCross = NAN;
    m_sumGauss = NAN;
    m_dtmFree = NAN;
    return;
  }

  m_sumCross = std::max(0.0, m_dtmFree);
  Float64 A = m_sumCross / m_sumGauss;
  m_fittedData[eltIndex]->m_fitAmplitude =
      A; // todo: warning m_fitAmplitude should be updated when
         // modifying sub-elements amplitudes: ex. rules.
  // Float64 A = std::max(0.0, m_sumCross / m_sumGauss);

  for (Int32 k = 0; k < nLines; k++) {
    if (m_Elements[eltIndex]->IsOutsideLambdaRange(k)) {
      continue;
    }
    m_fittedData[eltIndex]->m_FittedAmplitudes[k] =
        A * m_fittedData[eltIndex]->m_NominalAmplitudes[k];

    // limit the absorption to 0.0-1.0, so that it's never <0
    //*
    if (m_Elements[eltIndex]->getSignFactor(k) == -1 && m_absLinesLimit > 0.0 &&
        m_fittedData[eltIndex]->m_FittedAmplitudes[k] > m_absLinesLimit) {
      m_fittedData[eltIndex]->m_FittedAmplitudes[k] = m_absLinesLimit;
    }
    //*/

    //        if(A==0)
    //        {
    //            m_FittedAmplitudeErrorSigmas[k] = 0.0; //why would this be
    //            useful ?
    //        }
    //        else
    {
      m_fittedData[eltIndex]->m_FittedAmplitudeErrorSigmas[k] =
          m_fittedData[eltIndex]->m_NominalAmplitudes[k] * 1.0 /
          sqrt(m_sumGauss);
    }
  }
  return;
}

/**
 * @brief CLineModelElement::fitAmplitudeAndLambdaOffset
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
  Int32 nLines = m_Elements[eltIndex]->m_Lines.size();
  Int32 nSteps = int((max - min) / step + 0.5);

  bool atLeastOneOffsetToFit = false;
  if (enableOffsetFitting) {
    for (Int32 iR = 0; iR < nLines; iR++) {
      // check if the line is to be fitted
      if (m_Elements[eltIndex]->m_Lines[iR].GetOffsetFitEnabled()) {
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
    if (atLeastOneOffsetToFit) {
      Float64 offset = min + step * iO;
      for (Int32 iR = 0; iR < nLines; iR++) {
        if (m_Elements[eltIndex]->m_Lines[iR].GetOffsetFitEnabled()) {
          m_Elements[eltIndex]->m_Lines[iR].SetOffset(offset);
        }
      }
    }

    // fit for this offset
    fitAmplitude(eltIndex, spectralAxis, noContinuumfluxAxis, continuumfluxAxis,
                 redshift, lineIdx);

    // check fitting
    if (atLeastOneOffsetToFit) {
      Float64 dtm = m_Elements[eltIndex]->GetSumCross();
      Float64 mtm = m_Elements[eltIndex]->GetSumGauss();
      Float64 a = m_Elements[eltIndex]->GetFitAmplitude();
      Float64 term1 = a * a * mtm;
      Float64 term2 = -2. * a * dtm;
      Float64 fit = term1 + term2;
      if (fit < bestMerit) {
        bestMerit = fit;
        idxBestMerit = iO;
      }
    }
  }

  if (idxBestMerit >= 0 && atLeastOneOffsetToFit) {
    // set offset value
    if (atLeastOneOffsetToFit) {
      Float64 offset = min + step * idxBestMerit;
      Log.LogDebug("    multiline: offset best found=%f", offset);
      for (Int32 iR = 0; iR < nLines; iR++) {
        if (m_Elements[eltIndex]->m_Lines[iR].GetOffsetFitEnabled()) {
          m_Elements[eltIndex]->m_Lines[iR].SetOffset(offset);
        }
      }
    }
    // fit again for this offset
    fitAmplitude(eltIndex, spectralAxis, noContinuumfluxAxis, continuumfluxAxis,
                 redshift, lineIdx);
  }
}
