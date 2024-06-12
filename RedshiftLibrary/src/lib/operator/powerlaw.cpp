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

#include "RedshiftLibrary/operator/powerlaw.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/vectorOperations.h"
#include "RedshiftLibrary/operator/powerlawresult.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/spectrum/template/template.h"

#include <Eigen/Dense>
#include <cmath>
#include <execution>

using namespace NSEpic;
using namespace std;

// TODO Voir où les shared ptr sont à ajouter ici / syntaxe à adapter
// TODO adapt to multiobs
// TODO in tests, compare matrix inversion by formula / by lib (see
// CSvdlcFitter::fitAmplitudesLinesAndContinuumLinSolve)

COperatorPowerLaw::COperatorPowerLaw(const TFloat64List &redshifts,
                                     Float64 lambdaCut)
    : COperatorPowerLawBase(redshifts), m_lambdaCut(lambdaCut) {

  // Gets pixels of interest
  const CSpectrumSpectralAxis &spectrumLambda =
      Context.getSpectra()[0]->GetSpectralAxis();
  m_lambdaRanges[0]->getClosedIntervalIndices(spectrumLambda.GetSamplesVector(),
                                              m_firstPixelIdxInRange,
                                              m_lastPixelIdxInRange);
  m_nPixels = m_lastPixelIdxInRange - m_firstPixelIdxInRange + 1;
}

TPowerLawResult
COperatorPowerLaw::BasicFit(Float64 redshift, bool opt_extinction,
                            bool opt_dustFitting, Float64 nullFluxThreshold,
                            Int32 nLogSamplesMin, std::string method) {

  // Igm ism curves
  if (opt_extinction) {
    m_igmCorrectionMeiksin = Context.getFluxCorrectionMeiksin();
    m_nIgmCurves = m_igmCorrectionMeiksin->getIdxCount();
  } else {
    m_nIgmCurves = 1;
  }
  if (opt_dustFitting) {
    m_ismCorrectionCalzetti = Context.getFluxCorrectionCalzetti();
    m_nIsmCurves = m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs();
  } else {
    m_nIsmCurves = 1;
  }

  T3DCurve emittedCurve = computeEmittedCurve(redshift, nullFluxThreshold,
                                              opt_extinction, opt_dustFitting);

  // Calculates log / log
  T3DCurve lnCurves = computeLnCurves(emittedCurve);

  // TODO continue here noise inclusion

  // Step 3. Compute power law coefs and chi2
  T2DPowerLawCoefsPair coefs = powerLawCoefs3D(
      lnCurves, TInt32Range(0, m_nPixels), nLogSamplesMin, method);
  TChi2Result chi2Result =
      findMinChi2OnIgmIsm({emittedCurve.lambda, emittedCurve.flux,
                           lnCurves.pixelsToUse, emittedCurve.fluxError},
                          coefs);

  // Step 4. Creates result

  // TODO add error here
  TPowerLawResult result;
  result.chiSquare = chi2Result.chi2;
  result.igmIdx = chi2Result.igmIdx;
  result.ismIdx = chi2Result.ismIdx;
  result.coefs = coefs[chi2Result.igmIdx][chi2Result.ismIdx];
  result.lambdaCut = m_lambdaCut;
  result.lambdaRest = emittedCurve.lambda;
  result.emittedFlux = emittedCurve.flux[chi2Result.igmIdx][chi2Result.ismIdx];
  result.fluxError =
      emittedCurve.fluxError[chi2Result.igmIdx][chi2Result.ismIdx];
  result.pixelsToUse =
      lnCurves.pixelsToUse[chi2Result.igmIdx][chi2Result.ismIdx];
  return result;
};
// TODO see the different border cases (ex coefs set to 0)

T3DCurve COperatorPowerLaw::computeLnCurves(T3DCurve const &emittedCurve) {

  T3DCurve lnCurve = emittedCurve;
  lnCurve.lambda = lnLambda(emittedCurve.lambda);
  for (Int32 igmIdx = 0; igmIdx < m_nIgmCurves; igmIdx++) {
    for (Int32 ismIdx = 0; ismIdx < m_nIsmCurves; ismIdx++) {
      for (Int32 pixelIdx = m_firstPixelIdxInRange;
           pixelIdx <= m_lastPixelIdxInRange; pixelIdx++) {
        if (!emittedCurve.pixelsToUse[igmIdx][ismIdx][pixelIdx])
          continue;
        lnCurve.flux[igmIdx][ismIdx][pixelIdx] =
            std::log(emittedCurve.flux[igmIdx][ismIdx][pixelIdx]);
        // NB here we assume that too low values of flux have been masked (done
        // before with SNR threshold) We use sigma_log = sigma / flux *
        // lambdarest (see doc)
        lnCurve.fluxError[igmIdx][ismIdx][pixelIdx] =
            emittedCurve.fluxError[igmIdx][ismIdx][pixelIdx] /
            emittedCurve.flux[igmIdx][ismIdx][pixelIdx] *
            emittedCurve.lambda[pixelIdx];
      }
    }
  }

  return lnCurve;
}

// TODO rename curve3D
TChi2Result
COperatorPowerLaw::findMinChi2OnIgmIsm(T3DCurve const &curve3D,
                                       T2DPowerLawCoefsPair const &coefs) {
  T2DFloatList chi2 = computeChi2(curve3D, coefs);
  TInt32Pair minChi2Idxs = find2DVectorMinIndexes(chi2);
  return {minChi2Idxs.first, minChi2Idxs.second,
          chi2[minChi2Idxs.first][minChi2Idxs.second]};
}

Float64 COperatorPowerLaw::theoreticalFluxAtLambda(TPowerLawCoefsPair fullCoefs,
                                                   Float64 lambda) {
  Float64 theoreticalFlux = NAN;
  if (lambda < m_lambdaCut)
    theoreticalFlux = computePowerLaw(fullCoefs.first, lambda);
  else
    theoreticalFlux = computePowerLaw(fullCoefs.second, lambda);

  return theoreticalFlux;
}

Float64 COperatorPowerLaw::computePowerLaw(TPowerLawCoefs coefs,
                                           Float64 lambda) {
  return coefs.a * std::pow(lambda, coefs.b);
}

T2DFloatList COperatorPowerLaw::computeChi2(T3DCurve const &curve3D,
                                            T2DPowerLawCoefsPair const &coefs) {
  T2DFloatList chi2(m_nIgmCurves, std::vector<Float64>(m_nIsmCurves, 0));
  Float64 diff;
  for (Int16 igmIdx = 0; igmIdx < m_nIgmCurves; igmIdx++) {
    for (Int16 ismIdx = 0; ismIdx < m_nIsmCurves; ismIdx++) {
      for (Int32 pixelIdx = 0; pixelIdx < m_nPixels; pixelIdx++) {
        if (curve3D.pixelsToUse[igmIdx][ismIdx][pixelIdx]) {
          diff = curve3D.flux[igmIdx][ismIdx][pixelIdx] -
                 theoreticalFluxAtLambda(
                     TPowerLawCoefsPair(coefs[igmIdx][ismIdx].first,
                                        coefs[igmIdx][ismIdx].second),
                     curve3D.lambda[pixelIdx]);
          diff = diff / curve3D.fluxError[igmIdx][ismIdx][pixelIdx];
          chi2[igmIdx][ismIdx] += diff * diff;
        }
      }
    }
  }
  return chi2;
}

TAxisSampleList COperatorPowerLaw::lnLambda(TAxisSampleList const &lambda) {
  TAxisSampleList lnLambda(m_nPixels);
  std::transform(lambda.begin(), lambda.end(), lnLambda.begin(),
                 [](float value) { return std::log(value); });
  return lnLambda;
}

std::shared_ptr<COperatorResult> COperatorPowerLaw::Compute() {

  // Creates power law result
  std::shared_ptr<CPowerLawResult> result =
      std::make_shared<CPowerLawResult>(m_redshifts.size());

  // TODO see that later
  // // Loops on redshifts
  // for (Int32 i = 0; i < m_redshifts.size(); i++) {

  //   // Question : is there a diff in value with m_redshifts[i] ?
  //   Float64 redshift = result->Redshifts[i];
  // }

  return result;
};

T2DPowerLawCoefsPair COperatorPowerLaw::powerLawCoefs3D(
    T3DCurve const &lnCurves, TInt32Range const pixelsRange,
    Int32 nLogSamplesMin, std::string method) const {
  // All 3D vectors must be of the same size > 0

  Int32 nIgm = lnCurves.flux.size();
  Int32 nIsm = lnCurves.flux[0].size();

  T2DPowerLawCoefsPair powerLawsCoefs(nIgm,
                                      std::vector<TPowerLawCoefsPair>(nIsm));
  for (Int32 igmIdx = 0; igmIdx < nIgm; igmIdx++) {
    for (Int32 ismIdx = 0; ismIdx < nIsm; ismIdx++) {
      if (method == "full")
        powerLawsCoefs[igmIdx][ismIdx] = computeFullPowerLawCoefs(
            {lnCurves.lambda, lnCurves.flux[igmIdx][ismIdx],
             lnCurves.pixelsToUse[igmIdx][ismIdx],
             lnCurves.fluxError[igmIdx][ismIdx]},
            pixelsRange, nLogSamplesMin);
      else if (method == "simple")
        powerLawsCoefs[igmIdx][ismIdx] = {
            computeSimplePowerLawCoefs({lnCurves.lambda,
                                        lnCurves.flux[igmIdx][ismIdx],
                                        lnCurves.pixelsToUse[igmIdx][ismIdx],
                                        lnCurves.fluxError[igmIdx][ismIdx]},
                                       pixelsRange, nLogSamplesMin, false),
            {0, 0}};
      else if (method == "simpleWeighted")
        powerLawsCoefs[igmIdx][ismIdx] = {
            computeSimplePowerLawCoefs({lnCurves.lambda,
                                        lnCurves.flux[igmIdx][ismIdx],
                                        lnCurves.pixelsToUse[igmIdx][ismIdx],
                                        lnCurves.fluxError[igmIdx][ismIdx]},
                                       pixelsRange, nLogSamplesMin, true),
            {0, 0}};
      else
        THROWG(ErrorCode::INTERNAL_ERROR,
               Formatter() << "Unexepected methd " << method);
    }
  }
  return powerLawsCoefs;
}

void COperatorPowerLaw::checkInputPowerLawCoefs(
    TCurve const &lnCurve, TInt32Range const pixelsRange) const {
  Int32 n = lnCurve.lambda.size();

  // TODO refacto these errors
  if (n != lnCurve.flux.size() || n != lnCurve.fluxError.size())
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Cannot calculate power law coefs: both ln(lambda) "
                          "and ln(flux) and ln (fluxError)"
                          "vectors should have same size. Sizes : ln(lambda) "
                       << n << ", ln(flux) : " << lnCurve.flux.size()
                       << ", ln(fluxError) : " << lnCurve.fluxError.size());
  if (n != lnCurve.pixelsToUse.size())
    THROWG(ErrorCode::INTERNAL_ERROR,
           "Cannot calculate power law coefs: both ln(lambda) and pixelsToUse "
           "vectors should have same size");
  if (n == 0)
    THROWG(ErrorCode::INTERNAL_ERROR,
           "Cannot calculate power law coefs: empty input vector");
}

TPowerLawCoefs COperatorPowerLaw::computeSimplePowerLawCoefs(
    TCurve const &lnCurve, TInt32Range const pixelsRange, Int32 nLogSamplesMin,
    bool applyWeight) const {
  checkInputPowerLawCoefs(lnCurve, pixelsRange);

  Float64 SX = 0;
  Float64 SY = 0;
  Float64 SXY = 0;
  Float64 SX2 = 0;
  Float64 n = 0;
  Int32 N = 0;

  for (Int32 pixelIdx = pixelsRange.GetBegin(); pixelIdx < pixelsRange.GetEnd();
       pixelIdx++) {
    Float64 w = applyWeight ? lnCurve.fluxError[pixelIdx] : 1;
    Float64 X = lnCurve.lambda[pixelIdx] / w;
    Float64 Y = lnCurve.flux[pixelIdx] / w;
    SX += X / w;
    SY += Y / w;
    SXY += X * Y;
    SX2 += X * X;
    n += 1 / (w * w);
    N += 1;
  }

  if (N < nLogSamplesMin) {
    Flag.warning(WarningCode::FORCED_POWERLAW_TO_ZERO,
                 Formatter() << "COperatorPowerLaw::" << __func__ << ": only "
                             << N << " < " << nLogSamplesMin
                             << " samples with significant flux values. Power "
                                "law coefs are forced to zero.");
    return {0., 0.};
  }

  Float64 denom = n * SX2 - SX * SX;
  Float64 b = (n * SXY - SX * SY) / denom;
  Float64 a = std::exp((SY - b * SX) / n);
  return {a, b};
}

// TODO Split in several methods
// TODO attention où est le noise là dedans ?
TPowerLawCoefsPair
COperatorPowerLaw::computeFullPowerLawCoefs(TCurve const &lnCurve,
                                            TInt32Range const pixelsRange,
                                            Int32 nLogSamplesMin) const {

  checkInputPowerLawCoefs(lnCurve, pixelsRange);

  Float64 xc = std::log(m_lambdaCut);

  // 1 <-> power law first part / 2 <-> second part

  // n pixels
  Int32 N1 = 0;
  Int32 N2 = 0;
  // sum pixels weights
  Float64 n1 = 0;
  Float64 n2 = 0;
  // sum ln xi * wi
  Float64 sx1 = 0;
  Float64 sx2 = 0;
  // sum (ln xi)^2 * wi
  Float64 sxx1 = 0;
  Float64 sxx2 = 0;
  // sum ln yi * wi
  Float64 sy1 = 0;
  Float64 sy2 = 0;
  // sum ln xi * ln yi * wi
  Float64 sxy1 = 0;
  Float64 sxy2 = 0;

  // sum wi*(lnxc - lnxi)**2 (on second part only)
  Float64 sx2mc2 = 0;

  // Loop over all pixels to make necessary pre-calculations
  for (Int32 pixelIdx = pixelsRange.GetBegin(); pixelIdx < pixelsRange.GetEnd();
       pixelIdx++) {

    if (!lnCurve.pixelsToUse[pixelIdx])
      continue;

    Float64 xi = lnCurve.lambda[pixelIdx];
    Float64 yi = lnCurve.flux[pixelIdx];
    Float64 wi =
        1 / (lnCurve.fluxError[pixelIdx] * lnCurve.fluxError[pixelIdx]);

    if (xi < xc) {
      n1 += wi;
      sx1 += xi * wi;
      sxx1 += xi * xi * wi;
      sy1 += yi * wi;
      sxy1 += xi * yi * wi;
      N1 += 1;
    } else {
      n2 += wi;
      sx2 += xi * wi;
      sxx2 += xi * xi * wi;
      sy2 += yi * wi;
      sxy2 += xi * yi * wi;
      sx2mc2 += wi * (xc - xi) * (xc - xi);
      N2 += 1;
    }
  }

  if (N1 < nLogSamplesMin || N2 < nLogSamplesMin) {
    Flag.warning(WarningCode::FORCED_POWERLAW_TO_ZERO,
                 Formatter() << "COperatorPowerLaw::" << __func__ << ": only "
                             << N1 << " and " << N2 << " < " << nLogSamplesMin
                             << " samples with significant flux values. Power "
                                "law coefs are forced to zero.");
    return {{0, 0}, {0, 0}};
  }

  // Creates the Mc^TN^-1Mc matrix
  Float64 m11 = n1 + n2;
  Float64 m12 = n2 * xc + sx1;
  Float64 m13 = -n2 * xc + sx2;
  Float64 m21 = n2 * xc + sx1;
  Float64 m22 = n2 * xc * xc + sxx1;
  Float64 m23 = -n2 * xc * xc + xc * sx2;
  Float64 m31 = -n2 * xc + sx2;
  Float64 m32 = -n2 * xc * xc + xc * sx2;
  Float64 m33 = sx2mc2;

  // Creates the Mc^TY vector
  Float64 v1 = sy1 + sy2;
  Float64 v2 = sxy1 + xc * sy2;
  Float64 v3 = sxy2 - xc * sy2;

  Eigen::Matrix3d m;
  m << m11, m12, m13, m21, m22, m23, m31, m32, m33;

  Eigen::Vector3d v(v1, v2, v3);

  if (std::abs(m.determinant()) < DBL_MIN)
    THROWG(ErrorCode::INTERNAL_ERROR,
           "Cannot calculate power law coefs: division by zero");

  Eigen::Vector3d theta = m.inverse() * v;

  Float64 a1 = std::exp(theta(0));
  Float64 b1 = theta(1);
  Float64 b2 = theta(2);
  Float64 a2 = a1 * std::pow(xc, b1 - b2);

  return {{a1, b1}, {a2, b2}};
}

// TODO here initialize flux curve
T3DCurve
COperatorPowerLaw::initializeFluxCurve(Float64 redshift,
                                       Float64 nullFluxThreshold) const {
  // In order to take into account ism, igm we initialize a template with a flux
  // at 1 we then divide the initial flux by the template flux value once ism
  // igm has been applied

  const CSpectrumSpectralAxis &spectrumLambdaAxis =
      Context.getSpectra()[0]->GetSpectralAxis();
  // Flux
  const CSpectrumFluxAxis &spectrumFluxAxis =
      Context.getSpectra()[0]->GetFluxAxis();
  // Error
  const CSpectrumFluxAxis &spectrumFluxErrorAxis =
      Context.getSpectra()[0]->GetFluxAxis().GetError();

  // Step 2. Transform spectrum data
  // From lambda obs to lambda rest
  CSpectrumSpectralAxis spectrumLambdaRest =
      spectrumLambdaAxis.blueShift(redshift);

  // Initializes usable pixels
  T3DBoolList pixelsToUse(
      m_nIgmCurves, std::vector<std::vector<bool>>(
                        m_nIsmCurves, std::vector<bool>(m_nPixels, true)));

  TBoolList snrCompliantPixels = computeSNRCompliantPixels(
      spectrumFluxAxis.GetSamplesVector(),
      spectrumFluxErrorAxis.GetSamplesVector(), nullFluxThreshold);

  for (Int32 pixelIdx = 0; pixelIdx < m_nPixels; pixelIdx++) {
    if (!snrCompliantPixels[pixelIdx]) {
      for (Int32 igmIdx = 0; igmIdx < m_nIgmCurves; igmIdx++) {
        for (Int32 ismIdx = 0; ismIdx < m_nIsmCurves; ismIdx++) {
          pixelsToUse[igmIdx][ismIdx][pixelIdx] = false;
        }
      }
    }
  }

  T3DFloatList flux(m_nIgmCurves,
                    std::vector<std::vector<Float64>>(
                        m_nIsmCurves, spectrumFluxAxis.GetSamplesVector()));
  T3DFloatList fluxError(
      m_nIgmCurves,
      std::vector<std::vector<Float64>>(
          m_nIsmCurves, spectrumFluxErrorAxis.GetSamplesVector()));

  return {spectrumLambdaRest.GetSamplesVector(), flux, pixelsToUse, fluxError};
}

T3DFloatList COperatorPowerLaw::computeIsmIgmCorrections(
    Float64 redshift, CSpectrumSpectralAxis const &spectrumLambdaRest,
    bool opt_extinction, bool opt_dustFitting) const {
  // In order to access ism igm coefs, we initialize a template with a flux
  // at 1, and apply ism/igm on it

  CTemplate templateForCoefs("", "", spectrumLambdaRest,
                             std::vector<Float64>(m_nPixels, 1));
  templateForCoefs.InitIsmIgmConfig(redshift, m_ismCorrectionCalzetti,
                                    m_igmCorrectionMeiksin);

  T3DFloatList correctionCoefs(
      m_nIgmCurves, std::vector<std::vector<Float64>>(
                        m_nIsmCurves, std::vector<Float64>(m_nPixels, NAN)));

  Float64 igmCorrection = 1;
  Float64 ismCorrection = 1;
  for (Int32 igmIdx = 0; igmIdx < m_nIgmCurves; igmIdx++) {
    if (opt_extinction) { // igm
      templateForCoefs.ApplyMeiksinCoeff(igmIdx);
    }
    for (Int32 ismIdx = 0; ismIdx < m_nIsmCurves; ismIdx++) {
      if (opt_dustFitting) { // ism
        templateForCoefs.ApplyDustCoeff(ismIdx);
      }
      correctionCoefs[igmIdx][ismIdx] =
          templateForCoefs.GetFluxAxis().GetSamplesVector();
    }
  }
  return correctionCoefs;
}

T3DCurve COperatorPowerLaw::computeEmittedCurve(Float64 redshift,
                                                Float64 nullFluxThreshold,
                                                bool opt_extinction,
                                                bool opt_dustFitting) const {
  // In order to take into account ism, igm we initialize a template with a flux
  // at 1 we then divide the initial flux by the template flux value once ism
  // igm has been applied

  T3DCurve fluxCurve = initializeFluxCurve(redshift, nullFluxThreshold);

  if (opt_extinction || opt_dustFitting) {
    T3DFloatList correctionCoefs = computeIsmIgmCorrections(
        redshift, fluxCurve.lambda, opt_extinction, opt_dustFitting);
    for (Int32 igmIdx = 0; igmIdx < m_nIgmCurves; igmIdx++) {
      for (Int32 ismIdx = 0; ismIdx < m_nIsmCurves; ismIdx++) {
        for (Int32 pixelIdx = 0; pixelIdx < m_nPixels; pixelIdx++) {
          if (!fluxCurve.pixelsToUse[igmIdx][ismIdx][pixelIdx])
            continue;
          // If ism & igm completely absorb the flux, do not use this pixel
          Float64 correctionCoef = correctionCoefs[igmIdx][ismIdx][pixelIdx];
          if (correctionCoef < DBL_MIN) {
            fluxCurve.pixelsToUse[igmIdx][ismIdx][pixelIdx] = false;
            continue;
          }
          fluxCurve.flux[igmIdx][ismIdx][pixelIdx] =
              fluxCurve.flux[igmIdx][ismIdx][pixelIdx] / correctionCoef;
          fluxCurve.fluxError[igmIdx][ismIdx][pixelIdx] =
              fluxCurve.fluxError[igmIdx][ismIdx][pixelIdx] / correctionCoef;
        }
      }
    }
  }
  return fluxCurve;
}

TBoolList COperatorPowerLaw::computeSNRCompliantPixels(
    TFloat64List const &spectrumFlux, TFloat64List const &spectrumFluxError,
    Float64 nullFluxThreshold) const {
  // Masks pixels with insufficient SNR. In place modification of
  // curve.pixelsToUse
  TBoolList pixelsToUse = std::vector<bool>(m_nPixels, true);
  for (Int32 pixelIdx = m_firstPixelIdxInRange;
       pixelIdx <= m_lastPixelIdxInRange; pixelIdx++) {
    // Checks that signal / noise ratio is sufficient (> nullFluxThreshold)
    if (spectrumFlux[pixelIdx] <=
        nullFluxThreshold * spectrumFluxError[pixelIdx]) {
      pixelsToUse[pixelIdx] = false;
    }
  }
  return pixelsToUse;
}
