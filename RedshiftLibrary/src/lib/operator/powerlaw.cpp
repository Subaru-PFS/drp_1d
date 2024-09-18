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
#include "RedshiftLibrary/common/curve3d.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/vectorOperations.h"
#include "RedshiftLibrary/operator/continuumfitting.h"
#include "RedshiftLibrary/operator/powerlawresult.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/spectrum/template/template.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <execution>
#include <utility>

using namespace NSEpic;
using namespace std;

COperatorPowerLaw::COperatorPowerLaw(const TFloat64List &redshifts,
                                     Float64 lambdaCut)
    : COperatorContinuumFitting(redshifts), m_lambdaCut(lambdaCut) {

  // Gets pixels of interest and sets elements size
  m_spectra = Context.getSpectra();
  m_nSpectra = m_spectra.size();
  m_nPixels.resize(m_nSpectra);
  m_firstPixelIdxInRange.resize(m_nSpectra);
  m_lastPixelIdxInRange.resize(m_nSpectra);

  for (Int16 spectrumIdx = 0; spectrumIdx < m_nSpectra; spectrumIdx++) {
    const CSpectrumSpectralAxis &spectrumLambda =
        Context.getSpectra()[spectrumIdx]->GetSpectralAxis();
    m_lambdaRanges[spectrumIdx]->getClosedIntervalIndices(
        spectrumLambda.GetSamplesVector(), m_firstPixelIdxInRange[spectrumIdx],
        m_lastPixelIdxInRange[spectrumIdx]);
    m_nPixels[spectrumIdx] = m_lastPixelIdxInRange[spectrumIdx] -
                             m_firstPixelIdxInRange[spectrumIdx] + 1;
  }
  m_igmCorrectionMeiksin = Context.getFluxCorrectionMeiksin();
  m_ismCorrectionCalzetti = Context.getFluxCorrectionCalzetti();

  initIgmIsm(true, true);
}

void COperatorPowerLaw::initIgmIsm(bool opt_extinction, bool opt_dustFitting) {
  if (opt_extinction) {
    m_nIgmCurves = m_igmCorrectionMeiksin->getIdxCount();
  } else {
    m_nIgmCurves = 1;
  }
  if (opt_dustFitting) {
    m_nIsmCurves = m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs();
  } else {
    m_nIsmCurves = 1;
  }
}

TPowerLawResult COperatorPowerLaw::BasicFit(Float64 redshift,
                                            bool opt_extinction,
                                            bool opt_dustFitting,
                                            Float64 nullFluxThreshold,
                                            std::string method) {

  initIgmIsm(opt_extinction, opt_dustFitting);
  T3DCurve emittedCurve = computeEmittedCurve(redshift, nullFluxThreshold,
                                              opt_extinction, opt_dustFitting);
  T3DCurve lnCurve = computeLnCurve(emittedCurve);

  // Step 3. Compute power law coefs and chi2
  T2DPowerLawCoefsPair coefs = powerLawCoefs3D(lnCurve, method);
  std::vector<TFloat64List> ChiSquareInterm;
  TChi2Result chi2Result =
      findMinChi2OnIgmIsm(emittedCurve, coefs, ChiSquareInterm);
  // Step 4. Creates result
  TPowerLawResult result;
  result.chiSquare = chi2Result.chi2;
  result.coefs = coefs[chi2Result.igmIdx]
                      [chi2Result.ismIdx]; // TODO put indexes there later
  if (opt_dustFitting)
    result.ebmvCoef = m_ismCorrectionCalzetti->GetEbmvValue(chi2Result.ismIdx);
  if (opt_extinction)
    result.meiksinIdx = chi2Result.igmIdx;
  // Question: j'ai l'impression que ChiSquareInterm et IsmCalzettiCoeffInterm
  // ne sont utilisés null part. OK pour les supprimer du TPowerLawResult ?
  return result;
};

T3DCurve COperatorPowerLaw::computeLnCurve(T3DCurve const &emittedCurve) const {
  T3DCurve lnCurve = emittedCurve;
  lnCurve.setLambda(lnLambda(emittedCurve.getLambda()));
  for (Int32 igmIdx = 0; igmIdx < m_nIgmCurves; igmIdx++) {
    for (Int32 ismIdx = 0; ismIdx < m_nIsmCurves; ismIdx++) {
      for (Int32 pixelIdx = 0; pixelIdx < emittedCurve.size(); pixelIdx++) {
        Float64 fluxValue = NAN;
        if (!emittedCurve.getIsExtinctedAt(igmIdx, ismIdx, pixelIdx))
          fluxValue =
              std::log(emittedCurve.getFluxAt(igmIdx, ismIdx, pixelIdx));
        lnCurve.setFluxAt(igmIdx, ismIdx, pixelIdx, fluxValue);
      }
    }
  }

  return lnCurve;
}

TChi2Result COperatorPowerLaw::findMinChi2OnIgmIsm(
    T3DCurve const &curve3D, T2DPowerLawCoefsPair const &coefs,
    std::vector<TFloat64List> &ChiSquareInterm) {
  ChiSquareInterm = computeChi2(curve3D, coefs);
  TInt32Pair minChi2Idxs = find2DVectorMinIndexes(ChiSquareInterm);
  return {minChi2Idxs.first, minChi2Idxs.second,
          ChiSquareInterm[minChi2Idxs.first][minChi2Idxs.second]};
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

T2DList<Float64>
COperatorPowerLaw::computeChi2(T3DCurve const &curve3D,
                               T2DPowerLawCoefsPair const &coefs) {
  T2DList<Float64> chi2(m_nIgmCurves, TList<Float64>(m_nIsmCurves, 0));
  Float64 diff;

  for (Int16 igmIdx = 0; igmIdx < m_nIgmCurves; igmIdx++) {
    for (Int16 ismIdx = 0; ismIdx < m_nIsmCurves; ismIdx++) {
      for (Int32 pixelIdx = 0; pixelIdx < m_nPixels[0]; pixelIdx++) {
        if (curve3D.pixelIsChi2Valid(pixelIdx)) {
          diff = curve3D.getFluxAt(igmIdx, ismIdx, pixelIdx) -
                 theoreticalFluxAtLambda(
                     TPowerLawCoefsPair(coefs[igmIdx][ismIdx].first,
                                        coefs[igmIdx][ismIdx].second),
                     curve3D.getLambdaAt(pixelIdx));
          diff = diff / curve3D.getFluxErrorAt(igmIdx, ismIdx, pixelIdx);
          chi2[igmIdx][ismIdx] += diff * diff;
        }
      }
    }
  }
  return chi2;
}

TAxisSampleList
COperatorPowerLaw::lnLambda(TAxisSampleList const &lambda) const {
  TAxisSampleList lnLambda(lambda.size());
  std::transform(lambda.begin(), lambda.end(), lnLambda.begin(),
                 [](float value) { return std::log(value); });
  return lnLambda;
}

std::shared_ptr<COperatorResult>
COperatorPowerLaw::Compute(bool opt_extinction, bool opt_dustFitting,
                           Float64 nullFluxThreshold, std::string method) {

  // Creates power law result
  // TODO check that m_redshifts is set somewhere before
  std::shared_ptr<CPowerLawResult> result =
      std::make_shared<CPowerLawResult>(m_redshifts.size());

  result->Redshifts = m_redshifts;

  for (Int32 zIdx = 0; zIdx < m_redshifts.size(); zIdx++) {
    Float64 redshift = result->Redshifts[zIdx];
    TPowerLawResult result_z = BasicFit(
        redshift, opt_extinction, opt_dustFitting, nullFluxThreshold, method);

    result->set_at_redshift(zIdx, std::move(result_z));
  }

  result->CstLog = EstimateLikelihoodCstLog();

  return result;
}

void COperatorPowerLaw::addTooFewSamplesWarning(Int32 N, Int16 igmIdx,
                                                Int16 ismIdx,
                                                const char *funcName) const {
  Flag.warning(WarningCode::FORCED_POWERLAW_TO_ZERO,
               Formatter() << "COperatorPowerLaw::" << funcName << ": only "
                           << N << " < " << m_nLogSamplesMin
                           << " samples with significant flux values. Power "
                              "law coefs are forced to zero. igmIdx = "
                           << igmIdx << ", "
                           << "ismIdx = " << ismIdx);
}

T2DPowerLawCoefsPair
COperatorPowerLaw::powerLawCoefs3D(T3DCurve const &lnCurve,
                                   std::string method) const {

  T2DPowerLawCoefsPair powerLawsCoefs(
      m_nIgmCurves, std::vector<TPowerLawCoefsPair>(m_nIsmCurves));

  Float64 lnxc = std::log(m_lambdaCut);

  // Computes power law for each igm / ism depending on the method selected
  // above
  for (Int32 igmIdx = 0; igmIdx < m_nIgmCurves; igmIdx++) {
    for (Int32 ismIdx = 0; ismIdx < m_nIsmCurves; ismIdx++) {
      Int32 N1 = 0;
      Int32 N2 = 0;
      for (Int32 pixelIdx = 0; pixelIdx < lnCurve.size(); pixelIdx++) {
        if (lnCurve.pixelIsCoefValid(igmIdx, ismIdx, pixelIdx)) {
          if (lnCurve.getLambdaAt(pixelIdx) < lnxc) {
            N1 += 1;
          } else {
            N2 += 1;
          }
        }
      }
      Int32 N = N1 + N2;

      if (N < m_nLogSamplesMin) {
        addTooFewSamplesWarning(N, igmIdx, ismIdx, __func__);
        powerLawsCoefs[igmIdx][ismIdx] = {{0, 0}, {0, 0}};
      } else {
        TCurve curve = lnCurve.toCoefCurve(igmIdx, ismIdx);
        if (method == "full") {
          powerLawsCoefs[igmIdx][ismIdx] =
              computeFullPowerLawCoefs(N1, N2, curve);
        } else if (method == "simple") {
          TPowerLawCoefs coefs = computeSimplePowerLawCoefs(curve);
          powerLawsCoefs[igmIdx][ismIdx] = {coefs, coefs};
        } else if (method == "simpleWeighted") {
          TPowerLawCoefs coefs = compute2PassSimplePowerLawCoefs(curve);

          powerLawsCoefs[igmIdx][ismIdx] = {coefs, coefs};
        } else {
          THROWG(ErrorCode::INTERNAL_ERROR,
                 Formatter() << "Unexpected method " << method);
        }
      }
    }
  }
  return powerLawsCoefs;
}

TPowerLawCoefsPair
COperatorPowerLaw::computeFullPowerLawCoefs(Int32 N1, Int32 N2,
                                            TCurve const &lnCurve) const {
  // If one part of the curve has too little samples, calculate the coefs
  // with the other part, and set the same coefs on the small part
  Float64 lnxc = std::log(m_lambdaCut);

  TPowerLawCoefsPair powerLawsCoefs;
  TCurve lnPartCurve;
  lnPartCurve.reserve(N1 + N2);
  if (N1 < m_nLogSamplesMin) {
    for (Int32 pixelIdx = 0; pixelIdx < lnCurve.size(); pixelIdx++) {
      if (lnCurve.getLambdaAt(pixelIdx) > lnxc) {
        lnPartCurve.push_back(lnCurve.get_at_index(pixelIdx));
      }
    }
    TPowerLawCoefs coefs = compute2PassSimplePowerLawCoefs(lnPartCurve);
    powerLawsCoefs = {coefs, coefs};
  } else if (N2 < m_nLogSamplesMin) {
    for (Int32 pixelIdx = 0; pixelIdx < lnCurve.size(); pixelIdx++) {
      if (lnCurve.getLambdaAt(pixelIdx) < lnxc) {
        lnPartCurve.push_back(lnCurve.get_at_index(pixelIdx));
      }
    }
    TPowerLawCoefs coefs = compute2PassSimplePowerLawCoefs(lnPartCurve);
    powerLawsCoefs = {coefs, coefs};
  } else {
    powerLawsCoefs = compute2PassDoublePowerLawCoefs(lnCurve);
  }

  return powerLawsCoefs;
};

TPowerLawCoefs COperatorPowerLaw::compute2PassSimplePowerLawCoefs(
    TCurve const &lnCurves) const {
  TPowerLawCoefs coefsFirstEstim = computeSimplePowerLawCoefs(lnCurves);
  // TODO add an error here if is inf
  // TODO add here a condition on result ?
  TPowerLawCoefs coefsSecondEstim =
      computeSimplePowerLawCoefs(lnCurves, coefsFirstEstim);
  return coefsSecondEstim;
}

TPowerLawCoefs COperatorPowerLaw::computeSimplePowerLawCoefs(
    TCurve const &lnCurve,
    std::optional<TPowerLawCoefs> const &coefsFirstEstim) const {

  Float64 SX = 0;
  Float64 SY = 0;
  Float64 SXY = 0;
  Float64 SX2 = 0;
  Float64 n = 0;

  for (Int32 pixelIdx = 0; pixelIdx < lnCurve.size(); pixelIdx++) {
    Float64 w = 1;
    if (coefsFirstEstim.has_value()) {
      Float64 estimatedFlux = computeEstimatedFlux(
          coefsFirstEstim.value(), lnCurve.getLambdaAt(pixelIdx));
      w = lnCurve.getFluxErrorAt(pixelIdx) / estimatedFlux;
    }
    Float64 X = lnCurve.getLambdaAt(pixelIdx) / w;
    Float64 Y = lnCurve.getFluxAt(pixelIdx) / w;
    SX += X / w;
    SY += Y / w;
    SXY += X * Y;
    SX2 += X * X;
    n += 1 / (w * w);
  }

  Float64 denomInv = 1 / (n * SX2 - SX * SX);
  Float64 b = (n * SXY - SX * SY) * denomInv;
  Float64 a = std::exp((SY - b * SX) / n);

  Float64 sigmalna = sqrt(SX2 * denomInv);
  Float64 stda = a * sigmalna;
  Float64 sigmab = sqrt(n * denomInv);

  return {a, b, stda, sigmab};
}

TPowerLawCoefsPair COperatorPowerLaw::compute2PassDoublePowerLawCoefs(
    TCurve const &lnCurves) const {
  // Make a first calculation of power law coefficients without taking into
  // account the noise
  TPowerLawCoefsPair coefsFirstEstim = computeDoublePowerLawCoefs(lnCurves);
  TPowerLawCoefsPair coefsSecondEstim =
      computeDoublePowerLawCoefs(lnCurves, coefsFirstEstim);
  return coefsSecondEstim;
}

TPowerLawCoefsPair COperatorPowerLaw::computeDoublePowerLawCoefs(
    TCurve const &lnCurve,
    std::optional<TPowerLawCoefsPair> const &coefsFirstEstim) const {
  // NB : coefsFirstEstim are used to calculate the weights. If absent
  // calculations are made without taking error into account

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
  for (Int32 pixelIdx = 0; pixelIdx < lnCurve.size(); pixelIdx++) {
    Float64 xi = lnCurve.getLambdaAt(pixelIdx);
    Float64 yi = lnCurve.getFluxAt(pixelIdx);
    Float64 wi = 1;

    if (xi < xc) {
      if (coefsFirstEstim.has_value()) {
        Float64 estimatedFlux =
            computeEstimatedFlux(coefsFirstEstim.value().first, xi);
        wi = (estimatedFlux * estimatedFlux) /
             (lnCurve.getFluxErrorAt(pixelIdx) *
              lnCurve.getFluxErrorAt(pixelIdx));
      }
      n1 += wi;
      sx1 += xi * wi;
      sxx1 += xi * xi * wi;
      sy1 += yi * wi;
      sxy1 += xi * yi * wi;
      N1 += 1;
    } else {
      if (coefsFirstEstim.has_value()) {
        Float64 estimatedFlux =
            computeEstimatedFlux(coefsFirstEstim.value().second, xi);
        wi = (estimatedFlux * estimatedFlux) /
             (lnCurve.getFluxErrorAt(pixelIdx) *
              lnCurve.getFluxErrorAt(pixelIdx));
      }
      n2 += wi;
      sx2 += xi * wi;
      sxx2 += xi * xi * wi;
      sy2 += yi * wi;
      sxy2 += xi * yi * wi;
      sx2mc2 += wi * (xc - xi) * (xc - xi);
      N2 += 1;
    }
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

  if (std::abs(m.determinant()) < DBL_MIN)
    THROWG(ErrorCode::INTERNAL_ERROR,
           "Cannot calculate power law coefs: division by zero");

  Eigen::Matrix3d mInv = m.inverse();
  Eigen::Vector3d v(v1, v2, v3);

  Eigen::Vector3d theta = mInv * v;

  Float64 a1 = std::exp(theta(0));
  Float64 b1 = theta(1);
  Float64 b2 = theta(2);
  Float64 a2 = a1 * std::pow(std::exp(xc), b1 - b2);

  // Calculates var / covar
  Float64 varlna2 = mInv(0, 0) +
                    xc * xc * (mInv(1, 1) + mInv(2, 2) - 2 * mInv(1, 2)) +
                    2 * xc * (mInv(0, 1) - mInv(0, 2));
  Float64 sigmaa1 = a1 * std::sqrt(mInv(0, 0));
  Float64 sigmab1 = std::sqrt(mInv(1, 1));
  Float64 sigmaa2 = a2 * std::sqrt(varlna2);
  Float64 sigmab2 = std::sqrt(mInv(2, 2));

  return {{a1, b1, sigmaa1, sigmab1}, {a2, b2, sigmaa2, sigmab2}};
}

T3DCurve COperatorPowerLaw::initializeFluxCurve(Float64 redshift,
                                                Float64 nullFluxThreshold) {
  // In order to take into account ism, igm we initialize a template with a flux
  // at 1 we then divide the initial flux by the template flux value once ism
  // igm has been applied

  TList<Float64> spectrumLambda;
  TList<Float64> spectrumFlux;
  TList<Float64> spectrumFluxError;
  // Concatenates all curves
  // NB: at the end, lambda is not ordered anymore
  for (size_t spectrumIdx = 0; spectrumIdx < m_nSpectra; spectrumIdx++) {
    TList<Float64> tmpLambda = m_spectra[spectrumIdx]
                                   ->GetSpectralAxis()
                                   .extract(m_firstPixelIdxInRange[spectrumIdx],
                                            m_lastPixelIdxInRange[spectrumIdx])
                                   .GetSamplesVector();
    spectrumLambda.insert(spectrumLambda.end(),
                          std::make_move_iterator(tmpLambda.begin()),
                          std::make_move_iterator(tmpLambda.end()));

    TList<Float64> tmpFlux = m_spectra[spectrumIdx]
                                 ->GetFluxAxis()
                                 .extract(m_firstPixelIdxInRange[spectrumIdx],
                                          m_lastPixelIdxInRange[spectrumIdx])
                                 .GetSamplesVector();
    spectrumFlux.insert(spectrumFlux.end(),
                        std::make_move_iterator(tmpFlux.begin()),
                        std::make_move_iterator(tmpFlux.end()));

    TList<Float64> tmpError = m_spectra[spectrumIdx]
                                  ->GetFluxAxis()
                                  .GetError()
                                  .extract(m_firstPixelIdxInRange[spectrumIdx],
                                           m_lastPixelIdxInRange[spectrumIdx])
                                  .GetSamplesVector();
    spectrumFluxError.insert(spectrumFluxError.end(),
                             std::make_move_iterator(tmpError.begin()),
                             std::make_move_iterator(tmpError.end()));
  }

  CSpectrumSpectralAxis spectrumLambdaAxis(std::move(spectrumLambda));

  TList<uint8_t> maskedPixels =
      m_maskBuilder
          ->getMask(spectrumLambdaAxis, spectrumLambdaAxis.GetLambdaRange(),
                    redshift)
          .getMaskList();
  ;

  // Step 2. Transform spectrum data

  // Initializes usable pixels
  spectrumLambdaAxis.blueShiftInplace(redshift);
  TBoolList snrCompliantPixels = computeSNRCompliantPixels(
      spectrumFlux, spectrumFluxError, nullFluxThreshold);
  TCurve fluxCurve1D({spectrumLambdaAxis.GetSamplesVector(),
                      std::move(spectrumFlux), std::move(spectrumFluxError)});
  fluxCurve1D.sort();

  T3DCurve fluxCurve3D(m_nIgmCurves, m_nIsmCurves);
  fluxCurve3D.copyCurveAtAllIgmIsm(fluxCurve1D);
  fluxCurve3D.setIsSnrCompliant(std::move(snrCompliantPixels));
  fluxCurve3D.setMask(std::move(maskedPixels));

  return fluxCurve3D;
}

T3DList<Float64> COperatorPowerLaw::computeIsmIgmCorrections(
    Float64 redshift, CSpectrumSpectralAxis const &spectrumLambdaRest,
    bool opt_extinction, bool opt_dustFitting) const {
  // In order to access ism igm coefs, we initialize a template with a flux
  // at 1, and apply ism/igm on it
  Int32 n = spectrumLambdaRest.GetSamplesCount();
  CTemplate templateForCoefs(
      "", "", spectrumLambdaRest,
      std::vector<Float64>(spectrumLambdaRest.GetSamplesCount(), 1));
  templateForCoefs.InitIsmIgmConfig(redshift, m_ismCorrectionCalzetti,
                                    m_igmCorrectionMeiksin);

  T3DList<Float64> correctionCoefs(
      m_nIgmCurves,
      std::vector<std::vector<Float64>>(
          m_nIsmCurves,
          std::vector<Float64>(spectrumLambdaRest.GetSamplesCount(), NAN)));

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

TList<Float64> COperatorPowerLaw::computeIsmIgmCorrection(
    Float64 redshift, CSpectrumSpectralAxis const &spectrumLambdaRest,
    Int16 igmIdx, Float64 ismCoef) const {
  // In order to access ism igm coefs, we initialize a template with a flux
  // at 1, and apply ism/igm on it
  Int32 n = spectrumLambdaRest.GetSamplesCount();
  CTemplate templateForCoefs(
      "", "", spectrumLambdaRest,
      std::vector<Float64>(spectrumLambdaRest.GetSamplesCount(), 1));
  templateForCoefs.InitIsmIgmConfig(redshift, m_ismCorrectionCalzetti,
                                    m_igmCorrectionMeiksin);

  TList<Float64> correctionCoefs(spectrumLambdaRest.GetSamplesCount(), NAN);

  Float64 igmCorrection = 1;
  Float64 ismCorrection = 1;
  if (igmIdx > -1) {
    templateForCoefs.ApplyMeiksinCoeff(igmIdx);
  }
  if (ismCoef > 0) {
    Int32 ismIdx = -1;
    // TODO find what to take here
    ismIdx = m_ismCorrectionCalzetti->GetEbmvIndex(ismCoef);
    templateForCoefs.ApplyDustCoeff(ismIdx);
  }
  correctionCoefs = templateForCoefs.GetFluxAxis().GetSamplesVector();
  return correctionCoefs;
}

T3DCurve COperatorPowerLaw::computeEmittedCurve(Float64 redshift,
                                                Float64 nullFluxThreshold,
                                                bool opt_extinction,
                                                bool opt_dustFitting) {
  // In order to take into account ism, igm we initialize a template with a flux
  // at 1 we then divide the initial flux by the template flux value once ism
  // igm has been applied
  // Set m_ismIgmCorrections

  T3DCurve fluxCurve = initializeFluxCurve(redshift, nullFluxThreshold);

  T3DList<bool> isExtincted(
      fluxCurve.getNIgm(),
      T2DList<bool>(fluxCurve.getNIsm(), TList<bool>(fluxCurve.size(), false)));
  if (opt_extinction || opt_dustFitting) {
    T3DList<Float64> m_ismIgmCorrections = computeIsmIgmCorrections(
        redshift, fluxCurve.getLambda(), opt_extinction, opt_dustFitting);
    for (Int32 igmIdx = 0; igmIdx < m_nIgmCurves; igmIdx++) {
      for (Int32 ismIdx = 0; ismIdx < m_nIsmCurves; ismIdx++) {
        for (Int32 pixelIdx = 0; pixelIdx < fluxCurve.size(); pixelIdx++) {
          Float64 correctionCoef =
              m_ismIgmCorrections[igmIdx][ismIdx][pixelIdx];
          if (correctionCoef < DBL_MIN) {
            isExtincted[igmIdx][ismIdx][pixelIdx] = true;
            continue;
          }
          fluxCurve.setFluxAt(igmIdx, ismIdx, pixelIdx,
                              fluxCurve.getFluxAt(igmIdx, ismIdx, pixelIdx) /
                                  correctionCoef);
          fluxCurve.setFluxErrorAt(
              igmIdx, ismIdx, pixelIdx,
              fluxCurve.getFluxErrorAt(igmIdx, ismIdx, pixelIdx) /
                  correctionCoef);
        }
      }
    }
  }
  fluxCurve.setIsExtincted(std::move(isExtincted));
  return fluxCurve;
}

TBoolList COperatorPowerLaw::computeSNRCompliantPixels(
    TFloat64List const &spectrumFlux, TFloat64List const &spectrumFluxError,
    Float64 nullFluxThreshold) const {
  // Masks pixels with insufficient SNR. In place modification of
  // curve.pixelsToUse
  TBoolList snrCompliant = std::vector<bool>(spectrumFlux.size(), true);
  for (Int32 pixelIdx = 0; pixelIdx < spectrumFlux.size(); pixelIdx++) {
    // Checks that signal / noise ratio is sufficient (> nullFluxThreshold)
    if (spectrumFlux[pixelIdx] <
        nullFluxThreshold * spectrumFluxError[pixelIdx]) {
      snrCompliant[pixelIdx] = false;
    }
  }
  return snrCompliant;
}
Float64 COperatorPowerLaw::computeEstimatedFlux(TPowerLawCoefs const &coefs,
                                                Float64 x) const {
  // with x = ln(lambda)
  return coefs.a * std::pow(std::exp(x), coefs.b);
}

void COperatorPowerLaw::ComputeSpectrumModel(
    const std::shared_ptr<CContinuumModelSolution> &continuum, Int32 spcIndex,
    const std::shared_ptr<CModelSpectrumResult> &models) {

  // TODO here create a model with
  // lambda, flux created from input coefs
  auto lambdaObsAxis = m_spectra[spcIndex]->GetSpectralAxis();
  auto lambdaObs = lambdaObsAxis.GetSamplesVector();
  auto lambdaRestAxis = lambdaObsAxis.blueShift(continuum->redshift);

  // Calculates ism igm corrections for given lambdaRest / redshift
  TList<Float64> correctionCoefs =
      computeIsmIgmCorrection(continuum->redshift, lambdaRestAxis,
                              continuum->meiksinIdx, continuum->ebmvCoef);
  // TODO above check which one is igm / ism
  // Use lambda rest to calculate flux rest and apply ism igm on flux rest to
  // get flux obs
  TList<Float64> fluxObs(lambdaObs.size(), NAN);
  for (size_t pixelIdx = 0; pixelIdx < lambdaObs.size(); pixelIdx++) {
    fluxObs[pixelIdx] =
        theoreticalFluxAtLambda(
            {{continuum->a1, continuum->b1}, {continuum->a2, continuum->b2}},
            lambdaRestAxis.GetSamplesVector()[pixelIdx]) *
        correctionCoefs[pixelIdx];
  }

  Float64 overlapFraction = 0.0;
  TFloat64Range currentRange;

  models->addModel(CSpectrum(lambdaObsAxis, CSpectrumFluxAxis(fluxObs)),
                   m_spectra[spcIndex]->getObsID());
  return;
}
