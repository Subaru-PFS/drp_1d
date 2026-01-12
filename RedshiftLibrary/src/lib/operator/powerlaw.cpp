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
#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/common/vectorOperations.h"
#include "RedshiftLibrary/operator/continuumfitting.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/operator/powerlawresult.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/statistics/fitquality.h"

#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include <cmath>
#include <utility>

using namespace NSEpic;
using namespace std;

COperatorPowerLaw::COperatorPowerLaw(const TFloat64List &redshifts,
                                     Float64 lambdaCut)
    : COperatorContinuumFitting(), m_lambdaCut(lambdaCut) {

  // Gets pixels of interest and sets elements size
  SetRedshifts(redshifts);
  m_spectra = Context.getSpectra();
  m_nSpectra = m_spectra.size();
  m_nPixels.resize(m_nSpectra);
  m_kStart.resize(m_nSpectra);
  m_kEnd.resize(m_nSpectra);
  for (Int32 spectrumIdx = 0; spectrumIdx < m_nSpectra; spectrumIdx++) {
    const CSpectrumSpectralAxis &spectrumLambda =
        Context.getSpectra()[spectrumIdx]->GetSpectralAxis();
    std::tie(m_kStart[spectrumIdx], m_kEnd[spectrumIdx]) =
        m_lambdaRanges[spectrumIdx]->getClosestInnerIndices(
            spectrumLambda.GetSamplesVector());
    m_nPixels[spectrumIdx] = m_kEnd[spectrumIdx] - m_kStart[spectrumIdx] + 1;
  }
  m_igmCorrectionMeiksin = Context.getFluxCorrectionMeiksin();
  m_ismCorrectionCalzetti = Context.getFluxCorrectionCalzetti();
  m_powerCoefsLimits.first.min =
      Context.GetParameterStore()->GetScoped<Float64>(
          "powerLaw.firstPowerCoefMin");
  m_powerCoefsLimits.first.max =
      Context.GetParameterStore()->GetScoped<Float64>(
          "powerLaw.firstPowerCoefMax");
  m_powerCoefsLimits.second.min =
      Context.GetParameterStore()->GetScoped<Float64>(
          "powerLaw.secondPowerCoefMin");
  m_powerCoefsLimits.second.max =
      Context.GetParameterStore()->GetScoped<Float64>(
          "powerLaw.secondPowerCoefMax");
}

void COperatorPowerLaw::initIgmIsm(bool opt_extinction, bool opt_dustFitting,
                                   Int32 FitEbmvIdx, Int32 FitMeiksinIdx) {
  TIgmIsmIdxs igmIsmIdxs = Context.GetIsmIgmIdxList(
      opt_extinction, opt_dustFitting, FitEbmvIdx, FitMeiksinIdx);

  m_igmIdxList = igmIsmIdxs.igmIdxs;
  m_nIgmCurves = m_igmIdxList.size();
  m_ismIdxList = igmIsmIdxs.ismIdxs;
  m_nIsmCurves = m_ismIdxList.size();
}

TPowerLawResult COperatorPowerLaw::BasicFit(Float64 redshift,
                                            bool opt_extinction,
                                            bool opt_dustFitting,
                                            Float64 nullFluxThreshold) {

  TCurve curve = initializeFluxCurve(redshift, nullFluxThreshold);

  // handle null or negative spectrum
  auto N = std::count_if(boost::counting_iterator<Int32>(0),
                         boost::counting_iterator<Int32>(curve.size()),
                         [&curve](Int32 pixelIdx) {
                           return curve.pixelIsChi2AndSNRValid(pixelIdx);
                         });
  TPowerLawResult result;
  Int32 nUnmasked = 0;
  if (N < m_nSamplesMinForContinuumFit) {
    auto const flux = curve.computeUnmaskedFlux();
    nUnmasked = flux.size();
    if (nUnmasked < m_nSamplesMinForContinuumFit) {
      result.coefs = NULL_COEFS_PAIR;
      result.chiSquare = INFINITY;
      return result;
    }
    // If the number of valid pixels is too low, set igm / ism indexes to 0 and
    // constant power law
    auto const error = curve.computeUnmaskedFluxError();
    result.coefs = computeConstantLawCoefs(flux, error);
    T2DPowerLawCoefsPair coefs(1, TList<TPowerLawCoefsPair>(1, result.coefs));
    // Create a temporary 3D curve to compute chi2
    auto curve3D = T3DCurve(std::move(curve));
    auto const chi2 = computeChi2(curve3D, coefs);
    curve = TCurve(std::move(curve3D));
    result.chiSquare = chi2[0][0];
    if (opt_extinction)
      result.meiksinIdx = undefIdx;
    if (opt_dustFitting)
      result.ebmvCoef = 0.0;
  } else {
    T3DCurve emittedCurve = computeEmittedCurve(
        redshift, opt_extinction, opt_dustFitting, std::move(curve));

    // Step 3. Compute power law coefs and chi2
    T2DPowerLawCoefsPair coefs = powerLawCoefs3D(emittedCurve);
    TChi2Result chi2Result = findMinChi2OnIgmIsm(emittedCurve, coefs);
    // Step 4. Creates result
    result.chiSquare = chi2Result.chi2;
    result.coefs = coefs[chi2Result.igmIdx][chi2Result.ismIdx];

    // Adds number of pixels ised for the continuum fit info
    auto const igmIdx = chi2Result.igmIdx;
    auto const ismIdx = chi2Result.ismIdx;
    nUnmasked = std::count_if(
        boost::counting_iterator<Int32>(0),
        boost::counting_iterator<Int32>(emittedCurve.size()),
        [&emittedCurve, igmIdx, ismIdx](Int32 pixelIdx) {
          return emittedCurve.pixelIsCoefValid(igmIdx, ismIdx, pixelIdx);
        });

    if (opt_extinction)
      result.meiksinIdx = m_igmIdxList[chi2Result.igmIdx];
    if (opt_dustFitting)
      result.ebmvCoef = m_ismCorrectionCalzetti->GetEbmvValue(
          m_ismIdxList[chi2Result.ismIdx]);
    curve =
        TCurve(std::move(emittedCurve), chi2Result.igmIdx, chi2Result.ismIdx);
  }

  // Step 5. compute fit quality from residuals
  // compute power law model without isgm/igm since flux & error has been
  // inverse corrected
  auto modelFlux =
      computeModelFlux(CSpectrumSpectralAxis(curve.computeUnmaskedLambda()),
                       redshift, undefIdx, 0.0, result.coefs);
  T2DPowerLawCoefsPair coefs(1, TList<TPowerLawCoefsPair>(1, result.coefs));
  auto flux = curve.computeUnmaskedFlux();
  auto error = curve.computeUnmaskedFluxError();
  result.fitQuality = NSFitQuality::computeFitQuality(
      std::move(flux), std::move(modelFlux), std::move(error), NAN, undefIdx,
      nUnmasked);
  return result;
};
void COperatorPowerLaw::limitCoefs(TPowerLawCoefsPair &coefs) {
  const auto newb1 = limitCoef(coefs.first.b, m_powerCoefsLimits.first);
  const auto newb2 = limitCoef(coefs.second.b, m_powerCoefsLimits.second);
  auto b1_limited = newb1 != coefs.first.b;
  auto b2_limited = newb2 != coefs.second.b;

  if (b1_limited && b2_limited) {
    coefs = computeDoublePowerLawCoefs_b1_b2_fixed(newb1, newb2);
  } else if (b1_limited) {
    coefs = computeDoublePowerLawCoefs_b1_fixed(newb1);
  } else if (b2_limited) {
    coefs = computeDoublePowerLawCoefs_b2_fixed(newb2);
  }
}

Float64 COperatorPowerLaw::limitCoef(Float64 coef,
                                     TPowerCoefLimits limits) const {
  const auto newvalue = std::max(limits.min, std::min(limits.max, coef));
  return newvalue;
}

T3DCurve COperatorPowerLaw::computeLnCurve(T3DCurve const &emittedCurve) const {
  T3DCurve lnCurve = emittedCurve;
  lnCurve.setLambda(lnLambda(emittedCurve.getLambda()));
  for (Int32 igmIdx = 0; igmIdx < m_nIgmCurves; igmIdx++) {
    for (Int32 ismIdx = 0; ismIdx < m_nIsmCurves; ismIdx++) {
      for (Int32 pixelIdx = 0; pixelIdx < emittedCurve.size(); pixelIdx++) {
        Float64 fluxValue = NAN;
        if (emittedCurve.pixelIsCoefValid(igmIdx, ismIdx, pixelIdx))
          fluxValue =
              std::log(emittedCurve.getFluxAt(igmIdx, ismIdx, pixelIdx));
        lnCurve.setFluxAt(igmIdx, ismIdx, pixelIdx, fluxValue);
      }
    }
  }

  return lnCurve;
}

TChi2Result
COperatorPowerLaw::findMinChi2OnIgmIsm(T3DCurve const &curve3D,
                                       T2DPowerLawCoefsPair const &coefs) {
  T2DList<Float64> chiSquareInterm = computeChi2(curve3D, coefs);
  TInt32Pair minChi2Idxs = NSVectorOp::find2DVectorMinIndexes(chiSquareInterm);
  return {minChi2Idxs.first, minChi2Idxs.second,
          chiSquareInterm[minChi2Idxs.first][minChi2Idxs.second]};
}

T2DList<Float64>
COperatorPowerLaw::computeChi2(T3DCurve const &curve3D,
                               T2DPowerLawCoefsPair const &coefs) {
  std::function<bool(Int32)> considerPixel;
  considerPixel = [&curve3D](Int32 pixelIdx) {
    return curve3D.pixelIsChi2Valid(pixelIdx);
  };
  Int32 nIgmCurves = curve3D.getNIgm();
  Int32 nIsmCurves = curve3D.getNIsm();
  T2DList<Float64> chi2_all(nIgmCurves, TList<Float64>(nIsmCurves, INFINITY));

  for (Int32 igmIdx = 0; igmIdx < nIgmCurves; igmIdx++) {
    for (Int32 ismIdx = 0; ismIdx < nIsmCurves; ismIdx++) {
      Float64 chi2 = 0.0;
      Int32 nPixels = 0;
      for (Int32 pixelIdx = 0; pixelIdx < m_nPixels[0]; pixelIdx++) {
        if (considerPixel(pixelIdx)) {
          Float64 theoreticalFlux =
              curve3D.getIsExtinctedAt(igmIdx, ismIdx, pixelIdx)
                  ? 0
                  : computeDoublePowerLaw(coefs[igmIdx][ismIdx],
                                          curve3D.getLambdaAt(pixelIdx));
          Float64 diff =
              curve3D.getFluxAt(igmIdx, ismIdx, pixelIdx) - theoreticalFlux;
          diff = diff / curve3D.getFluxErrorAt(igmIdx, ismIdx, pixelIdx);
          chi2 += diff * diff;
          ++nPixels;
        }
      }
      if (nPixels > 0)
        chi2_all[igmIdx][ismIdx] = chi2;
    }
  }
  return chi2_all;
}

TAxisSampleList
COperatorPowerLaw::lnLambda(TAxisSampleList const &lambda) const {
  TAxisSampleList lnLambda(lambda.size());
  std::transform(lambda.begin(), lambda.end(), lnLambda.begin(),
                 [](float value) { return std::log(value); });
  return lnLambda;
}

std::shared_ptr<const COperatorResult>
COperatorPowerLaw::Compute(bool opt_extinction, bool opt_dustFitting,
                           Float64 nullFluxThreshold, Int32 FitEbmvIdx,
                           Int32 FitMeiksinIdx) {
  initIgmIsm(opt_extinction, opt_dustFitting, FitEbmvIdx, FitMeiksinIdx);

  // Creates power law result
  std::shared_ptr<CPowerLawResult> result =
      std::make_shared<CPowerLawResult>(m_redshifts.size());

  result->Redshifts = m_redshifts;
  for (Int32 zIdx = 0; zIdx < ssize(m_redshifts); zIdx++) {
    Float64 redshift = result->Redshifts[zIdx];
    TPowerLawResult result_z =
        BasicFit(redshift, opt_extinction, opt_dustFitting, nullFluxThreshold);

    result->set_at_redshift(zIdx, std::move(result_z));
  }

  result->CstLog = EstimateLikelihoodCstLog();

  return result;
}

void COperatorPowerLaw::addTooFewSamplesWarning(Int32 N, Int32 igmIdx,
                                                Int32 ismIdx,
                                                const char *funcName) const {
  Flag.warning(WarningCode::FORCED_POWERLAW_TO_ZERO,
               Formatter() << "COperatorPowerLaw::" << funcName << ": only "
                           << N << " < " << m_nSamplesMinForContinuumFit
                           << " samples with significant flux values. Power "
                              "law coefs are forced to zero. igmIdx = "
                           << igmIdx << ", "
                           << "ismIdx = " << ismIdx);
}

T2DPowerLawCoefsPair
COperatorPowerLaw::powerLawCoefs3D(T3DCurve const &emittedCurve) {

  T2DPowerLawCoefsPair powerLawsCoefs(
      emittedCurve.getNIgm(),
      std::vector<TPowerLawCoefsPair>(emittedCurve.getNIsm()));

  Float64 lnxc = std::log(m_lambdaCut);
  T3DCurve lnCurve = computeLnCurve(emittedCurve);

  // Computes power law for each igm / ism depending on the method selected
  // above
  for (Int32 igmIdx = 0; igmIdx < m_nIgmCurves; igmIdx++) {
    for (Int32 ismIdx = 0; ismIdx < m_nIsmCurves; ismIdx++) {
      auto const N = std::count_if(
          boost::counting_iterator<Int32>(0),
          boost::counting_iterator<Int32>(lnCurve.size()),
          [&lnCurve, igmIdx, ismIdx](Int32 pixelIdx) {
            return lnCurve.pixelIsCoefValid(igmIdx, ismIdx, pixelIdx);
          });
      auto const N1 = std::count_if(
          boost::counting_iterator<Int32>(0),
          boost::counting_iterator<Int32>(lnCurve.size()),
          [&lnCurve, igmIdx, ismIdx, lnxc](Int32 pixelIdx) {
            return lnCurve.pixelIsCoefValid(igmIdx, ismIdx, pixelIdx) &&
                   lnCurve.getLambdaAt(pixelIdx) < lnxc;
          });
      auto const N2 = N - N1;

      if (N < m_nSamplesMinForContinuumFit) {
        addTooFewSamplesWarning(N, igmIdx, ismIdx, __func__);
        powerLawsCoefs[igmIdx][ismIdx] = DEFAULT_COEFS_PAIR;
      } else {
        TCurve curve = lnCurve.toCoefCurve(igmIdx, ismIdx);
        powerLawsCoefs[igmIdx][ismIdx] =
            computeFullPowerLawCoefs(N1, N2, curve);
      }
    }
  }
  return powerLawsCoefs;
}

TPowerLawCoefsPair
COperatorPowerLaw::computeConstantLawCoefs(TFloat64List const &flux,
                                           TFloat64List const &error) const {
  // Computes a constant law. Use all unmasked pixels to compute the mean flux
  // (including the ones with low SNR)
  TFloat64List inverse_var(error.size());
  std::transform(error.cbegin(), error.cend(), inverse_var.begin(),
                 [](Float64 v) { return 1.0 / (v * v); });
  Float64 mean_amplitude = std::transform_reduce(flux.cbegin(), flux.cend(),
                                                 inverse_var.cbegin(), 0.0);
  Float64 const sum_inv_var =
      std::reduce(inverse_var.cbegin(), inverse_var.cend());
  mean_amplitude /= sum_inv_var;
  Float64 mean_amplitude_std = 1.0 / sqrt(sum_inv_var);
  TPowerLawCoefs coefs{mean_amplitude, 0.0, mean_amplitude_std, INFINITY};
  checkCoefsOrNull(coefs);
  return TPowerLawCoefsPair{coefs, coefs};
}

TPowerLawCoefsPair
COperatorPowerLaw::computeFullPowerLawCoefs(Int32 N1, Int32 N2,
                                            TCurve const &lnCurve) {
  // If one part of the curve has too little samples, calculate the coefs
  // with the other part, and set the same coefs on the small part
  Float64 lnxc = std::log(m_lambdaCut);

  TPowerLawCoefsPair powerLawsCoefs;
  TCurve lnPartCurve;
  lnPartCurve.reserve(N1 + N2);
  if (N1 < m_nSamplesMinForContinuumFit) {
    for (Int32 pixelIdx = 0; pixelIdx < lnCurve.size(); pixelIdx++) {
      if (lnCurve.getLambdaAt(pixelIdx) > lnxc) {
        lnPartCurve.push_back(lnCurve.get_at_index(pixelIdx));
      }
    }
    TPowerLawCoefs coefs = compute2PassSimplePowerLawCoefs(lnPartCurve);
    const auto newb = limitCoef(coefs.b, m_powerCoefsLimits.second);
    if (newb != coefs.b) {
      coefs = computeSimplePowerLawCoefs_b_fixed(newb);
    }
    powerLawsCoefs = {coefs, coefs};
  } else if (N2 < m_nSamplesMinForContinuumFit) {
    for (Int32 pixelIdx = 0; pixelIdx < lnCurve.size(); pixelIdx++) {
      if (lnCurve.getLambdaAt(pixelIdx) < lnxc) {
        lnPartCurve.push_back(lnCurve.get_at_index(pixelIdx));
      }
    }
    TPowerLawCoefs coefs = compute2PassSimplePowerLawCoefs(lnPartCurve);
    const auto newb = limitCoef(coefs.b, m_powerCoefsLimits.first);
    if (newb != coefs.b) {
      coefs = computeSimplePowerLawCoefs_b_fixed(newb);
    }
    powerLawsCoefs = {coefs, coefs};
  } else {
    powerLawsCoefs = compute2PassDoublePowerLawCoefs(lnCurve);
    limitCoefs(powerLawsCoefs);
  }

  checkCoefsOrNull(powerLawsCoefs);
  return powerLawsCoefs;
};

TPowerLawCoefs
COperatorPowerLaw::compute2PassSimplePowerLawCoefs(TCurve const &lnCurves) {

  TPowerLawCoefs coefs = computeSimplePowerLawCoefs(lnCurves);
  bool validCoefs = checkCoefsOrNull(coefs);
  if (validCoefs)
    coefs = computeSimplePowerLawCoefs(lnCurves, coefs);
  return coefs;
}

bool COperatorPowerLaw::checkCoefsOrNull(TPowerLawCoefs &coefs) const {
  if (coefs.a < DBL_MIN) {
    coefs = NULL_COEFS;
    return false;
  }
  return true;
}

bool COperatorPowerLaw::checkCoefsOrNull(TPowerLawCoefsPair &coefs) const {
  if (coefs.first.a < DBL_MIN || coefs.second.a < DBL_MIN) {
    coefs = NULL_COEFS_PAIR;
    return false;
  }
  return true;
}

void COperatorPowerLaw::updatePowerLawCalcStorageForSimple(
    TCurve const &lnCurve,
    std::optional<TPowerLawCoefs> const &coefsFirstEstim) {
  m_powerLawCalcStorage.xc = std::log(m_lambdaCut);

  m_powerLawCalcStorage.n1 = 0;
  m_powerLawCalcStorage.sx1 = 0;
  m_powerLawCalcStorage.sxx1 = 0;
  m_powerLawCalcStorage.sy1 = 0;
  m_powerLawCalcStorage.sxy1 = 0;
  m_powerLawCalcStorage.N1 = 0;
  for (Int32 pixelIdx = 0; pixelIdx < lnCurve.size(); pixelIdx++) {
    Float64 w = 1;
    if (coefsFirstEstim.has_value()) {
      Float64 estimatedFlux = computePowerLaw(
          coefsFirstEstim.value(), std::exp(lnCurve.getLambdaAt(pixelIdx)));
      w = estimatedFlux / lnCurve.getFluxErrorAt(pixelIdx);
    }
    Float64 X = lnCurve.getLambdaAt(pixelIdx);
    Float64 Y = lnCurve.getFluxAt(pixelIdx);
    Float64 w2 = w * w;
    m_powerLawCalcStorage.sx1 += X * w2;
    m_powerLawCalcStorage.sy1 += Y * w2;
    m_powerLawCalcStorage.sxy1 += X * Y * w2;
    m_powerLawCalcStorage.sxx1 += X * X * w2;
    m_powerLawCalcStorage.n1 += w2;
  }
}

TPowerLawCoefs COperatorPowerLaw::computeSimplePowerLawCoefs(
    TCurve const &lnCurve,
    std::optional<TPowerLawCoefs> const &coefsFirstEstim) {

  updatePowerLawCalcStorageForSimple(lnCurve, coefsFirstEstim);

  Float64 denomInv =
      1 / (m_powerLawCalcStorage.n1 * m_powerLawCalcStorage.sxx1 -
           m_powerLawCalcStorage.sx1 * m_powerLawCalcStorage.sx1);
  Float64 b = (m_powerLawCalcStorage.n1 * m_powerLawCalcStorage.sxy1 -
               m_powerLawCalcStorage.sx1 * m_powerLawCalcStorage.sy1) *
              denomInv;
  Float64 a =
      std::exp((m_powerLawCalcStorage.sy1 - b * m_powerLawCalcStorage.sx1) /
               m_powerLawCalcStorage.n1);

  Float64 sigmalna = sqrt(m_powerLawCalcStorage.sxx1 * denomInv);
  Float64 stda = a * sigmalna;
  Float64 stdb = sqrt(m_powerLawCalcStorage.n1 * denomInv);

  return {a, b, stda, stdb};
}

TPowerLawCoefs
COperatorPowerLaw::computeSimplePowerLawCoefs_b_fixed(Float64 b) const {

  Float64 a =
      std::exp(1 / m_powerLawCalcStorage.n1 *
               (m_powerLawCalcStorage.sy1 - b * m_powerLawCalcStorage.sx1));

  return {a, b, a * std::sqrt(1 / m_powerLawCalcStorage.n1), 0};
}

TPowerLawCoefsPair
COperatorPowerLaw::compute2PassDoublePowerLawCoefs(TCurve const &lnCurves) {
  // Make a first calculation of power law coefficients without taking into
  // account the noise
  TPowerLawCoefsPair coefs = computeDoublePowerLawCoefs(lnCurves);
  bool validCoefs = checkCoefsOrNull(coefs);
  if (validCoefs)
    coefs = computeDoublePowerLawCoefs(lnCurves, coefs);
  return coefs;
}

void COperatorPowerLaw::updatePowerLawCalcStorage(
    TCurve const &lnCurve,
    std::optional<TPowerLawCoefsPair> const &coefsFirstEstim) {
  m_powerLawCalcStorage.xc = std::log(m_lambdaCut);

  m_powerLawCalcStorage.n1 = 0;
  m_powerLawCalcStorage.sx1 = 0;
  m_powerLawCalcStorage.sxx1 = 0;
  m_powerLawCalcStorage.sy1 = 0;
  m_powerLawCalcStorage.sxy1 = 0;
  m_powerLawCalcStorage.N1 = 0;

  m_powerLawCalcStorage.n2 = 0;
  m_powerLawCalcStorage.sx2 = 0;
  m_powerLawCalcStorage.sxx2 = 0;
  m_powerLawCalcStorage.sy2 = 0;
  m_powerLawCalcStorage.sxy2 = 0;
  m_powerLawCalcStorage.sx2mc2 = 0;
  m_powerLawCalcStorage.N2 = 0;

  // Loop over all pixels to make necessary pre-calculations
  for (Int32 pixelIdx = 0; pixelIdx < lnCurve.size(); pixelIdx++) {
    Float64 xi = lnCurve.getLambdaAt(pixelIdx);
    Float64 yi = lnCurve.getFluxAt(pixelIdx);
    Float64 wi = 1;

    if (xi < m_powerLawCalcStorage.xc) {
      if (coefsFirstEstim.has_value()) {
        Float64 estimatedFlux =
            computePowerLaw(coefsFirstEstim.value().first, std::exp(xi));
        wi = (estimatedFlux * estimatedFlux) /
             (lnCurve.getFluxErrorAt(pixelIdx) *
              lnCurve.getFluxErrorAt(pixelIdx));
      }
      m_powerLawCalcStorage.n1 += wi;
      m_powerLawCalcStorage.sx1 += xi * wi;
      m_powerLawCalcStorage.sxx1 += xi * xi * wi;
      m_powerLawCalcStorage.sy1 += yi * wi;
      m_powerLawCalcStorage.sxy1 += xi * yi * wi;
      m_powerLawCalcStorage.N1 += 1;
    } else {
      if (coefsFirstEstim.has_value()) {
        Float64 estimatedFlux =
            computePowerLaw(coefsFirstEstim.value().second, std::exp(xi));
        wi = (estimatedFlux * estimatedFlux) /
             (lnCurve.getFluxErrorAt(pixelIdx) *
              lnCurve.getFluxErrorAt(pixelIdx));
      }
      m_powerLawCalcStorage.n2 += wi;
      m_powerLawCalcStorage.sx2 += xi * wi;
      m_powerLawCalcStorage.sxx2 += xi * xi * wi;
      m_powerLawCalcStorage.sy2 += yi * wi;
      m_powerLawCalcStorage.sxy2 += xi * yi * wi;
      m_powerLawCalcStorage.sx2mc2 += wi * (m_powerLawCalcStorage.xc - xi) *
                                      (m_powerLawCalcStorage.xc - xi);
      m_powerLawCalcStorage.N2 += 1;
    }
  }

  // Creates the Mc^TN^-1Mc matrix
  m_powerLawCalcStorage.m(0, 0) =
      m_powerLawCalcStorage.n1 + m_powerLawCalcStorage.n2;
  m_powerLawCalcStorage.m(0, 1) =
      m_powerLawCalcStorage.n2 * m_powerLawCalcStorage.xc +
      m_powerLawCalcStorage.sx1;
  m_powerLawCalcStorage.m(0, 2) =
      -m_powerLawCalcStorage.n2 * m_powerLawCalcStorage.xc +
      m_powerLawCalcStorage.sx2;
  m_powerLawCalcStorage.m(1, 0) =
      m_powerLawCalcStorage.n2 * m_powerLawCalcStorage.xc +
      m_powerLawCalcStorage.sx1;
  m_powerLawCalcStorage.m(1, 1) = m_powerLawCalcStorage.n2 *
                                      m_powerLawCalcStorage.xc *
                                      m_powerLawCalcStorage.xc +
                                  m_powerLawCalcStorage.sxx1;
  m_powerLawCalcStorage.m(1, 2) =
      -m_powerLawCalcStorage.n2 * m_powerLawCalcStorage.xc *
          m_powerLawCalcStorage.xc +
      m_powerLawCalcStorage.xc * m_powerLawCalcStorage.sx2;
  m_powerLawCalcStorage.m(2, 0) =
      -m_powerLawCalcStorage.n2 * m_powerLawCalcStorage.xc +
      m_powerLawCalcStorage.sx2;
  m_powerLawCalcStorage.m(2, 1) =
      -m_powerLawCalcStorage.n2 * m_powerLawCalcStorage.xc *
          m_powerLawCalcStorage.xc +
      m_powerLawCalcStorage.xc * m_powerLawCalcStorage.sx2;
  m_powerLawCalcStorage.m(2, 2) = m_powerLawCalcStorage.sx2mc2;

  m_powerLawCalcStorage.v(0) =
      m_powerLawCalcStorage.sy1 + m_powerLawCalcStorage.sy2;
  m_powerLawCalcStorage.v(1) =
      m_powerLawCalcStorage.sxy1 +
      m_powerLawCalcStorage.xc * m_powerLawCalcStorage.sy2;
  m_powerLawCalcStorage.v(2) =
      m_powerLawCalcStorage.sxy2 -
      m_powerLawCalcStorage.xc * m_powerLawCalcStorage.sy2;
}

TPowerLawCoefsPair COperatorPowerLaw::computeDoublePowerLawCoefs(
    TCurve const &lnCurve,
    std::optional<TPowerLawCoefsPair> const &coefsFirstEstim) {
  // NB : coefsFirstEstim are used to calculate the weights. If absent
  // calculations are made without taking error into account

  updatePowerLawCalcStorage(lnCurve, coefsFirstEstim);

  if (std::abs(m_powerLawCalcStorage.m.determinant()) < DBL_MIN)
    THROWG(ErrorCode::INTERNAL_ERROR,
           "Cannot calculate power law coefs: division by zero");

  Eigen::Matrix3d mInv = m_powerLawCalcStorage.m.inverse();

  Eigen::Vector3d theta = mInv * m_powerLawCalcStorage.v;

  Float64 a1 = std::exp(theta(0));
  Float64 b1 = theta(1);
  Float64 b2 = theta(2);
  auto [a2, sigmaa2] = computea2(a1, b1, b2, mInv(0, 0), mInv(1, 1), mInv(2, 2),
                                 mInv(1, 2), mInv(0, 1), mInv(0, 2));

  // Calculates var / covar

  Float64 sigmaa1 = stdExpA(a1, mInv(0, 0));
  Float64 sigmab1 = std::sqrt(mInv(1, 1));
  Float64 sigmab2 = std::sqrt(mInv(2, 2));

  return {{a1, b1, sigmaa1, sigmab1}, {a2, b2, sigmaa2, sigmab2}};
}

TPowerLawCoefsPair
COperatorPowerLaw::computeDoublePowerLawCoefs_b2_fixed(Float64 b2) {

  Eigen::Matrix2d m;
  m(0, 0) = m_powerLawCalcStorage.m(0, 0);
  m(0, 1) = m_powerLawCalcStorage.m(0, 1);
  m(1, 0) = m_powerLawCalcStorage.m(1, 0);
  m(1, 1) = m_powerLawCalcStorage.m(1, 1);

  Eigen::Vector2d v;
  auto g = m_powerLawCalcStorage.n2 * b2 * m_powerLawCalcStorage.xc -
           b2 * m_powerLawCalcStorage.sx2;
  v(0) = m_powerLawCalcStorage.v(0) + g;
  v(1) = m_powerLawCalcStorage.v(1) + m_powerLawCalcStorage.xc * g;

  Eigen::Matrix2d mInv = m.inverse();

  Eigen::Vector2d theta = mInv * v;

  Float64 a1 = std::exp(theta(0));
  Float64 b1 = theta(1);

  const Float64 varA1 = mInv(0, 0);
  const Float64 varb1 = mInv(1, 1);

  Float64 sigmaa1 = stdExpA(a1, varA1);
  Float64 sigmab1 = std::sqrt(varb1);

  auto [a2, sigmaa2] = computea2(a1, b1, b2, varA1, varb1, 0, 0, mInv(0, 1), 0);

  // If once recomputed b2 is out of bounds, we recompute a1 and a2 with both b1
  // and b2 fixed to limits
  const auto newb1 = limitCoef(b1, m_powerCoefsLimits.first);
  if (newb1 != b1) {
    return computeDoublePowerLawCoefs_b1_b2_fixed(newb1, b2);
  }

  return {{a1, b1, sigmaa1, sigmab1}, {a2, b2, sigmaa2, 0}};
}

TPowerLawCoefsPair
COperatorPowerLaw::computeDoublePowerLawCoefs_b1_fixed(Float64 b1) {

  Eigen::Matrix2d m;
  m(0, 0) = m_powerLawCalcStorage.m(0, 0);
  m(0, 1) = m_powerLawCalcStorage.n1 * m_powerLawCalcStorage.xc +
            m_powerLawCalcStorage.sx2;
  m(1, 0) = m(0, 1);
  m(1, 1) = m_powerLawCalcStorage.n1 * m_powerLawCalcStorage.xc *
                m_powerLawCalcStorage.xc +
            m_powerLawCalcStorage.sxx2;

  Eigen::Vector2d v;
  auto g = m_powerLawCalcStorage.n1 * b1 * m_powerLawCalcStorage.xc -
           b1 * m_powerLawCalcStorage.sx1;
  v(0) = m_powerLawCalcStorage.v(0) + g;
  v(1) = m_powerLawCalcStorage.xc * m_powerLawCalcStorage.sy1 +
         m_powerLawCalcStorage.sxy2 + m_powerLawCalcStorage.xc * g;

  Eigen::Matrix2d mInv = m.inverse();

  Eigen::Vector2d theta = mInv * v;

  Float64 a2 = std::exp(theta(0));

  const Float64 varA2 = mInv(0, 0);
  const Float64 varb2 = mInv(1, 1);
  Float64 b2 = theta(1);
  auto [a1, sigmaa1] = computea1(a2, b1, b2, varA2, 0, varb2, 0, 0, mInv(0, 1));

  Float64 sigmaa2 = stdExpA(a2, varA2);
  Float64 sigmab2 = std::sqrt(varb2);

  const auto newb2 = limitCoef(b2, m_powerCoefsLimits.second);
  // If once recomputed b2 is out of bounds, we recompute a1 and a2 with both b1
  // and b2 fixed to limits
  if (newb2 != b2) {
    return computeDoublePowerLawCoefs_b1_b2_fixed(b1, newb2);
  }

  return {{a1, b1, sigmaa1, 0}, {a2, b2, sigmaa2, sigmab2}};
}

TPowerLawCoefsPair
COperatorPowerLaw::computeDoublePowerLawCoefs_b1_b2_fixed(Float64 b1,
                                                          Float64 b2) {

  auto A1 = 1 / (m_powerLawCalcStorage.n1 + m_powerLawCalcStorage.n2) *
            (m_powerLawCalcStorage.sy1 + m_powerLawCalcStorage.sy2 -
             b1 * m_powerLawCalcStorage.sx1 - b2 * m_powerLawCalcStorage.sx2 +
             m_powerLawCalcStorage.n2 * (b2 - b1) * m_powerLawCalcStorage.xc);
  Float64 a1 = std::exp(A1);
  const Float64 varA1 =
      1 / (m_powerLawCalcStorage.n1 + m_powerLawCalcStorage.n2);
  auto [a2, sigmaa2] = computea2(a1, b1, b2, varA1, 0, 0, 0, 0, 0);
  Float64 sigmaa1 = stdExpA(a1, varA1);

  return {{a1, b1, sigmaa1, 0}, {a2, b2, sigmaa2, 0}};
}

TCurve COperatorPowerLaw::initializeFluxCurve(Float64 redshift,
                                              Float64 nullFluxThreshold) {
  // In order to take into account ism, igm we initialize a template with a flux
  // at 1 we then divide the initial flux by the template flux value once ism
  // igm has been applied

  TList<Float64> spectrumLambda;
  TList<Float64> spectrumFlux;
  TList<Float64> spectrumFluxError;
  // Concatenates all curves
  // NB at the end, lambda is not ordered anymore
  for (Int32 spectrumIdx = 0; spectrumIdx < m_nSpectra; spectrumIdx++) {
    TList<Float64> tmpLambda =
        m_spectra[spectrumIdx]
            ->GetSpectralAxis()
            .extract(m_kStart[spectrumIdx], m_kEnd[spectrumIdx])
            .GetSamplesVector();
    spectrumLambda.insert(spectrumLambda.end(),
                          std::make_move_iterator(tmpLambda.begin()),
                          std::make_move_iterator(tmpLambda.end()));

    TList<Float64> tmpFlux =
        m_spectra[spectrumIdx]
            ->GetFluxAxis()
            .extract(m_kStart[spectrumIdx], m_kEnd[spectrumIdx])
            .GetSamplesVector();
    spectrumFlux.insert(spectrumFlux.end(),
                        std::make_move_iterator(tmpFlux.begin()),
                        std::make_move_iterator(tmpFlux.end()));

    TList<Float64> tmpError =
        m_spectra[spectrumIdx]
            ->GetFluxAxis()
            .GetError()
            .extract(m_kStart[spectrumIdx], m_kEnd[spectrumIdx])
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
  TCurve fluxCurve1D(spectrumLambdaAxis.GetSamplesVector(),
                     std::move(spectrumFlux), std::move(spectrumFluxError),
                     std::move(maskedPixels),
                     TBoolList(spectrumLambdaAxis.GetSamplesCount(), false),
                     std::move(snrCompliantPixels));
  fluxCurve1D.sort();

  return fluxCurve1D;
}

T3DList<Float64> COperatorPowerLaw::computeIsmIgmCorrections(
    Float64 redshift, CSpectrumSpectralAxis const &spectrumLambdaRest,
    bool opt_extinction, bool opt_dustFitting) const {
  // In order to access ism igm coefs, we initialize a template with a flux
  // at 1, and apply ism/igm on it

  T3DList<Float64> correctionCoefs(
      m_nIgmCurves,
      std::vector<std::vector<Float64>>(
          m_nIsmCurves,
          std::vector<Float64>(spectrumLambdaRest.GetSamplesCount(), NAN)));
  CTemplate templateForCoefs(
      "", "", spectrumLambdaRest,
      std::vector<Float64>(spectrumLambdaRest.GetSamplesCount(), 1));
  templateForCoefs.InitIsmIgmConfig(redshift);
  for (Int32 igmIdx = 0; igmIdx < m_nIgmCurves; igmIdx++) {
    if (opt_extinction) { // igm
      templateForCoefs.ApplyMeiksinCoeff(m_igmIdxList[igmIdx]);
    }
    for (Int32 ismIdx = 0; ismIdx < m_nIsmCurves; ismIdx++) {
      if (opt_dustFitting) { // ism
        templateForCoefs.ApplyDustCoeff(m_ismIdxList[ismIdx]);
      }
      correctionCoefs[igmIdx][ismIdx] =
          templateForCoefs.GetFluxAxis().GetSamplesVector();
    }
  }
  return correctionCoefs;
}

TList<Float64> COperatorPowerLaw::computeIsmIgmCorrection(
    Float64 redshift, CSpectrumSpectralAxis const &spectrumLambdaRest,
    Int32 igmIdx, Float64 ismCoef) const {
  // In order to access ism igm coefs, we initialize a template with a flux
  // at 1, and apply ism/igm on it
  CTemplate templateForCoefs(
      "", "", spectrumLambdaRest,
      std::vector<Float64>(spectrumLambdaRest.GetSamplesCount(), 1));
  templateForCoefs.InitIsmIgmConfig(redshift);

  TList<Float64> correctionCoefs(spectrumLambdaRest.GetSamplesCount(), NAN);

  if (igmIdx > -1) {
    templateForCoefs.ApplyMeiksinCoeff(igmIdx);
  }
  if (ismCoef > 0) {
    Int32 ismIdx = -1;
    ismIdx = m_ismCorrectionCalzetti->GetEbmvIndex(ismCoef);
    templateForCoefs.ApplyDustCoeff(ismIdx);
  }
  correctionCoefs = templateForCoefs.GetFluxAxis().GetSamplesVector();
  return correctionCoefs;
}

T3DCurve COperatorPowerLaw::computeEmittedCurve(Float64 redshift,
                                                bool opt_extinction,
                                                bool opt_dustFitting,
                                                TCurve &&fluxCurve_) {
  // In order to take into account ism, igm we initialize a template with a flux
  // at 1 we then divide the initial flux by the template flux value once ism
  // igm has been applied
  // Set m_ismIgmCorrections

  T3DCurve fluxCurve(std::move(fluxCurve_), m_nIgmCurves, m_nIsmCurves);

  T3DList<bool> isExtincted(
      m_nIgmCurves,
      T2DList<bool>(m_nIsmCurves, TList<bool>(fluxCurve.size(), false)));
  if (opt_extinction || opt_dustFitting) {
    auto const ismIgmCorrections = computeIsmIgmCorrections(
        redshift, fluxCurve.getLambda(), opt_extinction, opt_dustFitting);
    for (Int32 igmIdx = 0; igmIdx < m_nIgmCurves; igmIdx++) {
      for (Int32 ismIdx = 0; ismIdx < m_nIsmCurves; ismIdx++) {
        for (Int32 pixelIdx = 0; pixelIdx < fluxCurve.size(); pixelIdx++) {
          Float64 correctionCoef = ismIgmCorrections[igmIdx][ismIdx][pixelIdx];
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
  for (Int32 pixelIdx = 0; pixelIdx < ssize(spectrumFlux); pixelIdx++) {
    // Checks that signal / noise ratio is sufficient (> nullFluxThreshold)
    if (spectrumFlux[pixelIdx] <
        nullFluxThreshold * spectrumFluxError[pixelIdx]) {
      snrCompliant[pixelIdx] = false;
    }
  }
  return snrCompliant;
}

TFloat64List COperatorPowerLaw::computeModelFlux(
    const CSpectrumSpectralAxis &lambdaRestAxis, const Float64 redshift,
    const Int32 meiksinIdx, const Float64 ebmvCoef,
    const TPowerLawCoefsPair &coefs) const {
  TList<Float64> correctionCoefs(lambdaRestAxis.GetSamplesCount(), 1.0);
  if (meiksinIdx || ebmvCoef)
    correctionCoefs =
        computeIsmIgmCorrection(redshift, lambdaRestAxis, meiksinIdx, ebmvCoef);
  TList<Float64> fluxObs(lambdaRestAxis.GetSamplesCount(), NAN);
  for (Int32 pixelIdx = 0; pixelIdx < lambdaRestAxis.GetSamplesCount();
       pixelIdx++) {
    fluxObs[pixelIdx] = computeDoublePowerLaw(coefs, lambdaRestAxis[pixelIdx]) *
                        correctionCoefs[pixelIdx];
  }
  return fluxObs;
}

CModelSpectrumResult COperatorPowerLaw::ComputeSpectrumModel(
    const CContinuumModelSolution &continuum, Int32 spcIndex) {

  auto const &lambdaObsAxis = m_spectra[spcIndex]->GetSpectralAxis();
  auto const &lambdaObs = lambdaObsAxis.GetSamplesVector();
  auto const lambdaRestAxis = lambdaObsAxis.blueShift(continuum.redshift);

  auto fluxObs = computeModelFlux(
      lambdaRestAxis, continuum.redshift, continuum.meiksinIdx,
      continuum.ebmvCoef,
      {{continuum.a1, continuum.b1}, {continuum.a2, continuum.b2}});

  return CModelSpectrumResult(lambdaObs, std::move(fluxObs),
                              m_spectra[spcIndex]->getObsID());
}

Float64 COperatorPowerLaw::stdExpA(Float64 a, Float64 varA) const {
  // From var(a) = var(exp(A)) ~ (exp(A))^2 * var(A) => std(a) = a *
  // sqrt(var(A))
  return a * std::sqrt(varA);
}

std::pair<Float64, Float64>
COperatorPowerLaw::computea2(const Float64 a1, const Float64 b1,
                             const Float64 b2, Float64 varA1, Float64 varb1,
                             Float64 varb2, Float64 covb1b2, Float64 covA1b1,
                             Float64 covA1b2) const {
  Float64 a2 = a1 * std::pow(m_lambdaCut, b1 - b2);
  Float64 varA2 = computeVarA2(varA1, varb1, varb2, covb1b2, covA1b1, covA1b2);
  Float64 stdaa2 = stdExpA(a2, varA2);
  return std::make_pair(a2, stdaa2);
}

std::pair<Float64, Float64>
COperatorPowerLaw::computea1(const Float64 a2, const Float64 b1,
                             const Float64 b2, Float64 varA2, Float64 varb1,
                             Float64 varb2, Float64 covb1b2, Float64 covA2b1,
                             Float64 covA2b2) const {
  Float64 a1 = a2 * std::pow(m_lambdaCut, b2 - b1);
  Float64 varA1 = computeVarA1(varA2, varb1, varb2, covb1b2, covA2b1, covA2b2);
  Float64 stdaa1 = stdExpA(a1, varA1);
  return std::make_pair(a1, stdaa1);
}

Float64 COperatorPowerLaw::computeVarA2(Float64 varA1, Float64 varb1,
                                        Float64 varb2, Float64 covb1b2,
                                        Float64 covA1b1,
                                        Float64 covA1b2) const {
  // Computes var(A2) from var and covar of A1, b1, b2. Can be used for var(a1)
  // exchanging all 2 with 1
  return varA1 +
         m_powerLawCalcStorage.xc * m_powerLawCalcStorage.xc *
             (varb1 + varb2 - 2 * covb1b2) +
         2 * m_powerLawCalcStorage.xc * (covA1b1 - covA1b2);
}

Float64 COperatorPowerLaw::computeVarA1(Float64 varA2, Float64 varb1,
                                        Float64 varb2, Float64 covb1b2,
                                        Float64 covA2b1,
                                        Float64 covA2b2) const {
  // Computes var(A1) from var and covar of A2, b1, b2.
  return computeVarA2(varA2, varb2, varb1, covb1b2, covA2b2, covA2b1);
}