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
#include <climits>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <boost/range/combine.hpp>

#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/operator/templatefittinglog.h"
#include "RedshiftLibrary/operator/templatefittingresult.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/spectrum/fullspectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/statistics/fitquality.h"

using namespace NSEpic;
using namespace std;

COperatorTemplateFittingLog::COperatorTemplateFittingLog(
    const TFloat64List &redshifts)
    : COperatorTemplateFittingBase(redshifts) {

  CheckRedshifts();
}

void COperatorTemplateFittingLog::SetRedshifts(const TFloat64List &redshifts) {
  COperatorTemplateFittingBase::SetRedshifts(redshifts);
  CheckRedshifts();
}

void COperatorTemplateFittingLog::CheckRedshifts() {
  if (m_redshifts.size() < 2) {
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Invalid redshift array: " << m_redshifts.size()
                       << " <2");
  }

  // check if spectrum sampling is a multiple of redshift step
  m_logstep = log((m_redshifts[1] + 1) / (m_redshifts[0] + 1));

  Float64 modulo;

  m_ssRatio = Context.getRebinnedFullSpectra()[0]
                  ->GetSpectralAxis()
                  .GetLogSamplingIntegerRatio(m_logstep, modulo);
  if (std::abs(modulo) > 1E-12)
    THROWG(ErrorCode::INTERNAL_ERROR, "spc and tpl do not have a lambdastep "
                                      "multiple of redshift step");
  for (auto it = Context.getRebinnedFullSpectra().cbegin() + 1,
            end = Context.getRebinnedFullSpectra().cend();
       it != end; ++it) {
    if ((*it)->GetSpectralAxis().GetLogSamplingIntegerRatio(
            m_logstep, modulo) != m_ssRatio)
      THROWG(ErrorCode::INTERNAL_ERROR,
             "spc multiobs do not share same lambda step");
    if (std::abs(modulo) > 1E-12)
      THROWG(ErrorCode::INTERNAL_ERROR, "spc and tpl do not have a lambdastep "
                                        "multiple of redshift step");
  }

  if (m_ssRatio == 1) {
    m_spectraFull.clear();
    for (const auto &spectrum : Context.getRebinnedFullSpectra())
      m_spectraFull.push_back(spectrum);
    m_lambdaRanges = Context.getRebinnedFullClampedLambdaRanges();
    return;
  }

  // else subsampling required, subsample each spectrum :
  //  (coarse redshift grid)
  m_spectraFull.clear();
  m_lambdaRanges.clear();
  m_spectraFull.reserve(Context.getSpectra().size());
  m_lambdaRanges.reserve(Context.getSpectra().size());

  for (auto const &[logSampledSpectrum_ptr, logSampledLambdaRange_ptr] :
       boost::combine(Context.getRebinnedFullSpectra(),
                      Context.getRebinnedFullClampedLambdaRanges())) {

    TMaskList mask_spc =
        logSampledSpectrum_ptr->GetSpectralAxis().GetSubSamplingMask(
            m_ssRatio, *logSampledLambdaRange_ptr);
    std::shared_ptr<CFullSpectrum> ssSpectrum =
        std::make_shared<CFullSpectrum>(*logSampledSpectrum_ptr, mask_spc);
    // scale the variance by ssratio
    CSpectrumNoiseAxis scaledNoise = ssSpectrum->GetErrorAxis();
    scaledNoise *= 1. / sqrt(m_ssRatio);
    ssSpectrum->SetErrorAxis(std::move(scaledNoise));

    // double make sure that subsampled spectrum is well sampled
    if (!ssSpectrum->GetSpectralAxis().IsLogSampled(m_logstep)) {
      THROWG(ErrorCode::INTERNAL_ERROR,
             "subsampled spectrum is not correctly log sampled "
             "at the redshift step");
    }

    // set the lambda range to the clamped subsampled spectrum
    std::shared_ptr<TFloat64Range> ssLambdaRange =
        std::make_shared<TFloat64Range>();
    ssSpectrum->GetSpectralAxis().ClampLambdaRange(*logSampledLambdaRange_ptr,
                                                   *ssLambdaRange);

    m_spectraFull.push_back(std::move(ssSpectrum));
    m_lambdaRanges.push_back(std::move(ssLambdaRange));
  }
}

// only works for mtm, Y=model^2, X=1.
Int32 COperatorTemplateFittingLog::EstimateMtMFast(const TFloat64List &X,
                                                   const TFloat64List &Y,
                                                   Int32 nShifts,
                                                   TFloat64List &XtY) {
  XtY.resize(nShifts);

  Int32 nX = X.size();
  Float64 xty = 0.0;
  for (Int32 j = 0; j < nX; j++) {
    xty += X[j] * Y[j];
  }
  XtY[0] = xty;

  for (Int32 k = 1; k < nShifts; k++) {
    xty = XtY[k - 1];
    xty -= X[0] * Y[0 + k - 1];
    xty += X[nX - 1] * Y[nX - 1 + k - 1];

    XtY[k] = xty;
  }
  return 0;
}

void COperatorTemplateFittingLog::EstimateXtY(
    const TFloat64List &X, const TFloat64List &Y, TFloat64List &XtY,
    FFTPlans &plans, EPrecomputedFFT fftX, EPrecomputedFFT fftY) {

  // Processing the FFT
  Int32 nX = X.size();
  Int32 nY = Y.size();
  Int32 nshifts = nY - nX + 1;
  Int32 nPadded = plans.nPaddedSamples;

  Int32 nPadBeforeSpc = nPadded - nX;

  Log.LogDebug(Formatter() << __func__
                           << ": Processing X-fft "
                              "with n="
                           << nX << ", padded to n=" << nPadded);

  bool computeXfft =
      (fftX == EPrecomputedFFT::none) || (!plans.isPrecomputed(fftX));
  if (computeXfft) {
    Log.LogDebug(Formatter()
                 << __func__
                 << ": Processing X-fft with nPadBeforeSpc=" << nPadBeforeSpc);
    for (Int32 k = 0; k < nPadBeforeSpc; k++)
      plans.inX[k] = 0.0;
    for (Int32 k = nPadBeforeSpc; k < nPadBeforeSpc + nX; k++)
      plans.inX[k] = X[nX - 1 - (k - nPadBeforeSpc)];

    fftw_execute(plans.pX);
    Log.LogDebug(Formatter() << __func__ << ": X-fft done");

    if (fftX != EPrecomputedFFT::none)
      plans.storeFFT(fftX, plans.outX);
  } else {
    plans.getFFT(fftX, plans.outX);
  }

  Log.LogDebug(Formatter() << __func__ << ": Processing X-fft with n=" << nY
                           << ", padded to n=" << nPadded);

  bool computeYfft =
      (fftY == EPrecomputedFFT::none) || (!plans.isPrecomputed(fftY));
  if (computeYfft) {
    Log.LogDebug(Formatter()
                 << __func__
                 << ": Processing Y-fft with nPadBeforeSpc=" << nPadBeforeSpc);
    for (Int32 k = 0; k < nY; k++)
      plans.inY[k] = Y[k];
    for (Int32 k = nY; k < nPadded; k++)
      plans.inY[k] = 0.0;

    fftw_execute(plans.pY);
    Log.LogDebug(Formatter() << __func__ << ": Y-fft done");

    if (fftY != EPrecomputedFFT::none)
      plans.storeFFT(fftY, plans.outY);
  } else {
    plans.getFFT(fftY, plans.outY);
  }

  // Multiplying the FFT outputs
  for (Int32 k = 0; k < nPadded; k++) {
    plans.outXY[k][0] = (plans.outY[k][0] * plans.outX[k][0] -
                         plans.outY[k][1] * plans.outX[k][1]);
    plans.outXY[k][1] = (plans.outY[k][0] * plans.outX[k][1] +
                         plans.outY[k][1] * plans.outX[k][0]);
  }

  fftw_execute(plans.pBackward);
  Log.LogDebug(Formatter() << __func__ << ": backward-fft done");

  XtY.resize(nshifts);
  XtY.front() = plans.inXY[nPadded - 1] / Float64(nPadded);
  for (Int32 k = 1; k < nshifts; k++)
    XtY[k] = plans.inXY[k - 1] / Float64(nPadded);
}

FFTPlans::FFTPlans(FFTPlans &&other) { *this = std::move(other); }

FFTPlans &FFTPlans::operator=(FFTPlans &&other) {
  freeFFTPlans();
  std::swap(nPaddedSamples, other.nPaddedSamples);
  std::swap(inX, other.inX);
  std::swap(outX, other.outX);
  std::swap(pX, other.pX);
  std::swap(inY, other.inY);
  std::swap(outY, other.outY);
  std::swap(pY, other.pY);
  std::swap(outXY, other.outXY);
  std::swap(inXY, other.inXY);
  std::swap(pBackward, other.pBackward);
  for (auto &[key, other_buffer] : other.precomputedFFT) {
    if (precomputedFFT.find(key) == precomputedFFT.end()) {
      precomputedFFT[key] = other_buffer;
      other.precomputedFFT.at(key) = nullptr;
    } else {
      std::swap(precomputedFFT.at(key), other.precomputedFFT.at(key));
    }
  }
  return *this;
}

void FFTPlans::initFFT(Int32 nPadded) {

  nPaddedSamples = nPadded;
  inX = (Float64 *)fftw_malloc(sizeof(Float64) * nPadded);
  outX = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPadded);
  pX = fftw_plan_dft_r2c_1d(nPadded, inX, outX, FFTW_ESTIMATE);
  if (inX == 0) {
    THROWG(ErrorCode::INTERNAL_ERROR, "Unable to allocate inSpc");
  }
  if (outX == 0) {
    THROWG(ErrorCode::INTERNAL_ERROR, "Unable to allocate outSpc");
  }

  inY = (Float64 *)fftw_malloc(sizeof(Float64) * nPadded);
  outY = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPadded);
  pY = fftw_plan_dft_r2c_1d(nPadded, inY, outY, FFTW_ESTIMATE);
  if (inY == 0) {
    THROWG(ErrorCode::INTERNAL_ERROR, "Unable to allocate inTpl");
  }
  if (outY == 0) {
    THROWG(ErrorCode::INTERNAL_ERROR, "Unable to allocate outTpl");
  }

  outXY = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPadded);
  inXY = (Float64 *)fftw_malloc(sizeof(Float64) * nPadded);
  pBackward = fftw_plan_dft_c2r_1d(nPadded, outXY, inXY, FFTW_ESTIMATE);
  if (outXY == 0) {
    THROWG(ErrorCode::INTERNAL_ERROR, "Unable to "
                                      "allocate outCombined");
  }
  if (inXY == 0) {
    THROWG(ErrorCode::INTERNAL_ERROR, "Unable to "
                                      "allocate inCombined");
  }
}

void FFTPlans::freeFFTPrecomputedBuffers() {
  for (auto &[_, fft] : precomputedFFT)
    if (fft != nullptr) {
      fftw_free(fft);
      fft = nullptr;
    }
  precomputedFFT.clear();
}

void FFTPlans::freeFFTPlans() {
  if (pX) {
    fftw_destroy_plan(pX);
    pX = nullptr;
  }
  if (inX) {
    fftw_free(inX);
    inX = nullptr;
  }
  if (outX) {
    fftw_free(outX);
    outX = nullptr;
  }
  if (pY) {
    fftw_destroy_plan(pY);
    pY = nullptr;
  }
  if (inY) {
    fftw_free(inY);
    inY = nullptr;
  }
  if (outY) {
    fftw_free(outY);
    outY = nullptr;
  }
  if (pBackward) {
    fftw_destroy_plan(pBackward);
    pBackward = nullptr;
  }
  if (inXY) {
    fftw_free(inXY);
    inXY = nullptr;
  }
  if (outXY) {
    fftw_free(outXY);
    outXY = nullptr;
  }

  freeFFTPrecomputedBuffers();
}

void FFTPlans::storeFFT(EPrecomputedFFT precomputed, fftw_complex *fft) {
  if (!isPrecomputed(precomputed))
    allocatePrecomputedFFT(precomputed);
  copyFFT(fft, precomputedFFT.at(precomputed));
}

void FFTPlans::getFFT(EPrecomputedFFT precomputed, fftw_complex *fft) {
  copyFFT(precomputedFFT.at(precomputed), fft);
}

void FFTPlans::copyFFT(fftw_complex *src, fftw_complex *dest) {
  for (Int32 k = 0; k < nPaddedSamples; k++) {
    dest[k][0] = src[k][0];
    dest[k][1] = src[k][1];
  }
}

bool FFTPlans::isPrecomputed(EPrecomputedFFT precomputed) {
  auto const it = precomputedFFT.find(precomputed);
  if (it == precomputedFFT.end())
    return false;
  if (it->second == nullptr)
    return false;
  return true;
}

void FFTPlans::allocatePrecomputedFFT(EPrecomputedFFT precomputed) {
  if (isPrecomputed(precomputed))
    THROWG(ErrorCode::INTERNAL_ERROR,
           "Not allocating precomputedFFT since alreay allocated");
  auto buffer =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPaddedSamples);
  if (buffer == nullptr)
    THROWG(ErrorCode::INTERNAL_ERROR, "Unable to allocate precomputedFFT");
  precomputedFFT[precomputed] = buffer;
}

TInt32RangeList
COperatorTemplateFittingLog::FindZRanges(const TFloat64List &redshifts) {
  TInt32RangeList izrangelist;
  TInt32List zsplit;
  if (m_enableIGM) {
    Float64 zmin_igm =
        GetIGMStartingRedshiftValue(m_spectraFull[0]->GetSpectralAxis()[0]);
    if (zmin_igm > redshifts.front() && zmin_igm < redshifts.back()) {
      Int32 i_zmin_igm = -1;
      TFloat64Index::getClosestLowerIndex(redshifts, zmin_igm, i_zmin_igm);
      zsplit.push_back(i_zmin_igm);
    }
    if (zmin_igm < redshifts.back()) {
      TFloat64List zbins_igm =
          m_templateRebined_bf[0].m_igmCorrectionMeiksin->getRedshiftBins();
      zbins_igm.pop_back(); // extend last bin up to the end
      for (const Float64 &z : zbins_igm) {
        if (z <= zmin_igm || z <= redshifts.front())
          continue;
        if (z >= redshifts.back())
          break;
        Int32 iz = -1;
        TFloat64Index::getClosestLowerIndex(redshifts, z, iz);
        zsplit.push_back(iz);
      }
    }

    Log.LogDebug(Formatter() << "FitAllz: indexes for full "
                                "LstSquare calculation, count = "
                             << zsplit.size());
    for (Int32 k = 0; k < ssize(zsplit); k++) {
      Log.LogDebug(Formatter() << "FitAllz: indexes ranges: "
                                  "for i="
                               << k << ", zsplit=" << zsplit[k]);
    }
  }

  Int32 izmin = 0;
  for (const Int32 &izmax : zsplit) {
    izrangelist.push_back(TInt32Range(izmin, izmax));
    izmin = izmax + 1;
  }
  Int32 izmax = redshifts.size() - 1;
  izrangelist.push_back(TInt32Range(izmin, izmax));

  Log.LogDebug(Formatter() << "FitAllz: indexes - "
                              "izrangelist calculation, count = "
                           << izrangelist.size());
  for (Int32 k = 0; k < ssize(izrangelist); k++) {
    Log.LogDebug(Formatter()
                 << "FitAllz: indexes ranges: "
                    "for i="
                 << k << ", zmin=" << redshifts[izrangelist[k].GetBegin()]
                 << ", zmax=" << redshifts[izrangelist[k].GetEnd()]);
  }

  return izrangelist;
}

/**
 * @brief COperatorTemplateFittingLog::FitAllz
 *
 * @param lambdaRange
 * @param result
 * @param opt_extinction
 * @param opt_dustFitting
 * @param spcMaskAdditional
 * @param logpriorze: if size=0, prior is deactivated
 * @return
 */
void COperatorTemplateFittingLog::FitAllz(
    std::shared_ptr<CTemplateFittingResult> result,
    const TInt32List &MeiksinList, const TInt32List &EbmvList,
    const CPriorHelper::TPriorZEList &logpriorze, CMask const &lineMask) {

  // prepare list of redshifts segments that keep the same IGM extinction curves
  TInt32RangeList izrangelist = FindZRanges(result->Redshifts);
  Int32 nzranges = izrangelist.size();

  // since dtd is cte, better compute it here
  const TAxisSampleList &error =
      m_spectraFull[0]->GetFluxAxis().GetError().GetSamplesVector();

  const Int32 nRedshifts = result->Redshifts.size();
  const TAxisSampleList &spectrumRebinedFluxRaw =
      m_spectraFull[0]->GetFluxAxis().GetSamplesVector();
  const Int32 nSpcPixels = spectrumRebinedFluxRaw.size();
  auto const &mask = m_spectraFull[0]->getMask();
  Float64 dtd = 0.0;
  TFloat64List inv_err2(nSpcPixels);
  TFloat64List inv_err(nSpcPixels);
  TFloat64List spcRebinedFluxOverErr2(nSpcPixels);
  TFloat64List spcRebinedFlux2OverErr2(nSpcPixels);
  for (Int32 j = 0; j < ssize(error); j++) {
    if (mask[j]) {
      inv_err[j] = 1.0 / error[j];
      inv_err2[j] = inv_err[j] * inv_err[j];
      spcRebinedFluxOverErr2[j] = spectrumRebinedFluxRaw[j] * inv_err2[j];
      spcRebinedFlux2OverErr2[j] =
          spcRebinedFluxOverErr2[j] * spectrumRebinedFluxRaw[j];
      dtd += spcRebinedFlux2OverErr2[j];
    }
  }

  for (Int32 k = 0; k < nzranges; k++) {
    // prepare the zrange-result container
    std::shared_ptr<CTemplateFittingResult> subresult;
    TFloat64Range zrange =
        TFloat64Range(result->Redshifts[izrangelist[k].GetBegin()],
                      result->Redshifts[izrangelist[k].GetEnd()]);
    if (m_enableIGM && nRedshifts > 1) {
      TFloat64List::const_iterator first = result->Redshifts.begin() +
                                           izrangelist[k].GetBegin(),
                                   last = result->Redshifts.begin() +
                                          izrangelist[k].GetEnd() + 1;
      TFloat64List subRedshifts(first, last);
      subresult = std::make_shared<CTemplateFittingResult>(
          subRedshifts.size(), EbmvList.size(), MeiksinList.size());

      subresult->Redshifts = subRedshifts;

    } else {
      subresult = std::make_shared<CTemplateFittingResult>(
          nRedshifts, EbmvList.size(), MeiksinList.size());
      subresult->Redshifts = result->Redshifts;
    }

    TInt32Range ilbda = FindTplSpectralIndex(zrange);
    Log.LogDebug(Formatter() << "FitAllz: zrange min=" << zrange.GetBegin()
                             << ", max=" << zrange.GetEnd());
    Log.LogDebug(Formatter()
                 << "FitAllz: full zmin=" << result->Redshifts[0]
                 << ", full zmax=" << result->Redshifts[nRedshifts - 1]);
    Log.LogDebug(Formatter() << "FitAllz: indexes tpl crop: "
                                "lbda min="
                             << ilbda.GetBegin() << ", max=" << ilbda.GetEnd());
    Log.LogDebug(Formatter() << "FitAllz: indexes tpl full: "
                                "lbda min="
                             << 0 << ", max="
                             << m_templateRebined_bf[0].GetSampleCount() - 1);

    CMask lineMaskatThisRange =
        lineMask.GetMasksCount()
            ? lineMask.extract(ilbda.GetBegin(), ilbda.GetEnd())
            : CMask();

    FitRangez(inv_err2, spcRebinedFluxOverErr2, spcRebinedFlux2OverErr2, ilbda,
              subresult, MeiksinList, EbmvList, dtd, lineMaskatThisRange);

    // copy subresults into global results
    updateGlobalResults(result, subresult, izrangelist[k].GetBegin());

    if (logpriorze.size() > 0)
      applyPrior(result, subresult, izrangelist[k].GetBegin(), logpriorze, dtd);

    // computeFitQuality(result, izrangelist[k].GetBegin(),
    //                   subresult->Redshifts.size(), ilbda.GetBegin(),
    //                   lineMaskatThisRange);
  }
}

void COperatorTemplateFittingLog::updateGlobalResults(
    const std::shared_ptr<CTemplateFittingResult> &result,
    const std::shared_ptr<const CTemplateFittingResult> &subResult,
    Int32 resultIdx) {
  for (Int32 isubz = 0, fullResultIdx = resultIdx;
       isubz < ssize(subResult->Redshifts); ++isubz, ++fullResultIdx) {
    if (fullResultIdx >= ssize(result->ChiSquare))
      THROWG(ErrorCode::INTERNAL_ERROR, "out-of-bound index");
    result->ChiSquare[fullResultIdx] = subResult->ChiSquare[isubz];
    result->FitQuality[fullResultIdx] = subResult->FitQuality[isubz];
    result->FitAmplitude[fullResultIdx] = subResult->FitAmplitude[isubz];
    result->FitAmplitudeError[fullResultIdx] =
        subResult->FitAmplitudeError[isubz];
    result->FitAmplitudeSigma[fullResultIdx] =
        subResult->FitAmplitudeSigma[isubz];
    result->FitDtM[fullResultIdx] = subResult->FitDtM[isubz];
    result->FitMtM[fullResultIdx] = subResult->FitMtM[isubz];
    result->SNR[fullResultIdx] = subResult->SNR[isubz];
    result->Overlap[fullResultIdx] = subResult->Overlap[isubz];
    result->FitEbmvCoeff[fullResultIdx] = subResult->FitEbmvCoeff[isubz];
    result->FitMeiksinIdx[fullResultIdx] = subResult->FitMeiksinIdx[isubz];

    for (Int32 kigm = 0;
         kigm < ssize(result->IgmMeiksinIdxIntermediate[fullResultIdx]); kigm++)
      result->IgmMeiksinIdxIntermediate[fullResultIdx][kigm] =
          subResult->IgmMeiksinIdxIntermediate[isubz][kigm];
    for (Int32 kism = 0;
         kism < ssize(result->ChiSquareIntermediate[fullResultIdx]); kism++) {
      result->IsmEbmvIdxIntermediate[fullResultIdx][kism] =
          subResult->IsmEbmvIdxIntermediate[isubz][kism];
      for (Int32 kigm = 0;
           kigm < ssize(result->ChiSquareIntermediate[fullResultIdx][kism]);
           kigm++) {
        result->ChiSquareIntermediate[fullResultIdx][kism][kigm] =
            subResult->ChiSquareIntermediate[isubz][kism][kigm];
      }
    }
  }
}

void COperatorTemplateFittingLog::applyPrior(
    const std::shared_ptr<CTemplateFittingResult> &result,
    const std::shared_ptr<const CTemplateFittingResult> &subResult,
    Int32 resultIdx, const CPriorHelper::TPriorZEList &logpriorze,
    Float64 dtd) {
  for (Int32 isubz = 0, fullResultIdx = resultIdx;
       isubz < ssize(subResult->Redshifts); ++isubz, ++fullResultIdx) {
    if (fullResultIdx >= ssize(result->ChiSquare))
      THROWG(ErrorCode::INTERNAL_ERROR, "out-of-bound index");

    Float64 logprior = 0.;
    if (logpriorze.size() > 0) {
      Int32 kism_best = 0;
      if (subResult->FitEbmvCoeff[isubz] != -1.0)
        kism_best =
            m_templateRebined_bf[0].m_ismCorrectionCalzetti->GetEbmvIndex(
                subResult->FitEbmvCoeff[isubz]);

      const CPriorHelper::SPriorTZE &pTZE =
          logpriorze[fullResultIdx][kism_best];
      logprior += -2.0 * pTZE.betaTE * pTZE.logprior_precompTE;
      logprior += -2.0 * pTZE.betaA * pTZE.logprior_precompA;
      logprior += -2.0 * pTZE.betaZ * pTZE.logprior_precompZ;

      if (pTZE.A_sigma > 0.0 && pTZE.A_mean > 0.0) {
        // now update the amplitude if there is any constraints from the
        // priors
        Float64 ampl = result->FitAmplitude[fullResultIdx];
        Float64 ampl_err = result->FitAmplitudeError[fullResultIdx];
        Float64 ampl_sigma = result->FitAmplitudeSigma[fullResultIdx];
        if (pTZE.betaA > 0.0) {
          Float64 bss2 = pTZE.betaA / (pTZE.A_sigma * pTZE.A_sigma);
          ampl = (result->FitDtM[fullResultIdx] + pTZE.A_mean * bss2) /
                 (result->FitMtM[fullResultIdx] + bss2);
          ampl_err = sqrt(result->FitMtM[fullResultIdx]) /
                     (result->FitMtM[fullResultIdx] + bss2);

        } else {
          ampl = result->FitDtM[fullResultIdx] / result->FitMtM[fullResultIdx];
          ampl_err = sqrt(1. / result->FitMtM[fullResultIdx]);
        }

        Log.LogDebug(Formatter()
                     << "update the amplitude "
                        "(a_mean="
                     << pTZE.A_mean << ", a_sigma=" << pTZE.A_sigma);
        Log.LogDebug(Formatter() << "update the amplitude "
                                    "(ampl was = "
                                 << result->FitAmplitude[fullResultIdx]
                                 << ", updated to " << ampl);

        // check negative amplitude
        ampl_sigma = ampl / ampl_err;
        applyPositiveAndNonNullConstraint(ampl_sigma, ampl);

        result->FitAmplitude[fullResultIdx] = ampl;
        result->FitAmplitudeError[fullResultIdx] = ampl_err;
        result->FitAmplitudeSigma[fullResultIdx] = ampl_sigma;
        const Float64 chi2 = dtd + result->FitMtM[fullResultIdx] * ampl * ampl -
                             2. * ampl * result->FitDtM[fullResultIdx];
        result->ChiSquare[fullResultIdx] = chi2;
        Float64 logPa = pTZE.betaA * (ampl - pTZE.A_mean) *
                        (ampl - pTZE.A_mean) / (pTZE.A_sigma * pTZE.A_sigma);
        if (std::isnan(logPa) || logPa != logPa || std::isinf(logPa)) {
          THROWG(ErrorCode::INTERNAL_ERROR,
                 Formatter() << " Invalid logPa value (a_mean=" << pTZE.A_mean
                             << ", a_sigma=" << pTZE.A_sigma);
        }
        logprior += logPa;
      } else {
        Log.LogDebug(Formatter()
                     << "NOT updating the "
                        "amplitude (a_mean="
                     << pTZE.A_mean << ", a_sigma=" << pTZE.A_sigma);
      }
      if (std::isnan(logprior) || logprior != logprior ||
          std::isinf(logprior)) {
        THROWG(ErrorCode::INTERNAL_ERROR,
               Formatter() << "Invalid logPa value (a_mean=" << pTZE.A_mean
                           << ", a_sigma=" << pTZE.A_sigma
                           << ", precompA=" << pTZE.logprior_precompA);
      }
      result->ChiSquare[fullResultIdx] += logprior;
    }
    result->LogPrior[fullResultIdx] = logprior;
  }
}

void COperatorTemplateFittingLog::computeFitQuality(
    const std::shared_ptr<CTemplateFittingResult> &result, Int32 resultIdx,
    Int32 subResultSize, Int32 firstTplIdx, CMask const &lineMask) {

  const auto &spectrumRebinedFluxRaw =
      m_spectraFull[0]->GetFluxAxis().GetSamplesVector();
  const Int32 nSpcPixels = spectrumRebinedFluxRaw.size();
  auto const &mask = m_spectraFull[0]->getMask();
  const auto &error =
      m_spectraFull[0]->GetFluxAxis().GetError().GetSamplesVector();

  for (Int32 isubz = 0, fullResultIdx = resultIdx; isubz < subResultSize;
       ++isubz, ++fullResultIdx, ++firstTplIdx) {
    if (fullResultIdx >= ssize(result->ChiSquare))
      THROWG(ErrorCode::INTERNAL_ERROR, "out-of-bound index");

    if (m_enableIGM && result->FitMeiksinIdx[fullResultIdx] != -1)
      ApplyMeiksinCoeff(result->FitMeiksinIdx[fullResultIdx]);
    if (m_enableISM && result->FitEbmvCoeff[fullResultIdx] != -1)
      ApplyDustCoeff(
          m_templateRebined_bf.front().m_ismCorrectionCalzetti->GetEbmvIndex(
              result->FitEbmvCoeff[fullResultIdx]));

    // Compute model flux
    const auto &tplRebinedFluxRaw{
        m_templateRebined_bf[0].GetFluxAxis().GetSamplesVector()};

    TAxisSampleList modelFlux{tplRebinedFluxRaw.begin() + firstTplIdx,
                              tplRebinedFluxRaw.begin() + firstTplIdx +
                                  nSpcPixels};
    for (auto lambdaIdx = 0; lambdaIdx < ssize(modelFlux); ++lambdaIdx) {
      modelFlux[lambdaIdx] =
          modelFlux[lambdaIdx] * result->FitAmplitude[fullResultIdx];
    }

    CMask combinedMask;
    if (lineMask.GetMasksCount()) {
      CMask combinedMask(lineMask, isubz, isubz + nSpcPixels);
      combinedMask.IntersectWith(mask);
    } else {
      combinedMask = mask;
    }
    const Int32 nSpcUnmaskedPixels = combinedMask.GetUnMaskedSampleCount();

    TList<TFloat64List> spcFluxVect;
    TList<TFloat64List> modelFluxVect;
    TList<TFloat64List> spcFluxErrorVect;
    TList<CMask> maskVect;
    spcFluxVect.push_back(spectrumRebinedFluxRaw);
    modelFluxVect.push_back((std::move(modelFlux)));
    spcFluxErrorVect.push_back(error);
    maskVect.push_back(std::move(combinedMask));
    result->FitQuality[fullResultIdx] = NSFitQuality::computeFitQuality(
        spcFluxVect, modelFluxVect, spcFluxErrorVect,
        result->ChiSquare[fullResultIdx], nSpcUnmaskedPixels, maskVect);
  }
}

/**
  // TODO : many vectors allocated in this function. Check if the allocation
 time is significant, and eventually use preallocated member buffers...
 * @brief COperatorTemplateFittingLog::FitRangez
 * @param spectrumRebinedLambda
 * @param spectrumRebinedFluxRaw
 * @param error
 * @param tplRebinedLambda
 * @param tplRebinedFluxRaw
 * @param nSpc
 * @param nTpl
 * @param result
 * @param MeiksinList
 * @param EbmvList
 * @return
 */
void COperatorTemplateFittingLog::FitRangez(
    const TFloat64List &inv_err2, const TFloat64List &spcRebinedFluxOverErr2,
    const TFloat64List &spcRebinedFlux2OverErr2,
    const TInt32Range &currentRange,
    const std::shared_ptr<CTemplateFittingResult> &result,
    const TInt32List &MeiksinList, const TInt32List &EbmvList,
    const Float64 &dtd, CMask const &lineMask) {

  const TAxisSampleList &spectrumRebinedLambda =
      m_spectraFull[0]->GetSpectralAxis().GetSamplesVector();
  auto const &spcMask = m_spectraFull[0]->getMask();
  const Int32 nSpc = spectrumRebinedLambda.size();
  const Int32 nSpcUnmaskedPixels = spcMask.GetUnMaskedSampleCount();

  const TAxisSampleList &tplRebinedLambdaGlobal =
      m_templateRebined_bf[0].GetSpectralAxis().GetSamplesVector();

  Int32 kstart, kend;
  kstart = currentRange.GetBegin();
  kend = currentRange.GetEnd();
  Int32 nTpl = kend - kstart + 1;

  TFloat64List lineMaskFloat;
  TFloat64List spcMaskFloat;
  if (lineMask.GetMasksCount()) {
    lineMaskFloat = TFloat64List(lineMask.getMaskList().begin(),
                                 lineMask.getMaskList().end());

    spcMaskFloat = TFloat64List(spcMask.getMaskList().begin(),
                                spcMask.getMaskList().end());
  }

  Float64 redshiftValueMeiksin = result->Redshifts[0];

  Log.LogDebug(Formatter() << "FitRangez: redshiftValueMeiksin = "
                           << redshiftValueMeiksin);
  Log.LogDebug(Formatter() << "FitRangez: spc[0] = "
                           << spectrumRebinedLambda[0]);
  Log.LogDebug(Formatter() << "FitRangez: spc[max] = "
                           << spectrumRebinedLambda[nSpc - 1]);
  Log.LogDebug(Formatter()
               << "FitRangez: tpl[0]*zmax = "
               << tplRebinedLambdaGlobal[kstart] *
                      (1.0 + result->Redshifts[result->Redshifts.size() - 1]));
  Log.LogDebug(Formatter() << "FitRangez: tpl[max]*zmin = "
                           << tplRebinedLambdaGlobal[kstart + nTpl - 1] *
                                  (1 + result->Redshifts[0]));

  Int32 nshifts = nTpl - nSpc + 1;
  Int32 nPaddedSamples = ceil(nTpl / 2.0) * 2;

  Log.LogDetail(Formatter() << "Now fitting using the FFT on "
                               "nshifts="
                            << nshifts << " values, for Meiksin redshift="
                            << redshiftValueMeiksin);

  Log.LogDebug(Formatter() << "FitRangez: initializing FFT "
                              "with n = "
                           << nPaddedSamples << " points");
  FFTPlans fftPlans(nPaddedSamples);

  TFloat64List z_vect = result->Redshifts;
  std::reverse(z_vect.begin(), z_vect.end());
  // prepare z array
  TFloat64List z_vect_verif(nshifts, 0.0);
  Float64 relative_zgrid_error_max = 0.;
  for (Int32 t = 0; t < nshifts; t++) {
    z_vect_verif[t] =
        (spectrumRebinedLambda[0] - tplRebinedLambdaGlobal[t + kstart]) /
        tplRebinedLambdaGlobal[t + kstart];
    // compare with z_vect
    Float64 relative_zgrid_error =
        std::abs(z_vect[t] - z_vect_verif[t]) / (1 + z_vect_verif[t]);
    if (relative_zgrid_error > relative_zgrid_error_max)
      relative_zgrid_error_max = relative_zgrid_error;
  }
  Log.LogDebug(Formatter() << "FitRangez: max diff in zgrid="
                           << relative_zgrid_error_max);
  if (relative_zgrid_error_max > 5E-7) {
    THROWG(ErrorCode::INTERNAL_ERROR,
           "z_vect and z_vect_verification do not correspond.");
  }

  // check borders
  if (z_vect.size() != z_vect_verif.size())
    THROWG(ErrorCode::INTERNAL_ERROR,
           "z_vect size and z_vect_verification size do not match.");
  Int32 nISM = EbmvList.size();
  Int32 nIGM = MeiksinList.size();

  // disable IGM if the redshift range and lambda range do not make the IGM
  // wavelength appear
  Int32 enableIGM = m_enableIGM;
  Int32 overrideNIGMTobesaved = -1;
  if (tplRebinedLambdaGlobal[kstart] > RESTLAMBDA_LYA && nIGM > 1) {
    overrideNIGMTobesaved = nIGM;
    nIGM = 1;
    enableIGM = 0;
    Log.LogDebug(Formatter() << "FitRangez: IGM disabled, "
                                "min-tpl-lbda="
                             << tplRebinedLambdaGlobal[kstart]);
  }

  // prepare best fit data buffer
  TFloat64List bestChi2(nshifts, DBL_MAX);
  std::vector<TFitQuality> bestFitQuality(nshifts);
  TFloat64List bestFitAmp(nshifts, NAN);
  TFloat64List bestFitAmpErr(nshifts, NAN);
  TFloat64List bestFitAmpSigma(nshifts, NAN);
  TFloat64List bestFitDtm(nshifts, NAN);
  TFloat64List bestFitMtm(nshifts, NAN);
  TFloat64List bestFitSNR(nshifts, NAN);
  TFloat64List bestISMCoeff(nshifts, NAN);
  TInt32List bestIGMIdx(nshifts, undefIdx);

  // prepare intermediate fit data buffer
  Int32 nIGMFinal = nIGM;
  if (overrideNIGMTobesaved > nIGM) {
    nIGMFinal = overrideNIGMTobesaved;
  }
  TList<TList<TFloat64List>> intermediateChi2(
      nshifts, TList<TFloat64List>(nISM, TFloat64List(nIGMFinal, DBL_MAX)));
  TList<TInt32List> intermediateIsmEbmvIdx(nshifts, EbmvList);
  TList<TInt32List> intermediateIgmMeiksinIdx(
      nshifts, enableIGM ? MeiksinList : TInt32List(nIGMFinal, undefIdx));

  // precompute DtD and nValidSamples in case of lineMask
  // since constant for all ism/igm
  TFloat64List DtD_vec(nshifts, dtd);
  TInt32List nValidSamples_vec(nshifts, nSpcUnmaskedPixels);
  if (lineMask.GetMasksCount()) {
    // Estimate DtD if lineMask
    EstimateXtY(spcRebinedFlux2OverErr2, lineMaskFloat, DtD_vec, fftPlans,
                EPrecomputedFFT::none, EPrecomputedFFT::tplMask);
    TFloat64List nValidSamples_vec_float;
    EstimateXtY(spcMaskFloat, lineMaskFloat, nValidSamples_vec_float, fftPlans,
                EPrecomputedFFT::none, EPrecomputedFFT::tplMask);
    std::transform(nValidSamples_vec_float.cbegin(),
                   nValidSamples_vec_float.cend(), nValidSamples_vec.begin(),
                   [](Float64 v) { return std::round(v); });
  }

  // note that there is no need to copy the ism/igm cause they already exist in
  // the rebinned template
  if (m_enableIGM || m_enableISM) {
    m_templateRebined_bf[0].InitIsmIgmConfig(kstart, kend,
                                             redshiftValueMeiksin);
  }
  for (Int32 kIGM = 0; kIGM < nIGM; kIGM++) {
    if (enableIGM) {
      Log.LogDebug(Formatter() << __func__ << ": IGM index=" << kIGM);
    }

    if (enableIGM) {
      Int32 meiksinIdx = MeiksinList[kIGM];
      ApplyMeiksinCoeff(meiksinIdx);
    }

    for (Int32 kISM = 0; kISM < nISM; kISM++) {
      if (m_enableISM) {
        Log.LogDebug(Formatter() << __func__ << ": ISM index =" << kISM);
      }

      if (m_enableISM) {
        Int32 kDust = EbmvList[kISM];
        ApplyDustCoeff(kDust);
      }

      const TAxisSampleList &tplRebinedFluxcorr =
          m_templateRebined_bf[0].GetFluxAxis().GetSamplesVector();
      // extract only the relevant part
      TFloat64List::const_iterator first = tplRebinedFluxcorr.begin() + kstart,
                                   last = tplRebinedFluxcorr.begin() + kend + 1;
      const TAxisSampleList tplRebinedFluxcorr_cropped(first, last);
      TAxisSampleList tpl2RebinedFlux(nTpl);

      if (ssize(tplRebinedFluxcorr_cropped) != nTpl) {
        THROWG(ErrorCode::INTERNAL_ERROR, "vector sizes do not match");
      }
      // compute the square of the corrected flux
      for (Int32 j = 0; j < nTpl; j++) {
        tpl2RebinedFlux[j] =
            tplRebinedFluxcorr_cropped[j] * tplRebinedFluxcorr_cropped[j];
      }

      // Estimate DtM: sumCross
      TFloat64List dtm_vec;

      EstimateXtY(spcRebinedFluxOverErr2, tplRebinedFluxcorr_cropped, dtm_vec,
                  fftPlans, EPrecomputedFFT::spcFluxOverErr2);

      if (ssize(dtm_vec) != nshifts)
        THROWG(ErrorCode::INTERNAL_ERROR,
               Formatter() << "Wrong size of return cross product ("
                           << dtm_vec.size() << " instead of expected "
                           << nshifts << ")");

      // Estimate MtM: sumT
      TFloat64List mtm_vec;
      EstimateXtY(inv_err2, tpl2RebinedFlux, mtm_vec, fftPlans,
                  EPrecomputedFFT::spcOneOverErr2);

      Log.LogDebug(Formatter() << __func__ << ": dtd = " << dtd);

      // Estimate Chi2
      if (ssize(mtm_vec) != nshifts) {
        THROWG(ErrorCode::INTERNAL_ERROR,
               Formatter() << "xty vector size do not match: dtm size = "
                           << nshifts << ", mtm size =" << mtm_vec.size());
      }
      TFloat64List chi2(nshifts, DBL_MAX);
      TFloat64List amp(nshifts, DBL_MAX);
      TFloat64List amp_sigma(nshifts);
      TFloat64List amp_err(nshifts, DBL_MAX);
      for (Int32 k = 0; k < nshifts; k++) {
        if (mtm_vec[k] == 0.0) {
          amp[k] = 0.0;
          amp_err[k] = 0.0;
          amp_sigma[k] = 0.0;
          chi2[k] = dtd; // keep at maximum
        } else {
          amp[k] = dtm_vec[k] / mtm_vec[k];
          amp_err[k] = sqrt(1. / mtm_vec[k]);
          amp_sigma[k] = amp[k] / amp_err[k];
          applyPositiveAndNonNullConstraint(amp_sigma[k], amp[k]);
          chi2[k] = DtD_vec[k] - 2 * dtm_vec[k] * amp[k] +
                    mtm_vec[k] * amp[k] * amp[k];
        }
      }

      for (Int32 k = 0; k < nshifts; k++) {
        intermediateChi2[k][kISM][kIGM] = chi2[k];
        // in the case of 1215A is not in the range, no need to
        // recompute with varying IGM coeff.
        if (overrideNIGMTobesaved > 1 && kIGM == 0) {
          for (Int32 koigm = 1; koigm < overrideNIGMTobesaved; koigm++) {
            intermediateChi2[k][kISM][koigm] = chi2[k];
          }
        }

        if (bestChi2[k] > chi2[k]) {
          bestChi2[k] = chi2[k];
          bestFitQuality[k].reducedChiSquare =
              NSFitQuality::reducedChi2(chi2[k], nValidSamples_vec[k]);
          bestFitQuality[k].pValue =
              NSFitQuality::pValue(chi2[k], nValidSamples_vec[k]);
          bestFitAmp[k] = amp[k];
          bestFitAmpErr[k] = amp_err[k];
          bestFitAmpSigma[k] = amp_sigma[k];
          bestFitDtm[k] = dtm_vec[k];
          bestFitMtm[k] = mtm_vec[k];
          bestFitSNR[k] = -1.;
          if (bestFitMtm[k] > 0) {
            bestFitSNR[k] = bestFitDtm[k] / std::sqrt(bestFitMtm[k]);
          }
          bestISMCoeff[k] =
              m_enableISM
                  ? m_templateRebined_bf[0]
                        .m_ismCorrectionCalzetti->GetEbmvValue(EbmvList[kISM])
                  : undefIdx;
          bestIGMIdx[k] = enableIGM ? MeiksinList[kIGM] : undefIdx;
        }
      }

      Log.LogDebug(Formatter()
                   << __func__ << ": spc lbda 0 =" << spectrumRebinedLambda[0]);
      Log.LogDebug(Formatter() << __func__ << ": tpl lbda 0 ="
                               << tplRebinedLambdaGlobal[kstart]);
      Float64 z_O =
          (spectrumRebinedLambda[0] - tplRebinedLambdaGlobal[kstart]) /
          tplRebinedLambdaGlobal[kstart];
      Log.LogDebug(Formatter() << __func__ << ": z 0 =" << z_O);
    }
  }

  // reversing all vectors
  std::reverse(z_vect.begin(), z_vect.end());
  std::reverse(bestChi2.begin(), bestChi2.end());
  std::reverse(bestFitQuality.begin(), bestFitQuality.end());
  std::reverse(bestFitAmp.begin(), bestFitAmp.end());
  std::reverse(bestFitAmpErr.begin(), bestFitAmpErr.end());
  std::reverse(bestFitAmpSigma.begin(), bestFitAmpSigma.end());
  std::reverse(bestFitDtm.begin(), bestFitDtm.end());
  std::reverse(bestFitMtm.begin(), bestFitMtm.end());
  std::reverse(bestFitSNR.begin(), bestFitSNR.end());
  std::reverse(bestISMCoeff.begin(), bestISMCoeff.end());
  std::reverse(bestIGMIdx.begin(), bestIGMIdx.end());
  std::reverse(intermediateChi2.begin(), intermediateChi2.end());
  for (Int32 k = 0; k < ssize(result->Redshifts); k++)
    result->Overlap[k] = TFloat64List(1, 1.0);

  result->ChiSquare = bestChi2;
  result->FitQuality = bestFitQuality;
  result->FitAmplitude = bestFitAmp;
  result->FitAmplitudeError = bestFitAmpErr;
  result->FitAmplitudeSigma = bestFitAmpSigma;
  result->FitDtM = bestFitDtm;
  result->FitMtM = bestFitMtm;
  result->SNR = bestFitSNR;
  result->FitEbmvCoeff = bestISMCoeff;
  result->FitMeiksinIdx = bestIGMIdx;
  result->ChiSquareIntermediate = intermediateChi2;
  // no need to reverse the two next: all values identical along z
  result->IsmEbmvIdxIntermediate = intermediateIsmEbmvIdx;
  result->IgmMeiksinIdxIntermediate = intermediateIgmMeiksinIdx;
}

// find indexes in templateSpectra for which Z falls into the redshift range
/**
 * Method logic:
 * 1. count the number of steps (Integer) between start/end borders of
 * redshiftRange and the min/max Z values corresponding to the
 * templateSpectralAxis min/max borders
 * 2. this number of zsteps correspond exactly to the right/left offsets on the
 * spectralAxis in log
 */
TInt32Range COperatorTemplateFittingLog::FindTplSpectralIndex(
    const TFloat64Range &redshiftrange) const {
  return FindTplSpectralIndex(m_spectraFull[0]->GetSpectralAxis(),
                              m_templateRebined_bf[0].GetSpectralAxis(),
                              redshiftrange);
}

TInt32Range COperatorTemplateFittingLog::FindTplSpectralIndex(
    const CSpectrumSpectralAxis &spcSpectralAxis,
    const CSpectrumSpectralAxis &tplSpectralAxis,
    const TFloat64Range &redshiftrange) const {

  const TFloat64Range spcRange =
      spcSpectralAxis
          .GetLambdaRange(); // this is correct only if spcSpectralAxis is
                             // intersected with the input lambdaRange
  const TFloat64Range tplRange = tplSpectralAxis.GetLambdaRange();
  const Int32 tplsize = tplSpectralAxis.GetSamplesCount();
  const Float64 logstep = tplSpectralAxis.GetlogGridStep();

  Float64 zmax =
      (spcRange.GetBegin() - tplRange.GetBegin()) /
      tplRange
          .GetBegin(); // get maximum z reachable with given template & spectra
  Float64 offset_left =
      log((zmax + 1) / (redshiftrange.GetEnd() + 1)) /
      logstep; //  should be integer at numerical precision before round
  Int32 ilbdamin = round(offset_left); // deduce min lambda from max reachable z
                                       // and max z in current range.
  Float64 zmin =
      (spcRange.GetEnd() - tplRange.GetEnd()) /
      tplRange
          .GetEnd(); // get minimum reachable z with given template & spectra

  Float64 offset_right =
      log((redshiftrange.GetBegin() + 1) / (zmin + 1)) / logstep;
  Int32 ilbdamax =
      tplsize - 1 - round(offset_right); // deduce max lambda from min reachable
                                         // z and min min z in current range.

  if (ilbdamax >= tplsize || ilbdamin < 0)
    THROWG(ErrorCode::INTERNAL_ERROR, " Failed to find indexes");

  if (ilbdamin > ilbdamax) {
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Problem with "
                          "tpl indexes for zranges, found ilbdamin="
                       << ilbdamin << " > ilbdamax=" << ilbdamax);
  }

  return TInt32Range(ilbdamin, ilbdamax);
}
/**
 * \brief COperatorTemplateFittingLog::Compute
 *
 * This method computes the log_likelihood for the input spc and the tpl on a
 *given redshift range (should be a regular grid): 0. checks :
 *      - is the redshift list a regular grid ? (assumed in the rest of the
 *method)
 *      - is the overlap always >100% in the given redshift range ?
 * 1. resample the input spectrum/tpl on a loglambda regular grid (always
 *resampling as of 2017-06-13, option todo: use the already log-sampled input
 *spectrum grid)
 *
 * lambdaRange is not clamped
 **/
std::shared_ptr<CTemplateFittingResult> COperatorTemplateFittingLog::Compute(
    const CTemplate &logSampledTpl, Float64 overlapThreshold,
    std::string opt_interp, bool opt_extinction, bool opt_dustFitting,
    Float64 opt_continuum_null_amp_threshold,
    const CPriorHelper::TPriorZEList &logpriorze, Int32 FitEbmvIdx,
    Int32 FitMeiksinIdx, TInt32Range zIdxRangeToCompute,
    std::shared_ptr<CTemplateFittingResult> const &dummyResult) {
  Log.LogDetail(Formatter() << "starting computation for template: "
                            << logSampledTpl.GetName());

  if (opt_dustFitting && logSampledTpl.CalzettiInitFailed())
    THROWG(ErrorCode::INTERNAL_ERROR, "ISM is no initialized");

  if (opt_dustFitting &&
      FitEbmvIdx >=
          logSampledTpl.m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs()) {
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << " Invalid calzetti index: (FitEbmvIdx=" << FitEbmvIdx
                       << ",while NPrecomputedEbmvCoeffs="
                       << logSampledTpl.m_ismCorrectionCalzetti
                              ->GetNPrecomputedEbmvCoeffs()
                       << ")");
  }

  if (opt_extinction && logSampledTpl.MeiksinInitFailed()) {
    THROWG(ErrorCode::INTERNAL_ERROR, "IGM is not initialized");
  }

  if (!logSampledTpl.GetSpectralAxis().IsLogSampled()) {
    THROWG(ErrorCode::INTERNAL_ERROR, "template is not log sampled");
  }
  // check if spc and tpl have same step
  const Float64 epsilon = 1E-8;
  if (std::abs(m_spectraFull[0]->GetSpectralAxis().GetlogGridStep() -
               logSampledTpl.GetSpectralAxis().GetlogGridStep() * m_ssRatio) >
      epsilon)
    THROWG(ErrorCode::INTERNAL_ERROR,
           "tpl and spc are not sampled with the same step");

  m_continuum_null_amp_threshold = opt_continuum_null_amp_threshold;

  // subsample template if necessary
  if (m_ssRatio == 1) { // no required subsampling
    m_templateRebined_bf[0] = logSampledTpl;
  } else {
    TInt32Range ilbda = FindTplSpectralIndex(
        m_spectraFull[0]->GetSpectralAxis(), logSampledTpl.GetSpectralAxis(),
        TFloat64Range(m_redshifts));
    TMaskList mask_tpl =
        logSampledTpl.GetSpectralAxis().GetSubSamplingMask(m_ssRatio, ilbda);

    m_templateRebined_bf[0] = CTemplate(logSampledTpl, mask_tpl);
    // double make sure that subsampled spectrum is well sampled
    if (!m_templateRebined_bf[0].GetSpectralAxis().IsLogSampled(m_logstep)) {
      THROWG(ErrorCode::INTERNAL_ERROR,
             "subsampled template "
             "is not log sampled with the redshift step");
    }
  }

  //**************** Fitting at all redshifts ****************//
  // Note: below corresponds to ::BasicFit code except that redshift loop
  // belongs to ::compute
  // Optionally apply some IGM absorption
  TIgmIsmIdxs igmIsmIdxs = m_templateRebined_bf.front().GetIsmIgmIdxList(
      opt_extinction, opt_dustFitting, FitEbmvIdx, FitMeiksinIdx);
  Int32 nIGMCoeffs = igmIsmIdxs.igmIdxs.size();
  Int32 nISMCoeffs = igmIsmIdxs.ismIdxs.size();

  m_enableIGM = opt_extinction;
  m_enableISM = opt_dustFitting;
  auto const result = std::make_shared<CTemplateFittingResult>(
      m_redshifts.size(), nISMCoeffs, nIGMCoeffs);
  result->Redshifts = m_redshifts;

  if (logpriorze.size() > 0 && logpriorze.size() != m_redshifts.size()) {
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter()
               << "Vector size do not match between logpriorz and redshift: "
               << logpriorze.size() << " != " << m_redshifts.size());
  }

  CMask lineMask;
  if (!m_maskBuilder->isDefaultMask())
    lineMask = maskTemplate();

  FitAllz(result, igmIsmIdxs.igmIdxs, igmIsmIdxs.ismIdxs, logpriorze, lineMask);

  //**************** End Fitting at all redshifts ****************//

  // overlap warning
  Float64 overlapValidInfZ = -1;
  for (Int32 i = 0; i < ssize(m_redshifts); i++) {
    if (result->Overlap[i].front() >= overlapThreshold) {
      overlapValidInfZ = m_redshifts[i];
      break;
    }
  }
  Float64 overlapValidSupZ = -1;
  for (Int32 i = m_redshifts.size() - 1; i >= 0; i--) {
    if (result->Overlap[i].front() >= overlapThreshold) {
      overlapValidSupZ = m_redshifts[i];
      break;
    }
  }
  if (overlapValidInfZ != m_redshifts.front() ||
      overlapValidSupZ != m_redshifts.back()) {
    Log.LogInfo(Formatter() << "overlap warning for " << logSampledTpl.GetName()
                            << ": minz=" << overlapValidInfZ
                            << ", maxz=" << overlapValidSupZ);
  }

  // estimate CstLog for PDF estimation
  //   note: with linemask it is not a constant anymore, since the number of
  //   valid pixels depends on z, BUT this quantity is never used in linemodel.
  //    It is used in templatefittingSolve, which has no linemask.
  result->CstLog = EstimateLikelihoodCstLog();

  return result;
}

CMask COperatorTemplateFittingLog::maskTemplate() {
  auto const &spectralAxis = m_templateRebined_bf[0].GetSpectralAxis();
  auto fluxAxis = m_templateRebined_bf[0].GetFluxAxis();

  auto const range = spectralAxis.GetLambdaRange();
  auto const mask = m_maskBuilder->getMask(spectralAxis, range, 0.0, 0);
  for (Int32 i = 0; i != spectralAxis.GetSamplesCount(); ++i) {
    if (!mask[i])
      fluxAxis[i] = 0.0;
  }
  m_templateRebined_bf[0].SetFluxAxis(std::move(fluxAxis));
  return mask;
}

Float64 COperatorTemplateFittingLog::EstimateLikelihoodCstLog() const {
  Float64 cstLog = 0.0;
  for (auto const &[spectrum_ptr, lambdaRange_ptr] :
       boost::combine(m_spectraFull, m_lambdaRanges)) {
    const CSpectrumSpectralAxis &spcSpectralAxis =
        spectrum_ptr->GetSpectralAxis();
    const TFloat64List &error =
        spectrum_ptr->GetFluxAxis().GetError().GetSamplesVector();
    auto const &mask = spectrum_ptr->getMask();

    Int32 numDevs = 0;

    Float64 sumLogNoise = 0.0;

    Int32 imin;
    Int32 imax;
    lambdaRange_ptr->getClosedIntervalIndices(
        spcSpectralAxis.GetSamplesVector(), imin, imax);
    for (Int32 j = imin; j <= imax; j++) {
      if (mask[j]) {
        numDevs++;
        sumLogNoise += log(error[j]);
      }
    }
    cstLog += -numDevs * 0.5 * log(2 * M_PI) - sumLogNoise;
  }
  return cstLog;
}
