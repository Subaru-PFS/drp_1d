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
#include "RedshiftLibrary/operator/templatefittinglog.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/operator/templatefittingresult.h"
#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <climits>

using namespace NSEpic;
using namespace std;

COperatorTemplateFittingLog::COperatorTemplateFittingLog(
    const CSpectrum &spectrum, const CSpectrum &logSampledSpectrum,
    const TFloat64Range &lambdaRange, const TFloat64List &redshifts)
    : COperatorTemplateFittingBase(spectrum, lambdaRange, redshifts),
      m_logSampledSpectrum(logSampledSpectrum) {
  if (!m_logSampledSpectrum.GetSpectralAxis().IsLogSampled()) {
    THROWG(INTERNAL_ERROR, "input spectrum is not log sampled");
  }

  m_logSampledSpectrum.GetSpectralAxis().ClampLambdaRange(
      lambdaRange, m_logSampledLambdaRange);

  CheckRedshifts();
};

void COperatorTemplateFittingLog::SetRedshifts(TFloat64List redshifts) {
  COperatorTemplateFittingBase::SetRedshifts(std::move(redshifts));
  CheckRedshifts();
}

void COperatorTemplateFittingLog::CheckRedshifts() {
  if (m_redshifts.size() < 2) {
    THROWG(INTERNAL_ERROR, Formatter() << "Invalid redshift array: "
                                       << m_redshifts.size() << " <2");
  }

  // check if spectrum sampling is a multiple of redshift step
  m_logstep = log((m_redshifts[1] + 1) / (m_redshifts[0] + 1));
  Float64 modulo;
  m_ssRatio = m_logSampledSpectrum.GetSpectralAxis().GetLogSamplingIntegerRatio(
      m_logstep, modulo);
  if (std::abs(modulo) > 1E-12)
    THROWG(INTERNAL_ERROR, "spc and tpl do not have a lambdastep "
                           "multiple of redshift step");

  // subsample spectrum if necessary (coarse redshift grid)
  if (m_ssRatio == 1) { // no required subsampling
    m_ssSpectrum = m_logSampledSpectrum;
  } else {
    TFloat64List mask_spc =
        m_logSampledSpectrum.GetSpectralAxis().GetSubSamplingMask(
            m_ssRatio, m_logSampledLambdaRange);
    m_ssSpectrum = CSpectrum(m_logSampledSpectrum, mask_spc);
    // scale the variance by ssratio
    CSpectrumNoiseAxis scaledNoise = m_ssSpectrum.GetErrorAxis();
    scaledNoise *= 1. / sqrt(m_ssRatio);
    m_ssSpectrum.SetErrorAxis(std::move(scaledNoise));

    // double make sure that subsampled spectrum is well sampled
    if (!m_ssSpectrum.GetSpectralAxis().IsLogSampled(m_logstep)) {
      THROWG(INTERNAL_ERROR, "subsampled spectrum is not correctly log sampled "
                             "at the redshift step");
    }
  }

  // set the lambda range to the clamped subsampled spectrum
  m_ssSpectrum.GetSpectralAxis().ClampLambdaRange(m_logSampledLambdaRange,
                                                  m_ssLambdaRange);
}

COperatorTemplateFittingLog::~COperatorTemplateFittingLog() { freeFFTPlans(); }

Int32 COperatorTemplateFittingLog::EstimateXtYSlow(const TFloat64List &X,
                                                   const TFloat64List &Y,
                                                   Int32 nShifts,
                                                   TFloat64List &XtY) {
  XtY.resize(nShifts);

  Int32 nX = X.size();
  Float64 xty = 0.0;
  for (Int32 k = 0; k < nShifts; k++) {
    xty = 0.0;
    for (Int32 j = 0; j < nX; j++) {
      xty += X[j] * Y[j + k];
    }
    XtY[k] = xty;
  }
  return 0;
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

void COperatorTemplateFittingLog::EstimateXtY(const TFloat64List &X,
                                              const TFloat64List &Y,
                                              Int32 nshifts, TFloat64List &XtY,
                                              Int32 precomputedFFT) {
  // Processing the FFT
  Int32 nSpc = X.size();
  Int32 nTpl = Y.size();
  Int32 nPadded = m_nPaddedSamples;

  Int32 nPadBeforeSpc = nPadded - nSpc; //(Int32)nPadded/2.0;
  Int32 nPadBeforeTpl = 0;

  Log.LogDebug("FitAllz: Processing spc-fft "
               "with n=%d, padded to n=%d",
               nSpc, nPadded);

  bool computeSpcFFT = true;
  if ((precomputedFFT == 0 && precomputedFFT_spcFluxOverErr2 != 0) ||
      (precomputedFFT == 1 && precomputedFFT_spcOneOverErr2 != 0))
    computeSpcFFT = false;

  if (computeSpcFFT) {
    Log.LogDebug("FitAllz: Processing spc-fft with nPadBeforeSpc=%d",
                 nPadBeforeSpc);
    for (Int32 k = 0; k < nPadBeforeSpc; k++)
      inSpc[k] = 0.0;
    for (Int32 k = nPadBeforeSpc; k < nPadBeforeSpc + nSpc; k++)
      inSpc[k] = X[nSpc - 1 - (k - nPadBeforeSpc)]; // X[k-nPadBeforeSpc];
    for (Int32 k = nPadBeforeSpc + nSpc; k < nPadded; k++)
      inSpc[k] = 0.0;

    fftw_execute(pSpc);
    Log.LogDebug("FitAllz: spc-fft done");

    // save computed FFT into precomputed buffer
    if (precomputedFFT == 0) {
      precomputedFFT_spcFluxOverErr2 =
          (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPadded);
      if (precomputedFFT_spcFluxOverErr2 == 0)
        THROWG(INTERNAL_ERROR, "Unable to "
                               "allocate precomputedFFT_spcFluxOverErr2");
      for (Int32 k = 0; k < nPadded; k++) {
        precomputedFFT_spcFluxOverErr2[k][0] = outSpc[k][0];
        precomputedFFT_spcFluxOverErr2[k][1] = outSpc[k][1];
      }
    }
    if (precomputedFFT == 1) {
      precomputedFFT_spcOneOverErr2 =
          (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPadded);
      if (precomputedFFT_spcOneOverErr2 == 0)
        THROWG(INTERNAL_ERROR, "Unable to "
                               "allocate precomputedFFT_spcOneOverErr2");
      for (Int32 k = 0; k < nPadded; k++) {
        precomputedFFT_spcOneOverErr2[k][0] = outSpc[k][0];
        precomputedFFT_spcOneOverErr2[k][1] = outSpc[k][1];
      }
    }

  } else {
    if (precomputedFFT == 0)
      for (Int32 k = 0; k < nPadded; k++) {
        outSpc[k][0] = precomputedFFT_spcFluxOverErr2[k][0];
        outSpc[k][1] = precomputedFFT_spcFluxOverErr2[k][1];
      }
    if (precomputedFFT == 1)
      for (Int32 k = 0; k < nPadded; k++) {
        outSpc[k][0] = precomputedFFT_spcOneOverErr2[k][0];
        outSpc[k][1] = precomputedFFT_spcOneOverErr2[k][1];
      }
  }

  Log.LogDebug("FitAllz: Processing tpl-fft with n=%d, padded to n=%d", nTpl,
               nPadded);

  for (Int32 k = 0; k < nPadBeforeTpl; k++) {
    inTpl_padded[k] = 0.0;
  }
  for (Int32 k = nPadBeforeTpl; k < nPadBeforeTpl + nTpl; k++) {
    inTpl_padded[k] = Y[k - nPadBeforeTpl];
  }
  for (Int32 k = nPadBeforeTpl + nTpl; k < nPadded; k++) {
    inTpl_padded[k] = 0.0;
  }

  for (Int32 k = 0; k < nPadded; k++) {
    inTpl[k] = inTpl_padded[k]; // inTpl_padded[nPadded-1-k][0];//
  }

  fftw_execute(pTpl);
  Log.LogDebug("FitAllz: tpl-fft done");

  // Multiplying the FFT outputs
  for (Int32 k = 0; k < nPadded; k++) {
    outCombined[k][0] =
        (outTpl[k][0] * outSpc[k][0] - outTpl[k][1] * outSpc[k][1]);
    outCombined[k][1] =
        (outTpl[k][0] * outSpc[k][1] + outTpl[k][1] * outSpc[k][0]);
    // Y conjugate
    // outCombined[k][0] =
    // (outTpl[k][0]*outSpc[k][0]+outTpl[k][1]*outSpc[k][1]);
    // outCombined[k][1] =
    // (outTpl[k][0]*outSpc[k][1]-outTpl[k][1]*outSpc[k][0]);
  }

  fftw_execute(pBackward);
  Log.LogDebug("FitAllz: backward-fft done");

  XtY.resize(nshifts);
  Int32 offsetSamples = 0; // nPadded/2.0;//nSpc;//(nTpl+nSpc);
  for (Int32 k = 0; k < nshifts; k++) {
    // XtY[k] = inCombined[k+(Int32)(nPadded/2.0)][0]/nPadded;
    Int32 IndexCyclic = k + offsetSamples;
    if (IndexCyclic < 0)
      IndexCyclic += nPadded;
    else if (IndexCyclic >= nPadded)
      IndexCyclic -= nPadded;

    XtY[k] = inCombined[IndexCyclic] / (Float64)nPadded;
  }
}

Int32 COperatorTemplateFittingLog::InitFFT(Int32 nPadded) {
  freeFFTPlans();

  inSpc = (Float64 *)fftw_malloc(sizeof(Float64) * nPadded);
  outSpc = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPadded);
  pSpc = fftw_plan_dft_r2c_1d(nPadded, inSpc, outSpc, FFTW_ESTIMATE);
  if (inSpc == 0) {
    THROWG(INTERNAL_ERROR, "Unable to allocate inSpc");
  }
  if (outSpc == 0) {
    THROWG(INTERNAL_ERROR, "Unable to allocate outSpc");
  }

  inTpl_padded = (Float64 *)fftw_malloc(sizeof(Float64) * nPadded);
  inTpl = (Float64 *)fftw_malloc(sizeof(Float64) * nPadded);
  outTpl = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPadded);
  pTpl = fftw_plan_dft_r2c_1d(nPadded, inTpl, outTpl, FFTW_ESTIMATE);
  if (inTpl_padded == 0) {
    THROWG(INTERNAL_ERROR, "Unable to "
                           "allocate inTpl_padded");
  }
  if (inTpl == 0) {
    THROWG(INTERNAL_ERROR, "Unable to allocate inTpl");
  }
  if (outTpl == 0) {
    THROWG(INTERNAL_ERROR, "Unable to allocate outTpl");
  }

  outCombined = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPadded);
  inCombined = (Float64 *)fftw_malloc(sizeof(Float64) * nPadded);
  pBackward =
      fftw_plan_dft_c2r_1d(nPadded, outCombined, inCombined, FFTW_ESTIMATE);
  if (outCombined == 0) {
    THROWG(INTERNAL_ERROR, "Unable to "
                           "allocate outCombined");
  }
  if (inCombined == 0) {
    THROWG(INTERNAL_ERROR, "Unable to "
                           "allocate inCombined");
  }

  // reinit fft precomputed buffers
  freeFFTPrecomputedBuffers();

  return 0;
}

void COperatorTemplateFittingLog::freeFFTPrecomputedBuffers() {
  if (precomputedFFT_spcFluxOverErr2 != NULL) {
    fftw_free(precomputedFFT_spcFluxOverErr2);
    precomputedFFT_spcFluxOverErr2 = 0;
  }
  if (precomputedFFT_spcOneOverErr2 != NULL) {
    fftw_free(precomputedFFT_spcOneOverErr2);
    precomputedFFT_spcOneOverErr2 = 0;
  }
}

void COperatorTemplateFittingLog::freeFFTPlans() {
  if (pSpc) {
    fftw_destroy_plan(pSpc);
    pSpc = 0;
  }
  if (inSpc) {
    fftw_free(inSpc);
    inSpc = 0;
  }
  if (outSpc) {
    fftw_free(outSpc);
    outSpc = 0;
  }
  if (pTpl) {
    fftw_destroy_plan(pTpl);
    pTpl = 0;
  }
  if (inTpl_padded) {
    fftw_free(inTpl_padded);
    inTpl_padded = 0;
  }
  if (inTpl) {
    fftw_free(inTpl);
    inTpl = 0;
  }
  if (outTpl) {
    fftw_free(outTpl);
    outTpl = 0;
  }
  if (pBackward) {
    fftw_destroy_plan(pBackward);
    pBackward = 0;
  }
  if (inCombined) {
    fftw_free(inCombined);
    inCombined = 0;
  }
  if (outCombined) {
    fftw_free(outCombined);
    outCombined = 0;
  }

  freeFFTPrecomputedBuffers();
}

TInt32RangeList
COperatorTemplateFittingLog::FindZRanges(const TFloat64List &redshifts) {
  TInt32RangeList izrangelist;
  TInt32List zsplit;
  if (m_enableIGM) {
    Float64 zmin_igm =
        GetIGMStartingRedshiftValue(m_ssSpectrum.GetSpectralAxis()[0]);
    if (zmin_igm > redshifts.front() && zmin_igm < redshifts.back()) {
      Int32 i_zmin_igm = -1;
      TFloat64Index::getClosestLowerIndex(redshifts, zmin_igm, i_zmin_igm);
      zsplit.push_back(i_zmin_igm);
    }
    if (zmin_igm < redshifts.back()) {
      TFloat64List zbins_igm =
          m_templateRebined_bf.m_igmCorrectionMeiksin->getRedshiftBins();
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

    Log.LogDebug("FitAllz: indexes for full "
                 "LstSquare calculation, count = %d",
                 zsplit.size());
    for (Int32 k = 0; k < zsplit.size(); k++) {
      Log.LogDebug("FitAllz: indexes ranges: "
                   "for i=%d, zsplit=%d",
                   k, zsplit[k]);
    }
  }

  Int32 izmin = 0;
  for (const Int32 &izmax : zsplit) {
    izrangelist.push_back(TInt32Range(izmin, izmax));
    izmin = izmax + 1;
  }
  Int32 izmax = redshifts.size() - 1;
  izrangelist.push_back(TInt32Range(izmin, izmax));

  Log.LogDebug("FitAllz: indexes - "
               "izrangelist calculation, count = %d",
               izrangelist.size());
  for (Int32 k = 0; k < izrangelist.size(); k++) {
    Log.LogDebug("FitAllz: indexes ranges: "
                 "for i=%d, zmin=%f, zmax=%f",
                 k, redshifts[izrangelist[k].GetBegin()],
                 redshifts[izrangelist[k].GetEnd()]);
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
Int32 COperatorTemplateFittingLog::FitAllz(
    std::shared_ptr<CTemplateFittingResult> result,
    const TInt32List &MeiksinList, const TInt32List &EbmvList,
    const CPriorHelper::TPriorZEList &logpriorze) {

  // prepare list of redshifts segments that keep the same IGM extinction curves
  TInt32RangeList izrangelist = FindZRanges(result->Redshifts);
  Int32 nzranges = izrangelist.size();

  // since dtd is cte, better compute it here
  const TAxisSampleList &error =
      m_ssSpectrum.GetFluxAxis().GetError().GetSamplesVector();
  ;
  const TAxisSampleList &spectrumRebinedFluxRaw =
      m_ssSpectrum.GetFluxAxis().GetSamplesVector();
  Float64 dtd = 0.0;
  TFloat64List inv_err2(error.size());
  for (Int32 j = 0; j < error.size(); j++) {
    inv_err2[j] = 1.0 / (error[j] * error[j]);
    dtd += spectrumRebinedFluxRaw[j] * spectrumRebinedFluxRaw[j] * inv_err2[j];
  }

  for (Int32 k = 0; k < nzranges; k++) {
    // prepare the zrange-result container
    std::shared_ptr<CTemplateFittingResult> subresult;
    TFloat64Range zrange =
        TFloat64Range(result->Redshifts[izrangelist[k].GetBegin()],
                      result->Redshifts[izrangelist[k].GetEnd()]);
    if (m_enableIGM && result->Redshifts.size() > 1) {
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
          result->Redshifts.size(), EbmvList.size(), MeiksinList.size());
      subresult->Redshifts = result->Redshifts;
    }

    TInt32Range ilbda = FindTplSpectralIndex(zrange);
    Log.LogDebug("FitAllz: zrange min=%f, max=%f", zrange.GetBegin(),
                 zrange.GetEnd());
    Log.LogDebug("FitAllz: full zmin=%f, full zmax=%f", result->Redshifts[0],
                 result->Redshifts[result->Redshifts.size() - 1]);
    Log.LogDebug("FitAllz: indexes tpl crop: "
                 "lbda min=%d, max=%d",
                 ilbda.GetBegin(), ilbda.GetEnd());
    Log.LogDebug("FitAllz: indexes tpl full: "
                 "lbda min=%d, max=%d",
                 0, m_templateRebined_bf.GetSampleCount() - 1);
    /*Log.LogDebug("FitAllz: tpl lbda "
                 "min*zmax=%f, max*zmin=%f",
                 tplRebinedLambdaGlobal[ilbda.GetBegin()] * (1.0 +
       zrange.GetEnd()), tplRebinedLambdaGlobal[ilbda.GetEnd()] * (1.0 +
       zrange.GetBegin()));*/

    FitRangez(inv_err2, ilbda, subresult, MeiksinList, EbmvList, dtd);

    // copy subresults into global results
    for (Int32 isubz = 0; isubz < subresult->Redshifts.size(); isubz++) {
      Int32 fullResultIdx = isubz + izrangelist[k].GetBegin();
      if (fullResultIdx >= result->ChiSquare.size())
        THROWG(INTERNAL_ERROR, "out-of-bound index");
      result->ChiSquare[fullResultIdx] = subresult->ChiSquare[isubz];
      result->FitAmplitude[fullResultIdx] = subresult->FitAmplitude[isubz];
      result->FitAmplitudeError[fullResultIdx] =
          subresult->FitAmplitudeError[isubz];
      result->FitAmplitudeSigma[fullResultIdx] =
          subresult->FitAmplitudeSigma[isubz];
      result->FitDtM[fullResultIdx] = subresult->FitDtM[isubz];
      result->FitMtM[fullResultIdx] = subresult->FitMtM[isubz];

      Float64 logprior = 0.;
      if (logpriorze.size() > 0) {
        Int32 kism_best = 0;
        if (subresult->FitEbmvCoeff[isubz] != -1.0)
          kism_best =
              m_templateRebined_bf.m_ismCorrectionCalzetti->GetEbmvIndex(
                  subresult->FitEbmvCoeff[isubz]);

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
            ampl =
                result->FitDtM[fullResultIdx] / result->FitMtM[fullResultIdx];
            ampl_err = sqrt(1. / result->FitMtM[fullResultIdx]);
          }

          Log.LogDebug("update the amplitude "
                       "(a_mean=%e, a_sigma=%e)",
                       pTZE.A_mean, pTZE.A_sigma);
          Log.LogDebug("update the amplitude "
                       "(ampl was = %e, updated to %e)",
                       result->FitAmplitude[fullResultIdx], ampl);

          // check negative amplitude
          ampl_sigma = ampl / ampl_err;
          applyPositiveAndNonNullConstraint(ampl_sigma, ampl);

          result->FitAmplitude[fullResultIdx] = ampl;
          result->FitAmplitudeError[fullResultIdx] = ampl_err;
          result->FitAmplitudeSigma[fullResultIdx] = ampl_sigma;
          result->ChiSquare[fullResultIdx] =
              dtd + result->FitMtM[fullResultIdx] * ampl * ampl -
              2. * ampl * result->FitDtM[fullResultIdx];

          Float64 logPa = pTZE.betaA * (ampl - pTZE.A_mean) *
                          (ampl - pTZE.A_mean) / (pTZE.A_sigma * pTZE.A_sigma);
          if (std::isnan(logPa) || logPa != logPa || std::isinf(logPa)) {
            THROWG(INTERNAL_ERROR,
                   Formatter() << " Invalid logPa value (a_mean=" << pTZE.A_mean
                               << ", a_sigma=" << pTZE.A_sigma);
          }
          logprior += logPa;
        } else {
          Log.LogDebug("NOT updating the "
                       "amplitude (a_mean=%e, a_sigma=%e)",
                       pTZE.A_mean, pTZE.A_sigma);
        }
        if (std::isnan(logprior) || logprior != logprior ||
            std::isinf(logprior)) {
          THROWG(INTERNAL_ERROR,
                 Formatter() << "Invalid logPa value (a_mean=" << pTZE.A_mean
                             << ", a_sigma=" << pTZE.A_sigma
                             << ", precompA=" << pTZE.logprior_precompA);
        }
        result->ChiSquare[fullResultIdx] += logprior;
      }
      result->Overlap[fullResultIdx] = subresult->Overlap[isubz];
      result->LogPrior[fullResultIdx] = logprior;
      result->FitEbmvCoeff[fullResultIdx] = subresult->FitEbmvCoeff[isubz];
      result->FitMeiksinIdx[fullResultIdx] = subresult->FitMeiksinIdx[isubz];
      result->Status[fullResultIdx] = subresult->Status[isubz];

      for (Int32 kism = 0;
           kism < result->ChiSquareIntermediate[fullResultIdx].size(); kism++) {
        for (Int32 kigm = 0;
             kigm < result->ChiSquareIntermediate[fullResultIdx][kism].size();
             kigm++) {
          result->ChiSquareIntermediate[fullResultIdx][kism][kigm] =
              subresult->ChiSquareIntermediate[isubz][kism][kigm];
          if (logpriorze.size() > 0) {
            Float64 logprior =
                0.; // not implemented -> not a problem for fullmodel, but will
                    // be necessary for tplmodel method for example
            result->ChiSquareIntermediate[fullResultIdx][kism][kigm] +=
                logprior;
          }
        }
      }
    }
  }

  return 0;
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
Int32 COperatorTemplateFittingLog::FitRangez(
    const TFloat64List &inv_err2, const TInt32Range &currentRange,
    const std::shared_ptr<CTemplateFittingResult> &result,
    const TInt32List &MeiksinList, const TInt32List &EbmvList,
    const Float64 &dtd) {
  const TAxisSampleList &spectrumRebinedLambda =
      m_ssSpectrum.GetSpectralAxis().GetSamplesVector();
  const TAxisSampleList &spectrumRebinedFluxRaw =
      m_ssSpectrum.GetFluxAxis().GetSamplesVector();
  Int32 nSpc = spectrumRebinedLambda.size();

  const TAxisSampleList &tplRebinedLambdaGlobal =
      m_templateRebined_bf.GetSpectralAxis().GetSamplesVector();

  Int32 kstart, kend;
  kstart = currentRange.GetBegin();
  kend = currentRange.GetEnd();
  Int32 nTpl = kend - kstart + 1;

  TAxisSampleList spcRebinedFluxOverErr2(nSpc);
  for (Int32 j = 0; j < nSpc; j++) {
    spcRebinedFluxOverErr2[j] = spectrumRebinedFluxRaw[j] * inv_err2[j];
  }

  Float64 redshiftValueMeiksin = result->Redshifts[0];

  Log.LogDebug("FitRangez: redshiftValueMeiksin = %f", redshiftValueMeiksin);
  Log.LogDebug("FitRangez: spc[0] = %f", spectrumRebinedLambda[0]);
  Log.LogDebug("FitRangez: spc[max] = %f", spectrumRebinedLambda[nSpc - 1]);
  Log.LogDebug("FitRangez: tpl[0]*zmax = %f",
               tplRebinedLambdaGlobal[kstart] *
                   (1.0 + result->Redshifts[result->Redshifts.size() - 1]));
  Log.LogDebug("FitRangez: tpl[max]*zmin = %f",
               tplRebinedLambdaGlobal[kstart + nTpl - 1] *
                   (1 + result->Redshifts[0]));

  Int32 nshifts = nTpl - nSpc + 1;
  m_nPaddedSamples = nTpl;

  Log.LogDetail("Now fitting using the FFT on "
                "nshifts=%d values, for Meiksin redshift=%f",
                nshifts, redshiftValueMeiksin);

  Log.LogDebug("FitRangez: initializing FFT "
               "with n = %d points",
               m_nPaddedSamples);
  InitFFT(m_nPaddedSamples);

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
  Log.LogDebug("FitRangez: max diff in zgrid=%e", relative_zgrid_error_max);
  if (relative_zgrid_error_max > 5E-7) {
    THROWG(INTERNAL_ERROR, "z_vect and z_vect_verification do not correspond.");
  }

  // check borders
  if (z_vect.size() != z_vect_verif.size())
    THROWG(INTERNAL_ERROR,
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
    Log.LogDebug("FitRangez: IGM disabled, "
                 "min-tpl-lbda=%f",
                 tplRebinedLambdaGlobal[kstart]);
  }

  // prepare best fit data buffer
  TFloat64List bestChi2(nshifts, DBL_MAX);
  TFloat64List bestFitAmp(nshifts, NAN);
  TFloat64List bestFitAmpErr(nshifts, NAN);
  TFloat64List bestFitAmpSigma(nshifts, NAN);
  TFloat64List bestFitDtm(nshifts, NAN);
  TFloat64List bestFitMtm(nshifts, NAN);
  TFloat64List bestISMCoeff(nshifts, NAN);
  TInt32List bestIGMIdx(nshifts, undefIdx);

  // prepare intermediate fit data buffer
  std::vector<std::vector<TFloat64List>> intermediateChi2;
  Int32 nIGMFinal = nIGM;
  if (overrideNIGMTobesaved > nIGM) {
    nIGMFinal = overrideNIGMTobesaved;
  }
  for (Int32 k = 0; k < nshifts; k++) {
    std::vector<TFloat64List> _ChiSquareISMList;
    for (Int32 kism = 0; kism < nISM; kism++) {
      TFloat64List _chi2List(nIGMFinal, DBL_MAX);
      _ChiSquareISMList.push_back(_chi2List);
    }
    intermediateChi2.push_back(_ChiSquareISMList);
  }

  // note that there is no need to copy the ism/igm cause they already exist in
  // the rebinned template
  if (m_enableIGM || m_enableISM) {
    m_templateRebined_bf.InitIsmIgmConfig(kstart, kend, redshiftValueMeiksin);
  }
  for (Int32 kIGM = 0; kIGM < nIGM; kIGM++) {
    if (enableIGM) {
      Log.LogDebug("FitRangez: IGM index=%d", kIGM);
    }

    if (enableIGM) {
      Int32 meiksinIdx = MeiksinList[kIGM];
      m_templateRebined_bf.ApplyMeiksinCoeff(meiksinIdx);
    }

    for (Int32 kISM = 0; kISM < nISM; kISM++) {
      if (m_enableISM) {
        Log.LogDebug("FitRangez: ISM index =%d", kISM);
      }

      if (m_enableISM) {
        Int32 kDust = EbmvList[kISM];
        m_templateRebined_bf.ApplyDustCoeff(kDust);
      }

      const TAxisSampleList &tplRebinedFluxcorr =
          m_templateRebined_bf.GetFluxAxis().GetSamplesVector();
      // extract only the relevant part
      TFloat64List::const_iterator first = tplRebinedFluxcorr.begin() + kstart,
                                   last = tplRebinedFluxcorr.begin() + kend + 1;
      const TAxisSampleList tplRebinedFluxcorr_cropped(first, last);
      TAxisSampleList tpl2RebinedFlux(nTpl);

      if (tplRebinedFluxcorr_cropped.size() != nTpl) {
        THROWG(INTERNAL_ERROR, "vector sizes do not match");
      }
      // compute the square of the corrected flux
      for (Int32 j = 0; j < nTpl; j++) {
        tpl2RebinedFlux[j] =
            tplRebinedFluxcorr_cropped[j] * tplRebinedFluxcorr_cropped[j];
      }

      // Estimate DtM: sumCross
      TFloat64List dtm_vec;
      EstimateXtY(spcRebinedFluxOverErr2, tplRebinedFluxcorr_cropped, nshifts,
                  dtm_vec, 0);

      Int32 dtm_vec_size = dtm_vec.size();

      // Estimate MtM: sumT
      TFloat64List mtm_vec;
      EstimateXtY(inv_err2, tpl2RebinedFlux, nshifts, mtm_vec, 1);

      Log.LogDebug("FitRangez: dtd = %e", dtd);

      // Estimate Chi2
      if (mtm_vec.size() != dtm_vec_size) {
        freeFFTPlans();
        THROWG(INTERNAL_ERROR,
               Formatter() << "xty vector size do not match: dtm size = "
                           << dtm_vec_size << ", mtm size =" << mtm_vec.size());
      }
      // Log.LogDetail(FitRangez: kISM = %d,
      // kIGM = %d", kISM, kIGM);
      TFloat64List chi2(dtm_vec_size, DBL_MAX);
      TFloat64List amp(dtm_vec_size, DBL_MAX);
      TFloat64List amp_sigma(dtm_vec_size);
      TFloat64List amp_err(dtm_vec_size, DBL_MAX);
      for (Int32 k = 0; k < dtm_vec_size; k++) {
        if (mtm_vec[k] == 0.0) {
          amp[k] = 0.0;
          amp_err[k] = 0.0;
          amp_sigma[k] = 0.0;
          chi2[k] = dtd;
        } else {
          amp[k] = dtm_vec[k] / mtm_vec[k];
          amp_err[k] = sqrt(1. / mtm_vec[k]);
          amp_sigma[k] = amp[k] / amp_err[k];
          applyPositiveAndNonNullConstraint(amp_sigma[k], amp[k]);
          chi2[k] =
              dtd - 2 * dtm_vec[k] * amp[k] + mtm_vec[k] * amp[k] * amp[k];
        }
        // Log.LogDetail(FitRangez: chi2[%d] =
        // %f", k, chi2[k]);
      }

      for (Int32 k = 0; k < dtm_vec_size; k++) {
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
          bestFitAmp[k] = amp[k];
          bestFitAmpErr[k] = amp_err[k];
          bestFitAmpSigma[k] = amp_sigma[k];
          bestFitDtm[k] = dtm_vec[k];
          bestFitMtm[k] = mtm_vec[k];
          bestISMCoeff[k] =
              m_enableISM
                  ? m_templateRebined_bf.m_ismCorrectionCalzetti->GetEbmvValue(
                        EbmvList[kISM])
                  : -1.0;
          bestIGMIdx[k] = enableIGM ? MeiksinList[kIGM] : -1;
        }
      }

      Log.LogDebug("FitRangez: spc lbda 0 =%f", spectrumRebinedLambda[0]);
      Log.LogDebug("FitRangez: tpl lbda 0 =%f", tplRebinedLambdaGlobal[kstart]);
      Float64 z_O =
          (spectrumRebinedLambda[0] - tplRebinedLambdaGlobal[kstart]) /
          tplRebinedLambdaGlobal[kstart];
      Log.LogDebug("FitRangez: z 0 =%f", z_O);
    }
  }

  // reversing all vectors
  std::reverse(z_vect.begin(), z_vect.end());
  std::reverse(bestChi2.begin(), bestChi2.end());
  std::reverse(bestFitAmp.begin(), bestFitAmp.end());
  std::reverse(bestFitAmpErr.begin(), bestFitAmpErr.end());
  std::reverse(bestFitAmpSigma.begin(), bestFitAmpSigma.end());
  std::reverse(bestFitDtm.begin(), bestFitDtm.end());
  std::reverse(bestFitMtm.begin(), bestFitMtm.end());
  std::reverse(bestISMCoeff.begin(), bestISMCoeff.end());
  std::reverse(bestIGMIdx.begin(), bestIGMIdx.end());
  TFloat64List intermChi2BufferReversed_array(intermediateChi2.size());
  for (Int32 kism = 0; kism < nISM; kism++) {
    for (Int32 kigm = 0; kigm < nIGMFinal; kigm++) {
      for (Int32 t = 0; t < nshifts; t++)
        intermChi2BufferReversed_array[t] =
            intermediateChi2[nshifts - 1 - t][kism][kigm];
      for (Int32 t = 0; t < nshifts; t++)
        result->ChiSquareIntermediate[t][kism][kigm] =
            intermChi2BufferReversed_array[t];
    }
  }
  for (Int32 k = 0; k < result->Redshifts.size(); k++) {
    result->Overlap[k] = 1.0;
    result->Status[k] = nStatus_OK;
  }
  result->ChiSquare = bestChi2;
  result->FitAmplitude = bestFitAmp;
  result->FitAmplitudeError = bestFitAmpErr;
  result->FitAmplitudeSigma = bestFitAmpSigma;
  result->FitDtM = bestFitDtm;
  result->FitMtM = bestFitMtm;
  result->FitEbmvCoeff = bestISMCoeff;
  result->FitMeiksinIdx = bestIGMIdx;

  freeFFTPlans();
  return 0;
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
  return FindTplSpectralIndex(m_ssSpectrum.GetSpectralAxis(),
                              m_templateRebined_bf.GetSpectralAxis(),
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
    THROWG(INTERNAL_ERROR, " Failed to find indexes");

  if (ilbdamin > ilbdamax) {
    THROWG(INTERNAL_ERROR, Formatter()
                               << "Problem with "
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
 *spectrum grid) 2. input: if additional_spcMasks size is 0, no additional mask
 *will be used, otherwise its size should match the redshifts list size
 *
 * lambdaRange is not clamped
 **/
std::shared_ptr<COperatorResult> COperatorTemplateFittingLog::Compute(
    const std::shared_ptr<const CTemplate> &logSampledTpl,
    Float64 overlapThreshold, const std::vector<CMask> &additional_spcMasks,
    std::string opt_interp, bool opt_extinction, bool opt_dustFitting,
    Float64 opt_continuum_null_amp_threshold,
    const CPriorHelper::TPriorZEList &logpriorze, Int32 FitEbmvIdx,
    Int32 FitMeiksinIdx) {
  Log.LogDetail("starting computation for template: %s",
                logSampledTpl->GetName().c_str());

  if (opt_dustFitting && logSampledTpl->CalzettiInitFailed())
    THROWG(INTERNAL_ERROR, "ISM is no initialized");

  if (opt_dustFitting &&
      FitEbmvIdx >=
          logSampledTpl->m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs()) {
    THROWG(INTERNAL_ERROR,
           Formatter() << " Invalid calzetti index: (FitEbmvIdx=" << FitEbmvIdx
                       << ",while NPrecomputedEbmvCoeffs="
                       << logSampledTpl->m_ismCorrectionCalzetti
                              ->GetNPrecomputedEbmvCoeffs()
                       << ")");
  }

  if (opt_extinction && logSampledTpl->MeiksinInitFailed()) {
    THROWG(INTERNAL_ERROR, "IGM is not initialized");
  }

  if (!logSampledTpl->GetSpectralAxis().IsLogSampled()) {
    THROWG(INTERNAL_ERROR, "template is not log sampled");
  }
  // check if spc and tpl have same step
  const Float64 epsilon = 1E-8;
  if (std::abs(m_logSampledSpectrum.GetSpectralAxis().GetlogGridStep() -
               logSampledTpl->GetSpectralAxis().GetlogGridStep()) > epsilon)
    THROWG(INTERNAL_ERROR, "tpl and spc are not sampled with the same step");

  m_continuum_null_amp_threshold = opt_continuum_null_amp_threshold;

  // subsample template if necessary
  if (m_ssRatio == 1) { // no required subsampling
    m_templateRebined_bf = *logSampledTpl;
  } else {
    TInt32Range ilbda = FindTplSpectralIndex(m_ssSpectrum.GetSpectralAxis(),
                                             logSampledTpl->GetSpectralAxis(),
                                             TFloat64Range(m_redshifts));
    TFloat64List mask_tpl =
        logSampledTpl->GetSpectralAxis().GetSubSamplingMask(m_ssRatio, ilbda);

    m_templateRebined_bf = CTemplate(*logSampledTpl, mask_tpl);
    // double make sure that subsampled spectrum is well sampled
    if (!m_templateRebined_bf.GetSpectralAxis().IsLogSampled(m_logstep)) {
      THROWG(INTERNAL_ERROR, "subsampled template "
                             "is not log sampled with the redshift step");
    }
  }

  //**************** Fitting at all redshifts ****************//
  // Note: below corresponds to ::BasicFit code except that redshift loop
  // belongs to ::compute
  // Optionally apply some IGM absorption
  TInt32List MeiksinList;
  TInt32List EbmvList;
  m_templateRebined_bf.GetIsmIgmIdxList(opt_extinction, opt_dustFitting,
                                        MeiksinList, EbmvList, FitEbmvIdx,
                                        FitMeiksinIdx);
  Int32 nIGMCoeffs = MeiksinList.size();
  Int32 nISMCoeffs = EbmvList.size();

  m_enableIGM = opt_extinction;
  m_enableISM = opt_dustFitting;
  std::shared_ptr<CTemplateFittingResult> result =
      std::make_shared<CTemplateFittingResult>(m_redshifts.size(), nISMCoeffs,
                                               nIGMCoeffs);
  result->Redshifts = m_redshifts;

  // WARNING: no additional masks coded for use as of 2017-06-13
  if (additional_spcMasks.size() != 0)
    Flag.warning(WarningCode::TEMPLATEFITTINGLOG_NO_MASK,
                 Formatter()
                     << "COperatorTemplateFittingLog::" << __func__
                     << "No additional masks used. "
                        "Feature not coded for this log-lambda operator");

  if (logpriorze.size() > 0 && logpriorze.size() != m_redshifts.size()) {
    THROWG(INTERNAL_ERROR,
           Formatter()
               << "Vector size do not match between logpriorz and redshift: "
               << logpriorze.size() << " != " << m_redshifts.size());
  }

  Int32 retFit = FitAllz(result, MeiksinList, EbmvList, logpriorze);

  if (retFit != 0) {
    THROWG(INTERNAL_ERROR, Formatter()
                               << "FitAllz failed with error " << retFit);
  }

  //**************** End Fitting at all redshifts ****************//

  // overlap warning
  Float64 overlapValidInfZ = -1;
  for (Int32 i = 0; i < m_redshifts.size(); i++) {
    if (result->Overlap[i] >= overlapThreshold && overlapValidInfZ == -1) {
      overlapValidInfZ = m_redshifts[i];
      break;
    }
  }
  Float64 overlapValidSupZ = -1;
  for (Int32 i = m_redshifts.size() - 1; i >= 0; i--) {
    if (result->Overlap[i] >= overlapThreshold && overlapValidSupZ == -1) {
      overlapValidSupZ = m_redshifts[i];
      break;
    }
  }
  if (overlapValidInfZ != m_redshifts.front() ||
      overlapValidSupZ != m_redshifts.back()) {
    Log.LogInfo("overlap warning for %s: "
                "minz=%.3f, maxz=%.3f",
                logSampledTpl->GetName().c_str(), overlapValidInfZ,
                overlapValidSupZ);
  }

  // estimate CstLog for PDF estimation
  result->CstLog = EstimateLikelihoodCstLog(
      m_ssSpectrum, m_ssLambdaRange); // 0.0;//Todo: check how to estimate
                                      // that value for loglambda//

  return result;
}
