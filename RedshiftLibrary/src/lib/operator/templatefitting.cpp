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
#include <algorithm> // std::sort
#include <climits>
#include <cmath>
#include <iostream>
#include <sstream>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/numeric/conversion/bounds.hpp>

#include "RedshiftLibrary/common/size.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/operator/templatefitting.h"
#include "RedshiftLibrary/operator/templatefittingresult.h"
#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/statistics/fitquality.h"

using namespace NSEpic;
using namespace std;

/**
 * @brief COperatorTemplateFitting::BasicFit
 * @param tpl
 * @param redshift
 * @param overlapThreshold
 * @param forcedAmplitude
 * @param opt_extinction
 * @param opt_dustFitting
 * @param spcMaskAdditional
 * @param priorjoint_pISM_tpl_z : vector size = nISM, joint prior p(ISM, TPL, Z)
 */
TFittingIsmIgmResult COperatorTemplateFitting::BasicFit(
    const std::shared_ptr<const CTemplate> &tpl, Float64 redshift,
    Float64 overlapThreshold, bool opt_extinction, bool opt_dustFitting,
    const CPriorHelper::TPriorEList &logpriore, const TInt32List &MeiksinList,
    const TInt32List &EbmvList) {
  bool chisquareSetAtLeastOnce = false;

  Int32 EbmvListSize = EbmvList.size();
  Int32 MeiksinListSize = MeiksinList.size();

  bool apply_priore =
      !logpriore.empty() && !tpl->CalzettiInitFailed() &&
      (ssize(logpriore) ==
       tpl->m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs());

  TFittingIsmIgmResult result(EbmvListSize, MeiksinListSize, m_spectra.size());

  TFloat64RangeList currentRanges(m_spectra.size()); // restframe
  for (Int32 spcIndex = 0; spcIndex < ssize(m_spectra);
       spcIndex++) { // loop for multi obs
    // Adapts template bin
    RebinTemplate(tpl, redshift, currentRanges[spcIndex],
                  result.overlapFraction[spcIndex], overlapThreshold, spcIndex);
    // Sets template indexes between which data is of interest for this given
    // spectrum / redshift
    currentRanges[spcIndex].getClosedIntervalIndices(
        m_templateRebined_bf[spcIndex].GetSpectralAxis().GetSamplesVector(),
        m_kStart[spcIndex], m_kEnd[spcIndex]);
  }

  // get masks & determine number of samples actually used
  auto const &[mask_list, n_samples] = getMaskListAndNSamples(redshift);

  if (opt_extinction)
    opt_extinction = igmIsInRange(currentRanges);
  m_option_igmFastProcessing =
      (opt_extinction && (MeiksinListSize > 1)) ? true : false;

  if (opt_dustFitting || opt_extinction) {
    for (Int32 spcIndex = 0; spcIndex < ssize(m_spectra); spcIndex++) {
      InitIsmIgmConfig(redshift, m_kStart[spcIndex], m_kEnd[spcIndex],
                       spcIndex);
    }
  }

  if (m_option_igmFastProcessing)
    init_fast_igm_processing(EbmvListSize);

  CPriorHelper::SPriorTZE logpriorTZEempty = {};

  // Loop on the meiksin Idx
  bool igmNotapplicable = true;
  for (Int32 kM = 0; kM < MeiksinListSize; kM++) {
    if (kM > 0 && igmNotapplicable) {
      // Now copy from the already calculated k>0 igm values
      for (Int32 kigm = 1; kigm < ssize(result.IgmMeiksinIdxInterm); kigm++)
        result.IgmMeiksinIdxInterm[kigm] = result.IgmMeiksinIdxInterm[0];
      for (Int32 kism = 0; kism < ssize(result.ChiSquareInterm); kism++)
        for (Int32 kigm = 1; kigm < ssize(result.ChiSquareInterm[kism]); kigm++)
          result.ChiSquareInterm[kism][kigm] = result.ChiSquareInterm[kism][0];
      break;
    }
    Int32 meiksinIdx = undefIdx;

    // Meiksin IGM extinction
    if (opt_extinction) {
      meiksinIdx = MeiksinList[kM]; // index for the Meiksin curve (0-6; 3
                                    // being the median extinction value)
      for (Int32 spcIndex = 0; spcIndex < ssize(m_spectra); spcIndex++)
        if (ApplyMeiksinCoeff(meiksinIdx, spcIndex))
          igmNotapplicable = false;
      if (igmNotapplicable)
        meiksinIdx = undefIdx;
    }

    result.IgmMeiksinIdxInterm[kM] = meiksinIdx;

    // Loop on the EBMV dust coeff
    for (Int32 kEbmv_ = 0; kEbmv_ < EbmvListSize; kEbmv_++) {
      Int32 kEbmv = EbmvList[kEbmv_];
      Float64 coeffEBMV = -1.0; // no ism by default

      if (opt_dustFitting) {
        coeffEBMV = tpl->m_ismCorrectionCalzetti->GetEbmvValue(kEbmv);
        for (Int32 spcIndex = 0; spcIndex < ssize(m_spectra); spcIndex++)
          ApplyDustCoeff(kEbmv, spcIndex);
      }
      result.IsmCalzettiIdxInterm[kEbmv_] = kEbmv;

      // Priors
      const CPriorHelper::SPriorTZE &logpriorTZE =
          apply_priore ? logpriore[kEbmv] : logpriorTZEempty;

      TFittingResult fitRes;

      // Chi2 calculation
      for (Int32 spcIndex = 0; spcIndex < ssize(m_spectra); spcIndex++) {
        fitRes.cross_result += ComputeCrossProducts(
            kM, kEbmv_, redshift, mask_list[spcIndex], spcIndex);
      }
      ComputeAmplitudeAndChi2(fitRes, logpriorTZE);

      Log.LogDebug(Formatter() << "BasicFit: z=" << redshift << " fit="
                               << fitRes.chiSquare << " coeffEBMV=" << coeffEBMV
                               << " meiksinIdx=" << meiksinIdx);

      result.ChiSquareInterm[kEbmv_][kM] = fitRes.chiSquare;

      // if minimizes chi2, sets it as the result
      if (fitRes.chiSquare < result.chiSquare) {
        TFittingResult &result_base =
            result; // upcasting (slicing) slicing, preserving specific
                    // TFittingIsmIGmResult members
        result_base = fitRes;
        result.reducedChiSquare =
            NSFitQuality::reducedChi2(result.chiSquare, n_samples);
        result.pValue = NSFitQuality::pValue(result.chiSquare, n_samples);
        result.ebmvCoef = coeffEBMV;
        result.meiksinIdx = meiksinIdx;
        chisquareSetAtLeastOnce = true;
      }
    }
  }

  if (!chisquareSetAtLeastOnce) {
    THROWG(ErrorCode::INVALID_MERIT_VALUES,
           Formatter() << "Template " << tpl->GetName()
                       << ": Not even one single valid fit/merit value found");
  }

  return result;
}

std::pair<TList<CMask>, Int32>
COperatorTemplateFitting::getMaskListAndNSamples(Float64 redshift) const {
  // get masks & determine number of samples actually used
  TList<CMask> mask_list;
  Int32 n_samples = 0; // total number of samples
  mask_list.reserve(m_spectra.size());
  for (Int32 spcIndex = 0; spcIndex < ssize(m_spectra); spcIndex++) {
    const CMask &mask =
        m_maskBuilder->getMask(m_spectra[spcIndex]->GetSpectralAxis(),
                               *m_lambdaRanges[spcIndex], redshift, spcIndex);
    n_samples +=
        std::count(mask.getMaskList().begin() + m_kStart[spcIndex],
                   mask.getMaskList().begin() + m_kEnd[spcIndex] + 1, Mask(1));
    mask_list.push_back(std::move(mask));
  }

  return std::make_pair(std::move(mask_list), n_samples);
}

void COperatorTemplateFitting::init_fast_igm_processing(Int32 EbmvListSize) {

  m_sumCross_outsideIGM.assign(m_spectra.size(),
                               TFloat64List(EbmvListSize, 0.0));
  m_sumT_outsideIGM.assign(m_spectra.size(), TFloat64List(EbmvListSize, 0.0));
  m_sumS_outsideIGM.assign(m_spectra.size(), TFloat64List(EbmvListSize, 0.0));
}

bool COperatorTemplateFitting::igmIsInRange(
    const TFloat64RangeList &ranges) const {
  for (auto const &range : ranges) {
    if (range.GetBegin() <= RESTLAMBDA_LYA)
      return true;
  }
  return false;
}

TCrossProductResult COperatorTemplateFitting::ComputeCrossProducts(
    Int32 kM, Int32 kEbmv_, Float64 redshift, CMask const &mask,
    Int32 spcIndex) {
  const CSpectrumFluxAxis &spcFluxAxis = m_spectra[spcIndex]->GetFluxAxis();
  const TAxisSampleList &Yspc = spcFluxAxis.GetSamplesVector();
  const TAxisSampleList &Ytpl =
      m_templateRebined_bf[spcIndex].GetFluxAxis().GetSamplesVector();
  const TAxisSampleList &Xtpl =
      m_templateRebined_bf[spcIndex].GetSpectralAxis().GetSamplesVector();

  TCrossProductResult fitResult;

  Float64 &sumCross = fitResult.sumCross;
  Float64 &sumT = fitResult.sumT;
  Float64 &sumS = fitResult.sumS;

  Float64 sumCross_IGM = 0.0;
  Float64 sumT_IGM = 0.0;
  Float64 sumS_IGM = 0.0;
  Int32 sumsIgmSaved = 0;

  Float64 err2 = 0.0;
  const CSpectrumNoiseAxis &error = spcFluxAxis.GetError();

  Int32 kIgmEnd = m_option_igmFastProcessing
                      ? m_templateRebined_bf[spcIndex].GetIgmEndIndex()
                      : -1;
  Int32 kEndloop =
      m_option_igmFastProcessing && kM > 0 ? kIgmEnd : m_kEnd[spcIndex];
  for (Int32 j = m_kStart[spcIndex]; j <= kEndloop; j++) {
    if (m_option_igmFastProcessing && sumsIgmSaved == 0 && j > kIgmEnd) {
      // store intermediate sums for IGM range
      sumCross_IGM = sumCross;
      sumT_IGM = sumT;
      sumS_IGM = sumS;
      sumsIgmSaved = 1;
    }

    if (mask[j]) {

      err2 = 1.0 / (error[j] * error[j]);

      // Tonry&Davis formulation
      sumCross += Yspc[j] * Ytpl[j] * err2;
      sumT += Ytpl[j] * Ytpl[j] * err2;
      sumS += Yspc[j] * Yspc[j] * err2;

      if (std::isinf(err2) || std::isnan(err2)) {
        THROWG(ErrorCode::INTERNAL_ERROR,
               Formatter() << "found invalid inverse variance : err2=" << err2
                           << ", for index=" << j
                           << " at restframe wl=" << Xtpl[j]);
      }

      if (std::isinf(sumS) || std::isnan(sumS) || sumS != sumS)
        THROWG(ErrorCode::INTERNAL_ERROR,
               Formatter() << "Invalid dtd value: dtd=" << sumS
                           << ", Yspc=" << Yspc[j] << ", err2=" << err2
                           << ", error=" << error[j] << ", for index=" << j
                           << " at restframe wl=" << Xtpl[j]);

      if (std::isinf(sumT) || std::isnan(sumT)) {
        THROWG(ErrorCode::INTERNAL_ERROR,
               Formatter() << "Invalid mtm value:" << sumT << " for index=" << j
                           << " at restframe wl=" << Xtpl[j]);
      }
    }
  }

  if (m_option_igmFastProcessing) {
    if (kM == 0) {
      m_sumCross_outsideIGM[spcIndex][kEbmv_] = sumCross - sumCross_IGM;
      m_sumT_outsideIGM[spcIndex][kEbmv_] = sumT - sumT_IGM;
      m_sumS_outsideIGM[spcIndex][kEbmv_] = sumS - sumS_IGM;
    } else {
      sumCross += m_sumCross_outsideIGM[spcIndex][kEbmv_];
      sumT += m_sumT_outsideIGM[spcIndex][kEbmv_];
      sumS += m_sumS_outsideIGM[spcIndex][kEbmv_];
    }
  }

  if (sumT == 0.0) {
    THROWG(ErrorCode::INTERNAL_ERROR, "empty leastsquare sum");
  }

  return fitResult;
}

void COperatorTemplateFitting::ComputeAmplitudeAndChi2(
    TFittingResult &fitResult,
    const CPriorHelper::SPriorTZE &logpriorTZE) const {
  const Float64 &sumCross = fitResult.cross_result.sumCross;
  const Float64 &sumT = fitResult.cross_result.sumT;
  const Float64 &sumS = fitResult.cross_result.sumS;

  Float64 &ampl = fitResult.ampl;
  Float64 &ampl_err = fitResult.ampl_err;
  Float64 &ampl_sigma = fitResult.ampl_sigma;
  Float64 &fit = fitResult.chiSquare;

  bool apply_priore =
      (logpriorTZE.A_sigma > 0.0 && logpriorTZE.betaA > 0.0) ? true : false;

  if (sumT == 0) {
    ampl = 0.0;
    ampl_err = 0.0;
    fit = sumS;
    ampl_sigma = 0.0;
  } else {

    ampl = sumCross / sumT;
    ampl_err = sqrt(1. / sumT);

    if (apply_priore) {
      Float64 bss2 =
          logpriorTZE.betaA / (logpriorTZE.A_sigma * logpriorTZE.A_sigma);
      ampl = (sumCross + logpriorTZE.A_mean * bss2) / (sumT + bss2);
      ampl_err = sqrt(sumT) / (sumT + bss2);
    }

    ampl_sigma = ampl / ampl_err;

    applyPositiveAndNonNullConstraint(ampl_sigma, ampl);

    // Generalized method (ampl can take any value now) for chi2 estimate
    fit = sumS + sumT * ampl * ampl - 2. * ampl * sumCross;
  }

  Float64 &logprior = fitResult.logprior;
  if (apply_priore) {
    logprior += -2. * logpriorTZE.betaTE * logpriorTZE.logprior_precompTE;
    logprior += -2. * logpriorTZE.betaA * logpriorTZE.logprior_precompA;
    logprior += -2. * logpriorTZE.betaZ * logpriorTZE.logprior_precompZ;
    if (logpriorTZE.A_sigma > 0.0) {
      Float64 logPa = logpriorTZE.betaA * (ampl - logpriorTZE.A_mean) *
                      (ampl - logpriorTZE.A_mean) /
                      (logpriorTZE.A_sigma * logpriorTZE.A_sigma);
      logprior += logPa;
    }
    fit += logprior;
  }
}

std::shared_ptr<CTemplateFittingResult> COperatorTemplateFitting::Compute(
    const std::shared_ptr<const CTemplate> &tpl, Float64 overlapThreshold,
    std::string opt_interp, bool opt_extinction, bool opt_dustFitting,
    Float64 opt_continuum_null_amp_threshold,
    const CPriorHelper::TPriorZEList &logprior, Int32 FitEbmvIdx,
    Int32 FitMeiksinIdx, TInt32Range zIdxRangeToCompute,
    std::shared_ptr<CTemplateFittingResult> const &result) {
  Log.LogDetail(
      Formatter()
      << "  Operator-TemplateFitting: starting computation for template: "
      << tpl->GetName());

  if (opt_dustFitting && tpl->CalzettiInitFailed())
    THROWG(ErrorCode::INTERNAL_ERROR, "ISM is not initialized");

  if (opt_dustFitting &&
      FitEbmvIdx >= tpl->m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs())
    THROWG(
        ErrorCode::INTERNAL_ERROR,
        Formatter() << "Invalid calzetti index. (FitEbmvIdx=" << FitEbmvIdx
                    << ", while NPrecomputedEbmvCoeffs="
                    << tpl->m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs()
                    << ")");

  if (opt_extinction && tpl->MeiksinInitFailed()) {
    THROWG(ErrorCode::INTERNAL_ERROR, "IGM is not initialized");
  }
  m_continuum_null_amp_threshold = opt_continuum_null_amp_threshold;
  TIgmIsmIdxs igmIsmIdxs = tpl->GetIsmIgmIdxList(
      opt_extinction, opt_dustFitting, FitEbmvIdx, FitMeiksinIdx);

  std::shared_ptr<CTemplateFittingResult> const templateFittingResult =
      result == nullptr
          ? std::make_shared<CTemplateFittingResult>(m_redshifts.size())
          : result;

  templateFittingResult->Redshifts = m_redshifts;

  if (logprior.size() > 0 && logprior.size() != m_redshifts.size())
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "prior list size(" << logprior.size()
                       << ") does not match the input redshift-list size :"
                       << m_redshifts.size());
  bool hasLogPrior = !logprior.empty();

  TInt32List zIdxs;
  if (zIdxRangeToCompute == TInt32Range(undefIdx, undefIdx))
    zIdxRangeToCompute = TInt32Range(0, m_redshifts.size() - 1);
  for (auto zIdx : zIdxRangeToCompute) {
    Float64 redshift = templateFittingResult->Redshifts[zIdx];
    // TODO move a condition up loop
    const CPriorHelper::TPriorEList &logp =
        hasLogPrior ? logprior[zIdx] : CPriorHelper::TPriorEList();

    TFittingIsmIgmResult result_z =
        BasicFit(tpl, redshift, overlapThreshold, opt_extinction,
                 opt_dustFitting, logp, igmIsmIdxs.igmIdxs, igmIsmIdxs.ismIdxs);

    templateFittingResult->set_at_redshift(zIdx, std::move(result_z));
  }

  // Question what to do with overlap in second pass ?
  // overlap warning
  Float64 overlapValidInfZ = -1;
  for (Int32 i = 0; i < ssize(m_redshifts); i++) {
    bool ok = true;
    for (auto ov : templateFittingResult->Overlap[i])
      ok &= (ov >= overlapThreshold);
    if (ok) {
      overlapValidInfZ = m_redshifts[i];
      break;
    }
  }
  Float64 overlapValidSupZ = -1;
  for (Int32 i = m_redshifts.size() - 1; i >= 0; i--) {
    bool ok = true;
    for (auto ov : templateFittingResult->Overlap[i])
      ok &= (ov >= overlapThreshold);
    if (ok) {
      overlapValidSupZ = m_redshifts[i];
      break;
    }
  }
  if (overlapValidInfZ != m_redshifts[0] ||
      overlapValidSupZ != m_redshifts[m_redshifts.size() - 1]) {
    Log.LogInfo(Formatter()
                << "  Operator-TemplateFitting: overlap warning for "
                << tpl->GetName()
                << ": "
                   "minz="
                << overlapValidInfZ << ", maxz=" << overlapValidSupZ);
  }

  // estimate CstLog for PDF estimation
  templateFittingResult->CstLog = EstimateLikelihoodCstLog();

  return templateFittingResult;
}
