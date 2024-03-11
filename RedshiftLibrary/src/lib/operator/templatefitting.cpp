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

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/operator/templatefitting.h"
#include "RedshiftLibrary/operator/templatefittingresult.h"
#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"

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
  bool amplForcePositive = true;
  bool chisquareSetAtLeastOnce = false;

  Int32 EbmvListSize = EbmvList.size();
  Int32 MeiksinListSize = MeiksinList.size();
  m_option_igmFastProcessing = (MeiksinListSize > 1 ? true : false);

  bool apply_priore =
      !logpriore.empty() && !tpl->CalzettiInitFailed() &&
      (logpriore.size() ==
       tpl->m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs());

  Int32 iEbmvCoeffMin = EbmvList.front();
  Int32 iEbmvCoeffMax = EbmvList.back();

  TFittingIsmIgmResult result(EbmvListSize, MeiksinListSize, m_spectra.size());

  TFloat64RangeList currentRanges(m_spectra.size()); // restframe
  for (Int32 spcIndex = 0; spcIndex < m_spectra.size(); spcIndex++) {

    RebinTemplate(tpl, redshift, currentRanges[spcIndex],
                  result.overlapFraction[spcIndex], overlapThreshold, spcIndex);
    currentRanges[spcIndex].getClosedIntervalIndices(
        m_templateRebined_bf[spcIndex].GetSpectralAxis().GetSamplesVector(),
        m_kStart[spcIndex], m_kEnd[spcIndex]);
  }
  if (opt_dustFitting || opt_extinction)
    InitIsmIgmConfig(redshift, tpl->m_ismCorrectionCalzetti,
                     tpl->m_igmCorrectionMeiksin, EbmvListSize);

  CPriorHelper::SPriorTZE logpriorTZEempty = {};

  // Loop on the meiksin Idx
  bool igmLoopUseless_WavelengthRange = false;
  for (Int32 kM = 0; kM < MeiksinListSize; kM++) {
    if (igmLoopUseless_WavelengthRange) {
      // Now copy from the already calculated k>0 igm values
      for (Int32 kism = 0; kism < result.ChiSquareInterm.size(); kism++) {
        for (Int32 kigm = 1; kigm < result.ChiSquareInterm[kism].size();
             kigm++) {
          result.ChiSquareInterm[kism][kigm] = result.ChiSquareInterm[kism][0];
          result.IsmCalzettiCoeffInterm[kism][kigm] =
              result.IsmCalzettiCoeffInterm[kism][0];
          result.IgmMeiksinIdxInterm[kism][kigm] =
              result.IgmMeiksinIdxInterm[kism][0];
        }
      }
      break;
    }
    Int32 meiksinIdx = MeiksinList[kM]; // index for the Meiksin curve (0-6; 3
    // being the median extinction value)

    bool igmCorrectionAppliedOnce = false;
    // Meiksin IGM extinction
    if (opt_extinction) {
      for (Int32 spcIndex = 0; spcIndex < m_spectra.size(); spcIndex++) {

        // check if lya belongs to current range.
        if (CheckLyaIsInCurrentRange(currentRanges[spcIndex]))
          igmCorrectionAppliedOnce = false;
        else
          igmCorrectionAppliedOnce = ApplyMeiksinCoeff(meiksinIdx, spcIndex);
      }
      if (!igmCorrectionAppliedOnce)
        igmLoopUseless_WavelengthRange = true;
    }
    // Loop on the EBMV dust coeff
    for (Int32 kEbmv_ = 0; kEbmv_ < EbmvListSize; kEbmv_++) {
      Int32 kEbmv = EbmvList[kEbmv_];
      Float64 coeffEBMV = -1.0; // no ism by default

      if (opt_dustFitting) {
        coeffEBMV = tpl->m_ismCorrectionCalzetti->GetEbmvValue(kEbmv);
        for (Int32 spcIndex = 0; spcIndex < m_spectra.size(); spcIndex++)
          ApplyDustCoeff(kEbmv, spcIndex);
      }

      const CPriorHelper::SPriorTZE &logpriorTZE =
          apply_priore ? logpriore[kEbmv] : logpriorTZEempty;
      TFittingResult fitRes;
      for (Int32 spcIndex = 0; spcIndex < m_spectra.size(); spcIndex++) {
        fitRes.cross_result +=
            ComputeCrossProducts(kM, kEbmv_, redshift, spcIndex);
      }
      ComputeAmplitudeAndChi2(fitRes, logpriorTZE);
      Log.LogDebug(Formatter() << "BasicFit: z=" << redshift << " fit="
                               << fitRes.chiSquare << " coeffEBMV=" << coeffEBMV
                               << " meiksinIdx=" << meiksinIdx);

      result.ChiSquareInterm[kEbmv_][kM] = fitRes.chiSquare;
      result.IsmCalzettiCoeffInterm[kEbmv_][kM] = coeffEBMV;
      result.IgmMeiksinIdxInterm[kEbmv_][kM] = meiksinIdx;

      if (fitRes.chiSquare < result.chiSquare) {
        TFittingResult &result_base = result; // upcasting (slicing)
        result_base = fitRes;                 // slicing, preserving specific
        // TFittingIsmIGmResult members

        result.EbmvCoeff = coeffEBMV;
        result.MeiksinIdx =
            igmCorrectionAppliedOnce == true ? meiksinIdx : undefIdx;
        chisquareSetAtLeastOnce = true;
      }
    }
  }

  if (!chisquareSetAtLeastOnce) {
    THROWG(INVALID_MERIT_VALUES,
           Formatter() << "Template " << tpl->GetName().c_str()
                       << ": Not even one single valid fit/merit value found");
  }

  return result;
}

void COperatorTemplateFitting::InitIsmIgmConfig(
    Float64 redshift,
    const std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
        &ismCorrectionCalzetti,
    const std::shared_ptr<const CSpectrumFluxCorrectionMeiksin>
        &igmCorrectionMeiksin,
    Int32 EbmvListSize) {

  for (Int32 spcIndex = 0; spcIndex < m_spectra.size(); spcIndex++)
    m_templateRebined_bf[spcIndex].InitIsmIgmConfig(
        m_kStart[spcIndex], m_kEnd[spcIndex], redshift, ismCorrectionCalzetti,
        igmCorrectionMeiksin);

  m_sumCross_outsideIGM.assign(m_spectra.size(),
                               TFloat64List(EbmvListSize, 0.0));
  m_sumT_outsideIGM.assign(m_spectra.size(), TFloat64List(EbmvListSize, 0.0));
  m_sumS_outsideIGM.assign(m_spectra.size(), TFloat64List(EbmvListSize, 0.0));
}

TCrossProductResult COperatorTemplateFitting::ComputeCrossProducts(
    Int32 kM, Int32 kEbmv_, Float64 redshift, Int32 spcIndex) {
  const CSpectrumFluxAxis &spcFluxAxis = m_spectra[spcIndex]->GetFluxAxis();
  const TAxisSampleList &Yspc = spcFluxAxis.GetSamplesVector();
  const TAxisSampleList &Ytpl =
      m_templateRebined_bf[spcIndex].GetFluxAxis().GetSamplesVector();
  const TAxisSampleList &Xtpl =
      m_templateRebined_bf[spcIndex].GetSpectralAxis().GetSamplesVector();

  const CMask &spcMaskAdditional =
      m_maskBuilder->getMask(m_spectra[spcIndex]->GetSpectralAxis(),
                             *m_lambdaRanges[spcIndex], redshift);
  TCrossProductResult fitResult;

  Float64 &sumCross = fitResult.sumCross;
  Float64 &sumT = fitResult.sumT;
  Float64 &sumS = fitResult.sumS;

  Float64 sumCross_IGM = 0.0;
  Float64 sumT_IGM = 0.0;
  Float64 sumS_IGM = 0.0;
  Int32 sumsIgmSaved = 0;

  Float64 err2 = 0.0;
  Int32 numDevs = 0;
  Int32 numDevsFull = 0;
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

    numDevsFull++;

    if (spcMaskAdditional[j]) {

      numDevs++;
      err2 = 1.0 / (error[j] * error[j]);

      // Tonry&Davis formulation
      sumCross += Yspc[j] * Ytpl[j] * err2;
      sumT += Ytpl[j] * Ytpl[j] * err2;
      sumS += Yspc[j] * Yspc[j] * err2;

      if (std::isinf(err2) || std::isnan(err2)) {
        THROWG(INTERNAL_ERROR,
               Formatter() << "found invalid inverse variance : err2=" << err2
                           << ", for index=" << j
                           << " at restframe wl=" << Xtpl[j]);
      }

      if (std::isinf(sumS) || std::isnan(sumS) || sumS != sumS)
        THROWG(INTERNAL_ERROR,
               Formatter() << "Invalid dtd value: dtd=" << sumS
                           << ", Yspc=" << Yspc[j] << ", err2=" << err2
                           << ", error=" << error[j] << ", for index=" << j
                           << " at restframe wl=" << Xtpl[j]);

      if (std::isinf(sumT) || std::isnan(sumT)) {
        THROWG(INTERNAL_ERROR, Formatter() << "Invalid mtm value:" << sumT
                                           << " for index=" << j
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

  if (numDevs == 0) {
    THROWG(INTERNAL_ERROR, "empty leastsquare sum");
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

/**
 * \brief
 *
 * input: if additional_spcMasks size is 0, no additional mask will be used,
 *otherwise its size should match the redshifts list size
 *
 **/
std::shared_ptr<COperatorResult> COperatorTemplateFitting::Compute(
    const std::shared_ptr<const CTemplate> &tpl, Float64 overlapThreshold,
    std::string opt_interp, bool opt_extinction, bool opt_dustFitting,
    Float64 opt_continuum_null_amp_threshold,
    const CPriorHelper::TPriorZEList &logpriorze, Int32 FitEbmvIdx,
    Int32 FitMeiksinIdx) {
  Log.LogDetail(
      Formatter()
      << "  Operator-TemplateFitting: starting computation for template: "
      << tpl->GetName().c_str());

  if (opt_dustFitting && tpl->CalzettiInitFailed())
    THROWG(INTERNAL_ERROR, "ISM is not initialized");

  if (opt_dustFitting &&
      FitEbmvIdx >= tpl->m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs())
    THROWG(
        INTERNAL_ERROR,
        Formatter() << "Invalid calzetti index. (FitEbmvIdx=" << FitEbmvIdx
                    << ", while NPrecomputedEbmvCoeffs="
                    << tpl->m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs()
                    << ")");

  if (opt_extinction && tpl->MeiksinInitFailed()) {
    THROWG(INTERNAL_ERROR, "IGM is not initialized");
  }
  m_continuum_null_amp_threshold = opt_continuum_null_amp_threshold;

  std::shared_ptr<CTemplateFittingResult> result =
      std::make_shared<CTemplateFittingResult>(m_redshifts.size());
  TInt32List MeiksinList;
  TInt32List EbmvList;
  tpl->GetIsmIgmIdxList(opt_extinction, opt_dustFitting, MeiksinList, EbmvList,
                        FitEbmvIdx, FitMeiksinIdx);

  result->Redshifts = m_redshifts;

  if (logpriorze.size() > 0 && logpriorze.size() != m_redshifts.size())
    THROWG(INTERNAL_ERROR,
           Formatter() << "prior list size(" << logpriorze.size()
                       << ") does not match the input redshift-list size :"
                       << m_redshifts.size());

  for (Int32 i = 0; i < m_redshifts.size(); i++) {
    const CPriorHelper::TPriorEList &logp =
        logpriorze.size() > 0 && logpriorze.size() == m_redshifts.size()
            ? logpriorze[i]
            : CPriorHelper::TPriorEList();

    Float64 redshift = result->Redshifts[i];

    TFittingIsmIgmResult result_z =
        BasicFit(tpl, redshift, overlapThreshold, opt_extinction,
                 opt_dustFitting, logp, MeiksinList, EbmvList);

    result->set_at_redshift(i, std::move(result_z));
  }

  // overlap warning
  Float64 overlapValidInfZ = -1;
  for (Int32 i = 0; i < m_redshifts.size(); i++) {
    bool ok = true;
    for (auto ov : result->Overlap[i])
      ok &= (ov >= overlapThreshold);
    if (ok) {
      overlapValidInfZ = m_redshifts[i];
      break;
    }
  }
  Float64 overlapValidSupZ = -1;
  for (Int32 i = m_redshifts.size() - 1; i >= 0; i--) {
    bool ok = true;
    for (auto ov : result->Overlap[i])
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
                << tpl->GetName().c_str()
                << ": "
                   "minz="
                << overlapValidInfZ << ", maxz=" << overlapValidSupZ);
  }

  // estimate CstLog for PDF estimation
  result->CstLog = EstimateLikelihoodCstLog();

  return result;
}
