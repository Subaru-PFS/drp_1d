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
#include <numeric>
#include <sstream>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/chrono/thread_clock.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>

#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/operator/tplcombination.h"
#include "RedshiftLibrary/operator/tplcombinationresult.h"
#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"

namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;

void COperatorTplcombination::BasicFit_preallocateBuffers(
    const CSpectrum &spectrum, const TTemplateConstRefList &tplList) {
  Int32 componentCount = tplList.size();
  // Pre-Allocate the rebined template and mask with regard to the spectrum size
  m_templatesRebined_bf.resize(componentCount);
  m_masksRebined_bf.resize(componentCount);
  m_spcSpectralAxis_restframe.SetSize(spectrum.GetSampleCount());

  for (Int32 ktpl = 0; ktpl < componentCount; ktpl++) {
    m_templatesRebined_bf[ktpl].m_ismCorrectionCalzetti =
        tplList[ktpl]->m_ismCorrectionCalzetti;
    m_templatesRebined_bf[ktpl].m_igmCorrectionMeiksin =
        tplList[ktpl]->m_igmCorrectionMeiksin;
  }
}

// @BasicFit@ in its current content and the commented code corresponds to
// @::computeModelSpectrum@ since there is no loops over igm/ism corrections

void COperatorTplcombination::BasicFit(
    const CSpectrum &spectrum, const TTemplateConstRefList &tplList,
    const TFloat64Range &lambdaRange, Float64 redshift,
    Float64 overlapThreshold, STplcombination_basicfitresult &fittingResults,
    Float64 forcedAmplitude, bool opt_extinction, bool opt_dustFitting,
    CMask spcMaskAdditional, const CPriorHelper::TPriorEList &logpriore,
    const TInt32List &MeiksinList, const TInt32List &EbmvList) {
  Log.LogDebug(Formatter() << " BasicFit - for z=" << redshift);

  boost::chrono::thread_clock::time_point start_prep =
      boost::chrono::thread_clock::now();

  bool chisquareSetAtLeastOnce = false;

  const CSpectrumSpectralAxis &spcSpectralAxis = spectrum.GetSpectralAxis();
  const CSpectrumFluxAxis &spcFluxAxis = spectrum.GetFluxAxis();
  const CSpectrumNoiseAxis &spcError = spcFluxAxis.GetError();

  if (spcMaskAdditional.GetMasksCount() != spcFluxAxis.GetSamplesCount())
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "spcMaskAdditional does not "
                          "have the same size as the spectrum flux vector... ("
                       << spcMaskAdditional.GetMasksCount() << " vs "
                       << spcFluxAxis.GetSamplesCount() << ")");

  TFloat64Range currentRange;
  RebinTemplate(spectrum, tplList, redshift, lambdaRange, currentRange,
                fittingResults.overlapFraction.front(), overlapThreshold);

  Int32 kStart = -1, kEnd = -1, kIgmEnd = -1;
  // I consider here that all templates share the same spectralAxis
  currentRange.getClosedIntervalIndices(
      m_templatesRebined_bf.front().GetSpectralAxis().GetSamplesVector(),
      kStart, kEnd);

  Int32 kStart_model =
      kStart; // mainly used at high redshifts, when desextincting spectrum is
              // happenning with null coeffs
  if (opt_extinction || opt_dustFitting)
    for (auto &tpl : m_templatesRebined_bf)
      tpl.InitIsmIgmConfig(kStart, kEnd, redshift);

  if (opt_extinction)
    kIgmEnd = m_templatesRebined_bf.front().GetIgmEndIndex();

  // determine min and max value of ebmv coeff
  Int32 nISM = EbmvList.size();
  Int32 nIGM = MeiksinList.size();

  Int32 iEbmvCoeffMin = EbmvList.front();
  Int32 iEbmvCoeffMax = EbmvList.back();

  // Linear fit
  Int32 n = kEnd - kStart + 1;
  Log.LogDebug(
      Formatter() << " prep. linear fitting with n=" << n
                  << " "
                     "samples in the clamped lambdarange spectrum (imin="
                  << kStart
                  << ", "
                     "lbda_min="
                  << spcSpectralAxis[kStart] << " - imax=" << kEnd
                  << ", lbda_max=" << spcSpectralAxis[kEnd] << ")");

  gsl_matrix *X, *cov;
  gsl_vector *y, *w, *c;

  Int32 nddl = tplList.size();
  X = gsl_matrix_alloc(n, nddl);
  y = gsl_vector_alloc(n);
  w = gsl_vector_alloc(n);
  c = gsl_vector_alloc(nddl);
  cov = gsl_matrix_alloc(nddl, nddl);

  // Normalizing factor
  Float64 normFactor = GetNormFactor(spcFluxAxis, kStart, n);

  Log.LogDetail(Formatter() << " Linear fitting, found "
                               "normalization Factor="
                            << normFactor);

  bool option_igmFastProcessing =
      (MeiksinList.size() == 1 ? false : true); // TODO
  bool igmLoopUseless_WavelengthRange = false;
  fittingResults.chiSquare = INFINITY; // final best Xi2 value
  Float64 chisq, SNR;
  Float64 dtd_complement =
      0.; // value mainly relevant with DisextinctData method

  // create a template with cte flux = 1, to be used only when disextincting
  // data and noise
  CTemplate identityTemplate(
      "identity", "idle", m_templatesRebined_bf.front().GetSpectralAxis(),
      std::move(CSpectrumFluxAxis(
          m_templatesRebined_bf.front().GetSampleCount(), 1)));
  // Prepare the fit data, once for all
  Float64 yi;
  Float64 ei;
  for (Int32 i = 0; i < n; i++) {
    yi = spcFluxAxis[i + kStart] / normFactor;
    ei = spcError[i + kStart] / normFactor;

    gsl_vector_set(y, i, yi);              // y[i] = yi
    gsl_vector_set(w, i, 1.0 / (ei * ei)); // w[i] = 1/(ei*ei)
  }
  for (Int32 kigm = 0; kigm < nIGM; kigm++) {
    if (igmLoopUseless_WavelengthRange) {
      // Now copy from the already calculated k>0 igm values
      for (auto &Xi2 : fittingResults.ChiSquareInterm)
        Xi2 = TFloat64List(nIGM, Xi2.front());
      for (auto &ismCoeffs : fittingResults.IsmCalzettiCoeffInterm)
        ismCoeffs = TFloat64List(nIGM, ismCoeffs.front());
      for (auto &igmCoeffs : fittingResults.IgmMeiksinIdxInterm)
        igmCoeffs = TInt32List(nIGM, igmCoeffs.front());
      break;
    }
    Int32 meiksinIdx = MeiksinList[kigm];

    bool igmCorrectionAppliedOnce = false;
    // applyMeiksin on all templates
    if (opt_extinction) {
      if (currentRange.GetBegin() / (1 + redshift) >
          RESTLAMBDA_LYA) // bug??: why dividing by (1+z) if currentRange is
                          // already in restframe??
        // if(COperatorTemplateFitting::CheckLyaIsInCurrentRange(currentRange))
        igmCorrectionAppliedOnce = false;
      else {
        igmCorrectionAppliedOnce =
            true; // since we are multiplying the bools, rather init to true
        for (Int32 iddl = 0; iddl < nddl; iddl++) {
          igmCorrectionAppliedOnce &=
              m_templatesRebined_bf[iddl].ApplyMeiksinCoeff(meiksinIdx);
        }
      }
      if (!igmCorrectionAppliedOnce)
        igmLoopUseless_WavelengthRange = true;
    }

    for (Int32 kEbmv_ = 0; kEbmv_ < nISM; kEbmv_++) {
      Int32 kEbmv = EbmvList[kEbmv_];
      Float64 coeffEBMV = -1.0; // no ism by default
      // apply ism on all templates, once for all
      if (opt_dustFitting) {
        coeffEBMV =
            m_templatesRebined_bf.front().m_ismCorrectionCalzetti->GetEbmvValue(
                kEbmv);
        for (Int32 iddl = 0; iddl < nddl; iddl++)
          m_templatesRebined_bf[iddl].ApplyDustCoeff(kEbmv);
      }

      // preparing the template matrix for computing the Xi2 value
      // the fastigm has its effect mainly on the number of spectrum to consider
      // for computing the amplitudes and then the fit
      for (Int32 i = 0; i < n; i++) {
        for (Int32 iddl = 0; iddl < nddl; iddl++) {
          Float64 fval = m_templatesRebined_bf[iddl].GetFluxAxis()[i + kStart];
          gsl_matrix_set(X, i, iddl, fval); // i.e., X[i,iddl]=fval -> X is an
                                            // extract of m_templatesRebinned_bf
        }
      }
      // Now fitting
      boost::chrono::thread_clock::time_point stop_prep =
          boost::chrono::thread_clock::now();
      Float64 duration_prep =
          boost::chrono::duration_cast<boost::chrono::microseconds>(stop_prep -
                                                                    start_prep)
              .count();
      Log.LogDebug(Formatter() << " Linear fitting, preparation "
                                  "time = "
                               << duration_prep << " microsec");
      boost::chrono::thread_clock::time_point start_fit =
          boost::chrono::thread_clock::now();
      {
        gsl_multifit_linear_workspace *work =
            gsl_multifit_linear_alloc(n, nddl);
        gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
        gsl_multifit_linear_free(work);
      }
      chisq += dtd_complement; // complement with the dtd value, only in the
                               // case of DisextinctData.
      //
      boost::chrono::thread_clock::time_point stop_fit =
          boost::chrono::thread_clock::now();
      Float64 duration_fit =
          boost::chrono::duration_cast<boost::chrono::microseconds>(stop_fit -
                                                                    start_fit)
              .count();
      Log.LogDebug(Formatter()
                   << " Linear fitting, fit = " << duration_fit << " microsec");
      boost::chrono::thread_clock::time_point start_postprocess =
          boost::chrono::thread_clock::now();

#define C(i) (gsl_vector_get(c, (i)))
#define COV(i, j) (gsl_matrix_get(cov, (i), (j)))
      Log.LogDebug(Formatter() << "# best fit: Y = " << C(0) << " X1 + " << C(1)
                               << " X2 ...");
      Log.LogDebug(Formatter() << "# covariance matrix:");
      Log.LogDebug(Formatter() << "[");
      Log.LogDebug(Formatter() << "  " << COV(0, 0) << ", " << COV(0, 1));
      Log.LogDebug(Formatter() << "  " << COV(1, 0) << "," << COV(1, 1));
      Log.LogDebug(Formatter() << "]");
      Log.LogDebug(Formatter() << "# chisq/n = " << chisq / n);

      for (Int32 iddl = 0; iddl < nddl; iddl++) {
        Float64 a = gsl_vector_get(c, iddl) * normFactor;
        Log.LogDebug(Formatter() << "# Found amplitude " << iddl << ": " << a
                                 << " +- " << COV(iddl, iddl) * normFactor);
      }

      // save the fitted amps and fitErrors, etc...
      Float64 a, err2, err, sA = 0, sE = 0;
      for (Int32 iddl = 0; iddl < nddl; iddl++) {
        a = gsl_vector_get(c, iddl);
        fittingResults.fittingAmplitudes[iddl] = a * normFactor;
        err2 = COV(iddl, iddl);
        err = std::sqrt(err2);
        fittingResults.fittingAmplitudeErrors[iddl] = err * normFactor;
        fittingResults.fittingAmplitudeSigmas[iddl] = a / err;
        sA += a * a;
        sE += err2;
      }
      SNR = std::sqrt(sA / sE);
      fittingResults.fittingAmplitudesInterm[kEbmv_][kigm] =
          fittingResults.fittingAmplitudes; // saving

      // save covariance matrix into MtM
      for (Int32 iddl = 0; iddl < nddl; iddl++) {
        for (Int32 iddc = 0; iddc < nddl; iddc++) {
          fittingResults.COV[iddl][iddc] = COV(iddl, iddc);
        }
      }
      if (fittingResults.fittingAmplitudes.size() != nddl) {
        Log.LogDebug(Formatter() << " Found nfittedamps(="
                                 << fittingResults.fittingAmplitudes.size()
                                 << ") "
                                    "different than nddl(="
                                 << nddl << ")");
      }

      if (chisq < fittingResults.chiSquare) {
        fittingResults.chiSquare = chisq;
        fittingResults.SNR = SNR;
        fittingResults.meiksinIdx =
            igmCorrectionAppliedOnce ? meiksinIdx : undefIdx;
        fittingResults.ebmvCoef = coeffEBMV;
        chisquareSetAtLeastOnce = true;
      }

      // save the interm chisquares in the intermediate vector
      fittingResults.ChiSquareInterm[kEbmv_][kigm] = chisq;
      fittingResults.IsmCalzettiCoeffInterm[kEbmv_][kigm] = coeffEBMV;
      fittingResults.IgmMeiksinIdxInterm[kEbmv_][kigm] =
          igmCorrectionAppliedOnce ? meiksinIdx : undefIdx;

      boost::chrono::thread_clock::time_point stop_postprocess =
          boost::chrono::thread_clock::now();
      Float64 duration_postprocess =
          boost::chrono::duration_cast<boost::chrono::microseconds>(
              stop_postprocess - start_postprocess)
              .count();
      Log.LogDebug(Formatter() << " Linear fitting, postprocess = "
                                  ""
                               << duration_postprocess << " microsec");

    } // end iterating over ISM
  }   // end iterating over IGM

  // fittingResults.modelSpectrum =
  // CSpectrum(CSpectrumSpectralAxis(std::move(spc_extract)),
  // CSpectrumFluxAxis(std::move(modelFlux)));

  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_vector_free(w);
  gsl_vector_free(c);
  gsl_matrix_free(cov);

  if (!chisquareSetAtLeastOnce) {
    THROWG(ErrorCode::INVALID_MERIT_VALUES,
           Formatter() << "Not even one single valid fit/merit value found");
  }
}

void COperatorTplcombination::RebinTemplate(
    const CSpectrum &spectrum, const TTemplateConstRefList &tplList,
    Float64 redshift, const TFloat64Range &lambdaRange,
    TFloat64Range &currentRange, Float64 &overlapFraction,
    const Float64 overlapThreshold) {
  Float64 onePlusRedshift = 1.0 + redshift;

  // shift lambdaRange backward to be in restframe
  TFloat64Range spcLambdaRange_restframe;
  TFloat64Range lambdaRange_restframe(lambdaRange.GetBegin() / onePlusRedshift,
                                      lambdaRange.GetEnd() / onePlusRedshift);

  // redshift in restframe the tgtSpectralAxis,
  m_spcSpectralAxis_restframe = spectrum.GetSpectralAxis().ShiftByWaveLength(
      onePlusRedshift, CSpectrumSpectralAxis::nShiftBackward);
  m_spcSpectralAxis_restframe.ClampLambdaRange(lambdaRange_restframe,
                                               spcLambdaRange_restframe);

  TFloat64Range intersectedAllLambdaRange(spcLambdaRange_restframe);

  // Now interpolating all the templates
  Log.LogDebug(" BasicFit - interpolating");
  for (Int32 ktpl = 0; ktpl < tplList.size(); ktpl++) {
    const CSpectrumSpectralAxis &tplSpectralAxis =
        tplList[ktpl]->GetSpectralAxis();
    const CSpectrumFluxAxis &tplFluxAxis = tplList[ktpl]->GetFluxAxis();

    // Compute clamped lambda range over template
    TFloat64Range tplLambdaRange;
    tplSpectralAxis.ClampLambdaRange(lambdaRange_restframe, tplLambdaRange);

    // if there is any intersection between the lambda range of the spectrum and
    // the lambda range of the template Compute the intersected range
    TFloat64Range intersectedLambdaRange(0.0, 0.0);
    TFloat64Range::Intersect(tplLambdaRange, spcLambdaRange_restframe,
                             intersectedLambdaRange);

    // find lambda range intersection common to all templates
    intersectedAllLambdaRange.IntersectWith(intersectedLambdaRange);

    CTemplate &itplTplSpectrum = m_templatesRebined_bf[ktpl];
    CMask &itplMask = m_masksRebined_bf[ktpl];

    tplList[ktpl]->Rebin(intersectedLambdaRange, m_spcSpectralAxis_restframe,
                         itplTplSpectrum, itplMask);

    const CSpectrumSpectralAxis &itplTplSpectralAxis =
        itplTplSpectrum.GetSpectralAxis();
    Log.LogDebug(
        Formatter()
        << " Rebinned template #" << ktpl
        << " has n=" << itplTplSpectralAxis.GetSamplesCount()
        << " samples in "
           "lambdaRange: "
        << itplTplSpectralAxis[0] << "  - "
        << itplTplSpectralAxis[itplTplSpectralAxis.GetSamplesCount() - 1]);

    overlapFraction =
        m_spcSpectralAxis_restframe.IntersectMaskAndComputeOverlapFraction(
            lambdaRange_restframe, itplMask);

    // Check for overlap rate
    if (overlapFraction < overlapThreshold || overlapFraction <= 0.0) {
      THROWG(ErrorCode::OVERLAPFRACTION_NOTACCEPTABLE,
             Formatter() << "overlapFraction of " << overlapFraction);
    }
  }
  currentRange = intersectedAllLambdaRange;
  return;
}

/**
 * \brief
 *
 * input: if additional_spcMasks size is 0, no additional mask will be used,
 *otherwise its size should match the redshifts list size
 * @lambdaRange is not clamped
 **/
std::shared_ptr<COperatorResult> COperatorTplcombination::Compute(
    const CSpectrum &spectrum, const TTemplateConstRefList &tplList,
    const TFloat64Range &lambdaRange, const TFloat64List &redshifts,
    Float64 overlapThreshold, const std::vector<CMask> &additional_spcMasks,
    const std::string &opt_interp, bool opt_extinction, bool opt_dustFitting,
    const CPriorHelper::TPriorZEList &logpriorze, Int32 FitEbmvIdx,
    Int32 FitMeiksinIdx) {
  Int32 componentCount = tplList.size();
  Log.LogInfo(Formatter() << " starting computation with N-template = "
                          << componentCount);

  for (Int32 ktpl = 0; ktpl < componentCount; ktpl++) {
    if (opt_dustFitting && tplList[ktpl]->CalzettiInitFailed()) {
      THROWG(ErrorCode::INTERNAL_ERROR, "ISM is not initialized");
    }
    if (opt_extinction && tplList[ktpl]->MeiksinInitFailed()) {
      THROWG(ErrorCode::INTERNAL_ERROR, "IGM is not initialized");
    }
  }

  Log.LogDebug(Formatter() << " allocating memory for buffers (N = "
                           << componentCount << ")");

  BasicFit_preallocateBuffers(spectrum, tplList);

  // sort the redshift and keep track of the indexes
  TFloat64List sortedRedshifts;
  TFloat64List sortedIndexes;
  // This is a vector of {value,index} pairs
  vector<pair<Float64, Int32>> vp;
  vp.reserve(redshifts.size());
  for (Int32 i = 0; i < redshifts.size(); i++) {
    vp.push_back(make_pair(redshifts[i], i));
  }
  std::sort(vp.begin(), vp.end());
  for (Int32 i = 0; i < vp.size(); i++) {
    sortedRedshifts.push_back(vp[i].first);
    sortedIndexes.push_back(vp[i].second);
  }

  Log.LogDebug(" prepare the results");

  TInt32List MeiksinList;
  TInt32List EbmvList;
  TIgmIsmIdxs igmIsmIdxs;
  m_templatesRebined_bf.front().GetIsmIgmIdxList(
      opt_extinction, opt_dustFitting, igmIsmIdxs, FitEbmvIdx, FitMeiksinIdx);
  Int32 MeiksinListSize = igmIsmIdxs.igmIdxs.size();
  Int32 EbmvListSize = igmIsmIdxs.ismIdxs.size();
  Log.LogDebug(Formatter() << " prepare N ism coeffs = " << EbmvListSize);
  Log.LogDebug(Formatter() << " prepare N igm coeffs = " << MeiksinListSize);
  std::shared_ptr<CTplCombinationResult> result =
      make_shared<CTplCombinationResult>(sortedRedshifts.size(), EbmvListSize,
                                         MeiksinListSize, componentCount);
  result->Redshifts = sortedRedshifts;

  // default mask
  bool useDefaultMask =
      additional_spcMasks.size() != sortedRedshifts.size() ? true : false;
  CMask default_spcMask(spectrum.GetSampleCount());
  if (useDefaultMask)
    for (Int32 km = 0; km < default_spcMask.GetMasksCount(); km++)
      default_spcMask[km] = 1.0;

  if (additional_spcMasks.size() != sortedRedshifts.size() &&
      additional_spcMasks.size() != 0)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "masks-list and redshift size do not match: "
                       << additional_spcMasks.size()
                       << "!=" << sortedRedshifts.size());

  TFloat64Range clampedlambdaRange;
  spectrum.GetSpectralAxis().ClampLambdaRange(lambdaRange, clampedlambdaRange);

  for (Int32 i = 0; i < sortedRedshifts.size(); i++) {
    const CMask &additional_spcMask =
        useDefaultMask ? default_spcMask
                       : additional_spcMasks[sortedIndexes[i]];

    const CPriorHelper::TPriorEList &logp =
        logpriorze.size() > 0 && logpriorze.size() == sortedRedshifts.size()
            ? logpriorze[i]
            : CPriorHelper::TPriorEList();

    Float64 redshift = result->Redshifts[i];

    STplcombination_basicfitresult fittingResults(EbmvListSize, MeiksinListSize,
                                                  componentCount);

    BasicFit(spectrum, tplList, clampedlambdaRange, redshift, overlapThreshold,
             fittingResults, -1, opt_extinction, opt_dustFitting,
             additional_spcMask, logp, igmIsmIdxs.igmIdxs, igmIsmIdxs.ismIdxs);

    result->ChiSquare[i] = fittingResults.chiSquare;
    result->Overlap[i] = fittingResults.overlapFraction;
    result->FitAmplitude[i] = fittingResults.fittingAmplitudes;
    result->FitAmplitudeSigma[i] = fittingResults.fittingAmplitudeSigmas;
    result->FitAmplitudeError[i] = fittingResults.fittingAmplitudeErrors;
    result->SNR[i] = fittingResults.SNR;
    result->FitCOV[i] = fittingResults.COV;
    // result->LogPrior[i]=NAN: //not yet calculated
    result->FitEbmvCoeff[i] = fittingResults.ebmvCoef;
    result->FitMeiksinIdx[i] = fittingResults.meiksinIdx;
    result->ChiSquareIntermediate[i] = fittingResults.ChiSquareInterm;
    result->IsmEbmvCoeffIntermediate[i] = fittingResults.IsmCalzettiCoeffInterm;
    result->IgmMeiksinIdxIntermediate[i] = fittingResults.IgmMeiksinIdxInterm;
  }

  // overlap warning
  Float64 overlapValidInfZ = -1;
  for (Int32 i = 0; i < sortedRedshifts.size(); i++) {
    if (result->Overlap[i].front() >= overlapThreshold) {
      overlapValidInfZ = sortedRedshifts[i];
      break;
    }
  }
  Float64 overlapValidSupZ = -1;
  for (Int32 i = sortedRedshifts.size() - 1; i >= 0; i--) {
    if (result->Overlap[i].front() >= overlapThreshold) {
      overlapValidSupZ = sortedRedshifts[i];
      break;
    }
  }
  if (overlapValidInfZ != sortedRedshifts[0] ||
      overlapValidSupZ != sortedRedshifts[sortedRedshifts.size() - 1]) {
    Log.LogInfo(Formatter() << " overlap warning for: minz=" << overlapValidInfZ
                            << ", maxz=" << overlapValidSupZ);
  }

  // estimate CstLog for PDF estimation
  result->CstLog = EstimateLikelihoodCstLog(spectrum, clampedlambdaRange);

  // Deallocate the rebined template and mask buffers
  m_templatesRebined_bf.clear();
  m_masksRebined_bf.clear();

  return result;
}

std::shared_ptr<CModelSpectrumResult>
COperatorTplcombination::ComputeSpectrumModel(
    const CSpectrum &spectrum, const TTemplateConstRefList &tplList,
    Float64 redshift, Float64 ebmvCoef, Int32 meiksinIdx,
    const TFloat64List &amplitudes, const TFloat64Range &lambdaRange,
    const Float64 overlapThreshold) {
  Log.LogDetail(
      Formatter()
      << "  Operator-COperatorTplCombination: building spectrum model "
         "tptCombination for candidate Zcand="
      << redshift);

  BasicFit_preallocateBuffers(spectrum, tplList);
  Int32 nddl = tplList.size();

  // Estatus status;
  Float64 overlapFraction = 0.0;
  TFloat64Range currentRange;
  RebinTemplate(spectrum, tplList, redshift, lambdaRange, currentRange,
                overlapFraction, overlapThreshold);
  Int32 kStart = -1, kEnd = -1, kIgmEnd = -1;

  currentRange.getClosedIntervalIndices(
      m_templatesRebined_bf.front().GetSpectralAxis().GetSamplesVector(),
      kStart, kEnd);

  // create identityTemplate on which we apply meiksin and ism, once for all
  // tpllist
  CTemplate identityTemplate(
      "identity", "idle", m_templatesRebined_bf.front().GetSpectralAxis(),
      CSpectrumFluxAxis(m_templatesRebined_bf.front().GetSampleCount(), 1));

  if ((ebmvCoef > 0.) || (meiksinIdx > -1)) {
    identityTemplate.InitIsmIgmConfig(kStart, kEnd, redshift);
  }

  if (ebmvCoef > 0.) {
    Int32 idxEbmv = -1;
    idxEbmv = identityTemplate.m_ismCorrectionCalzetti->GetEbmvIndex(ebmvCoef);

    if (idxEbmv != -1)
      identityTemplate.ApplyDustCoeff(idxEbmv);
  }

  if (meiksinIdx > -1) {
    identityTemplate.ApplyMeiksinCoeff(meiksinIdx);
  }

  const CSpectrumFluxAxis &extinction = identityTemplate.GetFluxAxis();

  Int32 modelSize = spectrum.GetSampleCount();

  CSpectrumFluxAxis modelFlux(modelSize, 0.0);
  for (Int32 iddl = 0; iddl < nddl; iddl++) {
    const CSpectrumFluxAxis &tmp = m_templatesRebined_bf[iddl].GetFluxAxis();
    for (Int32 k = 0; k < modelSize; k++) {
      modelFlux[k] += amplitudes[iddl] * tmp[k] * extinction[k];
    }
  }
  CSpectrumSpectralAxis modelSpcAxis =
      spectrum.GetSpectralAxis().ShiftByWaveLength(
          (1.0 + redshift), CSpectrumSpectralAxis::nShiftForward);

  // Deallocate the rebined template and mask buffers
  m_templatesRebined_bf.clear();
  m_masksRebined_bf.clear();
  std::shared_ptr<CModelSpectrumResult> ret =
      std::make_shared<CModelSpectrumResult>();
  ret->addModel(CSpectrum(std::move(modelSpcAxis), std::move(modelFlux)), "");
  return ret;
}

/**
 * \brief this function estimates the likelihood_cstLog term withing the
 *wavelength range
 **/
Float64 COperatorTplcombination::EstimateLikelihoodCstLog(
    const CSpectrum &spectrum, const TFloat64Range &lambdaRange) {
  const CSpectrumSpectralAxis &spcSpectralAxis = spectrum.GetSpectralAxis();
  const TFloat64List &error =
      spectrum.GetFluxAxis().GetError().GetSamplesVector();

  Int32 numDevs = 0;
  Float64 cstLog = 0.0;
  Float64 sumLogNoise = 0.0;

  Int32 imin;
  Int32 imax;
  lambdaRange.getClosedIntervalIndices(spcSpectralAxis.GetSamplesVector(), imin,
                                       imax);
  for (Int32 j = imin; j <= imax; j++) {
    numDevs++;
    sumLogNoise += log(error[j]);
  }

  cstLog = -numDevs * 0.5 * log(2 * M_PI) - sumLogNoise;

  return cstLog;
}

Float64
COperatorTplcombination::GetNormFactor(const CSpectrumFluxAxis spcFluxAxis,
                                       Int32 kStart, Int32 n) {
  Float64 maxabsval = DBL_MIN;
  for (Int32 k = 0; k < n; k++) {
    if (maxabsval < std::abs(spcFluxAxis[k + kStart])) {
      maxabsval = std::abs(spcFluxAxis[k + kStart]);
    }
  }
  return maxabsval;
}
