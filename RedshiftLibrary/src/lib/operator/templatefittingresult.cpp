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
#include "RedshiftLibrary/operator/templatefittingresult.h"
#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/operator/templatefitting.h"
#include "RedshiftLibrary/operator/twopass.h"

using namespace NSEpic;

CTemplateFittingResult::CTemplateFittingResult(Int32 n)
    : CTwoPassResult("CTemplateFittingResult"), ChiSquare(n),
      ReducedChiSquare(n), ChiSquarePhot(n), FitAmplitude(n),
      FitAmplitudeError(n), FitAmplitudeSigma(n), FitEbmvCoeff(n),
      FitMeiksinIdx(n), FitDtM(n), FitMtM(n), LogPrior(n),
      ChiSquareIntermediate(n), IsmEbmvCoeffIntermediate(n),
      IgmMeiksinIdxIntermediate(n), SNR(n), Overlap(n) {
  Redshifts.resize(n);
}

CTemplateFittingResult::CTemplateFittingResult(Int32 n, Int32 EbmvListSize,
                                               Int32 MeiksinListSize)
    : CTemplateFittingResult::CTemplateFittingResult(n) {
  ChiSquareIntermediate.assign(
      n, std::vector<TFloat64List>(EbmvListSize,
                                   TFloat64List(MeiksinListSize, DBL_MAX)));
  IsmEbmvCoeffIntermediate.assign(
      n, std::vector<TFloat64List>(EbmvListSize,
                                   TFloat64List(MeiksinListSize, NAN)));
  IgmMeiksinIdxIntermediate.assign(
      n, std::vector<TInt32List>(EbmvListSize,
                                 TInt32List(MeiksinListSize, undefIdx)));
}

void CTemplateFittingResult::set_at_redshift(Int32 i, TFittingIsmIgmResult val,
                                             Int32 igmIdx, Int32 ismIdx) {
  ChiSquare[i] = val.chiSquare;
  ReducedChiSquare[i] = val.reducedChisquare;
  ChiSquarePhot[i] = val.chiSquare_phot;
  FitAmplitude[i] = val.ampl;
  FitAmplitudeError[i] = val.ampl_err;
  FitAmplitudeSigma[i] = val.ampl_sigma;
  FitEbmvCoeff[i] = val.ebmvCoef;
  FitMeiksinIdx[i] = val.meiksinIdx;
  FitDtM[i] = val.cross_result.sumCross;
  FitMtM[i] = val.cross_result.sumT;

  SNR[i] = SNRCalculation(FitDtM[i], FitMtM[i]);
  LogPrior[i] = val.logprior;
  Overlap[i] = val.overlapFraction;

  bool isNotFromLineModelSolve =
      ismIdx != -1 && igmIdx != -1 &&
      ssize(ChiSquareIntermediate[i]) > ismIdx &&
      ssize(ChiSquareIntermediate[i][ismIdx]) > igmIdx;
  bool isSecondPass = igmIdx != undefIdx && ismIdx != undefIdx;

  if (isSecondPass && isNotFromLineModelSolve) {
    ChiSquareIntermediate[i][ismIdx][igmIdx] =
        std::move(val.ChiSquareInterm[0][0]);
    IsmEbmvCoeffIntermediate[i][ismIdx][igmIdx] =
        std::move(val.IsmCalzettiCoeffInterm[0][0]);
    IgmMeiksinIdxIntermediate[i][ismIdx][igmIdx] =
        std::move(val.IgmMeiksinIdxInterm[0][0]);
  } else {
    ChiSquareIntermediate[i] = std::move(val.ChiSquareInterm);
    IsmEbmvCoeffIntermediate[i] =
        std::move(val.IsmCalzettiCoeffInterm); // TODO see if useful ?
    IgmMeiksinIdxIntermediate[i] =
        std::move(val.IgmMeiksinIdxInterm); // TODO see if useful ?
  }
}

Float64 CTemplateFittingResult::SNRCalculation(Float64 dtm, Float64 mtm) {
  Float64 snr;
  if (mtm <= 0) {
    snr = NAN;
  } else if (dtm <= 0) {
    snr = 0;
  } else {
    snr = dtm / std::sqrt(mtm);
  }
  return snr;
}

void CTemplateFittingResult::updateVectors(
    TsecondPassIndices const &secondPassIndices) {
  auto const &[insertionIdx, initialInsertionIdx, overwrittenSourceIndices,
               count] = secondPassIndices;
  Int32 ndup = overwrittenSourceIndices.size();
  insertWithDuplicates<Float64>(ChiSquare, insertionIdx, count, NAN, ndup);
  insertWithDuplicates<Float64>(ReducedChiSquare, insertionIdx, count, NAN,
                                ndup);
  insertWithDuplicates<Float64>(ChiSquarePhot, insertionIdx, count, NAN, ndup);
  insertWithDuplicates<Float64>(FitAmplitude, insertionIdx, count, NAN, ndup);
  insertWithDuplicates<Float64>(FitAmplitudeError, insertionIdx, count, NAN,
                                ndup);
  insertWithDuplicates<Float64>(FitAmplitudeSigma, insertionIdx, count, NAN,
                                ndup);
  insertWithDuplicates<Float64>(FitEbmvCoeff, insertionIdx, count, NAN, ndup);
  insertWithDuplicates<Int32>(FitMeiksinIdx, insertionIdx, count, -1, ndup);
  insertWithDuplicates<Float64>(FitDtM, insertionIdx, count, NAN, ndup);
  insertWithDuplicates<Float64>(FitMtM, insertionIdx, count, NAN, ndup);
  insertWithDuplicates<Float64>(LogPrior, insertionIdx, count, NAN, ndup);
  insertWithDuplicates<Float64>(SNR, insertionIdx, count, NAN, ndup);

  insertWithDuplicates<std::vector<Float64>>(
      Overlap, insertionIdx, count,
      std::vector<Float64>(Overlap[0].size(), NAN), ndup);

  auto const chi2ToInsert =
      interpolatedChiSquareIntermediate(secondPassIndices);
  insertWithDuplicates(ChiSquareIntermediate, insertionIdx, chi2ToInsert, ndup);

  size_t nIsm = 0;
  size_t nIgm = 0;
  if (ChiSquareIntermediate.size() > 0) {
    nIsm = ChiSquareIntermediate[0].size();
    if (ChiSquareIntermediate[0].size() > 0)
      nIgm = ChiSquareIntermediate[0][0].size();
  }

  T3DList<Float64> toInsertFloat(
      count, T2DList<Float64>(nIsm, std::vector<Float64>(nIgm, NAN)));
  insertWithDuplicates<std::vector<std::vector<Float64>>>(
      IsmEbmvCoeffIntermediate, insertionIdx, toInsertFloat, ndup);

  T3DList<Int32> toInsertInt(
      count, T2DList<Int32>(nIsm, std::vector<Int32>(nIgm, -1)));
  insertWithDuplicates<std::vector<std::vector<Int32>>>(
      IgmMeiksinIdxIntermediate, insertionIdx, toInsertInt, ndup);
}

// interpolated values from ChiSquareIntermediate
T3DList<Float64> CTemplateFittingResult::interpolatedChiSquareIntermediate(
    TsecondPassIndices const &secondPassIndices) {
  auto const &[_, initialInsertionIdx, overwrittenSourceIndices, count] =
      secondPassIndices;

  Int32 ndup = overwrittenSourceIndices.size();
  auto const &[nIsm, nIgm] = getIsmIgmSizes(ChiSquareIntermediate);

  T3DList<Float64> interpolatedToInsert(
      count, T2DList<Float64>(nIsm, TList<Float64>(nIgm)));

  // left
  Int32 dupStart = 0;
  Int32 const firstSourceIdx = 0;
  Int32 const firstZidx = 0;
  if (overwrittenSourceIndices[dupStart] != firstSourceIdx) {
    Int32 idx = std::max(initialInsertionIdx - 1, firstZidx);
    Int32 interp_idx1 = firstSourceIdx;
    Int32 interp_idx2 = overwrittenSourceIndices[dupStart];
    interpolatedBetweenTwoSamples(
        ChiSquareIntermediate[idx], ChiSquareIntermediate[idx + 1],
        TInt32Range(interp_idx1, interp_idx2), interpolatedToInsert);
    ++dupStart;
  }

  // right
  // TODO, handle case ndup==1 ?
  Int32 dupEnd = ndup - 1;
  Int32 const lastSourceIdx = count - 1;
  Int32 const lastZidx = Int32(ChiSquareIntermediate.size()) - 1;
  if (overwrittenSourceIndices[dupEnd] != lastSourceIdx) {
    Int32 idx = std::min(initialInsertionIdx + dupEnd + 1, lastZidx);
    Int32 interp_idx1 = overwrittenSourceIndices[dupEnd];
    Int32 interp_idx2 = lastSourceIdx;
    interpolatedBetweenTwoSamples(
        ChiSquareIntermediate[idx], ChiSquareIntermediate[idx + 1],
        TInt32Range(interp_idx1, interp_idx2), interpolatedToInsert);
    --dupEnd;
  }

  // middle
  for (Int32 dup = dupStart; dup < dupEnd; ++dup) {
    Int32 idx = initialInsertionIdx + dup;
    Int32 interp_idx1 = overwrittenSourceIndices[dup];
    Int32 interp_idx2 = overwrittenSourceIndices[dup + 1];
    interpolatedBetweenTwoSamples(
        ChiSquareIntermediate[idx], ChiSquareIntermediate[idx + 1],
        TInt32Range(interp_idx1, interp_idx2), interpolatedToInsert);
  }

  return interpolatedToInsert;
}

void CTemplateFittingResult::interpolatedBetweenTwoSamples(
    T2DList<Float64> &array1, T2DList<Float64> &array2,
    TInt32Range const &destIdxRange, T3DList<Float64> dest) {
  auto const &[nIsm, nIgm] = getIsmIgmSizes(ChiSquareIntermediate);

  Int32 nInterpolate = destIdxRange.GetLength() + 1;
  for (Int32 ismIdx = 0; ismIdx < nIsm; ++ismIdx) {
    for (Int32 igmIdx = 0; igmIdx < nIgm; ++igmIdx) {
      Float64 v1 = array1[ismIdx][igmIdx];
      Float64 v2 = array2[ismIdx][igmIdx];
      TFloat64List linearVect = createLinearInterpVector(v1, v2, nInterpolate);
      for (Int32 zIdx : destIdxRange)
        dest[zIdx][ismIdx][igmIdx] = linearVect[zIdx];
    }
  }
}

std::pair<Int32, Int32> CTemplateFittingResult::getIsmIgmSizes(
    T3DList<Float64> const &ChiSquareIntermediate) {
  size_t nIsm = 0;
  size_t nIgm = 0;
  if (ChiSquareIntermediate.size() > 0) {
    nIsm = ChiSquareIntermediate[0].size();
    if (ChiSquareIntermediate[0].size() > 0)
      nIgm = ChiSquareIntermediate[0][0].size();
  }
  return std::make_pair(nIsm, nIgm);
}
