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
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;

CTemplateFittingResult::CTemplateFittingResult(Int32 n)
    : CTwoPassResult("CTemplateFittingResult"), ChiSquare(n),
      ReducedChiSquare(n), pValue(n), ChiSquarePhot(n), FitAmplitude(n),
      FitAmplitudeError(n), FitAmplitudeSigma(n), FitEbmvCoeff(n),
      FitMeiksinIdx(n), FitDtM(n), FitMtM(n), LogPrior(n), SNR(n),
      ChiSquareIntermediate(n), IsmEbmvIdxIntermediate(n),
      IgmMeiksinIdxIntermediate(n), Overlap(n) {
  Redshifts.resize(n);
}

CTemplateFittingResult::CTemplateFittingResult(Int32 n, Int32 EbmvListSize,
                                               Int32 MeiksinListSize)
    : CTemplateFittingResult::CTemplateFittingResult(n) {
  ChiSquareIntermediate.assign(
      n, std::vector<TFloat64List>(EbmvListSize,
                                   TFloat64List(MeiksinListSize, DBL_MAX)));
  IsmEbmvIdxIntermediate.assign(n, TInt32List(EbmvListSize, undefIdx));
  IgmMeiksinIdxIntermediate.assign(n, TInt32List(MeiksinListSize, undefIdx));
}

Int32 CTemplateFittingResult::getIsmIndexInIntermediate(Int32 ebmvIdx) const {
  return CIndexing<Int32>::getIndex(IsmEbmvIdxIntermediate.front(), ebmvIdx);
}

Int32 CTemplateFittingResult::getIgmIndexInIntermediate(Int32 zIdx,
                                                        Int32 igmIdx) const {
  return CIndexing<Int32>::getIndex(IgmMeiksinIdxIntermediate[zIdx], igmIdx);
}

void CTemplateFittingResult::set_at_redshift(Int32 i,
                                             TFittingIsmIgmResult val) {
  ChiSquare[i] = val.chiSquare;
  ReducedChiSquare[i] = val.reducedChiSquare;
  pValue[i] = val.pValue;
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

  auto const sizes = getIsmIgmSizes();
  auto const val_sizes = std::pair<Int32, Int32>(
      val.ChiSquareInterm.size(), val.ChiSquareInterm.front().size());
  bool has_values = !ChiSquareIntermediate[i].empty();
  bool insert = has_values && (sizes != val_sizes);
  if (insert) {

    // here we have only one (ism,igm) computed in val.
    Int32 ismIdx = getIsmIndexInIntermediate(val.IsmCalzettiIdxInterm.front());
    Int32 igmIdx =
        getIgmIndexInIntermediate(i, val.IgmMeiksinIdxInterm.front());

    ChiSquareIntermediate[i][ismIdx][igmIdx] =
        std::move(val.ChiSquareInterm[0][0]);
    IsmEbmvIdxIntermediate[i][ismIdx] =
        std::move(val.IsmCalzettiIdxInterm.front());
    IgmMeiksinIdxIntermediate[i][igmIdx] =
        std::move(val.IgmMeiksinIdxInterm.front());
  } else {
    ChiSquareIntermediate[i] = std::move(val.ChiSquareInterm);
    IsmEbmvIdxIntermediate[i] = std::move(val.IsmCalzettiIdxInterm);
    IgmMeiksinIdxIntermediate[i] = std::move(val.IgmMeiksinIdxInterm);
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
  auto const &[insertionIdx, overwrittenSourceIndices, count, largeStepFactor] =
      secondPassIndices;
  Int32 ndup = overwrittenSourceIndices.size();
  insertWithDuplicates<Float64>(ChiSquare, insertionIdx, count, NAN, ndup);
  insertWithDuplicates<Float64>(ReducedChiSquare, insertionIdx, count, NAN,
                                ndup);
  insertWithDuplicates<Float64>(pValue, insertionIdx, count, NAN, ndup);
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

  auto const chi2ToInsert = interpolateBetweenDuplicates(
      ChiSquareIntermediate, insertionIdx, overwrittenSourceIndices, count,
      largeStepFactor);
  insertWithDuplicates(ChiSquareIntermediate, insertionIdx, chi2ToInsert, ndup);

  auto const ismEbmvIdxToInsert = interpolateBetweenDuplicates(
      IsmEbmvIdxIntermediate, insertionIdx, overwrittenSourceIndices, count,
      largeStepFactor);
  insertWithDuplicates(IsmEbmvIdxIntermediate, insertionIdx, ismEbmvIdxToInsert,
                       ndup);

  auto const igmMeiksinIdxToInsert = interpolateBetweenDuplicates(
      IgmMeiksinIdxIntermediate, insertionIdx, overwrittenSourceIndices, count,
      largeStepFactor);
  insertWithDuplicates<TInt32List>(IgmMeiksinIdxIntermediate, insertionIdx,
                                   igmMeiksinIdxToInsert, ndup);
}

std::pair<Int32, Int32> CTemplateFittingResult::getIsmIgmSizes() const {
  size_t nIsm = 0;
  size_t nIgm = 0;
  if (!ChiSquareIntermediate.empty()) {
    nIsm = ChiSquareIntermediate.front().size();
    if (!ChiSquareIntermediate.front().empty())
      nIgm = ChiSquareIntermediate.front().front().size();
  }
  return std::make_pair(nIsm, nIgm);
}
