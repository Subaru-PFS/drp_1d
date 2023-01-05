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
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/operator/templatefitting.h"

using namespace NSEpic;

CTemplateFittingResult::CTemplateFittingResult(Int32 n)
    : COperatorResult("CTemplateFittingResult"), Redshifts(n), ChiSquare(n),
      ChiSquarePhot(n), FitAmplitude(n), FitAmplitudeError(n),
      FitAmplitudeSigma(n), FitEbmvCoeff(n), FitMeiksinIdx(n), FitDtM(n),
      FitMtM(n), LogPrior(n), ChiSquareIntermediate(n),
      IsmEbmvCoeffIntermediate(n), IgmMeiksinIdxIntermediate(n), SNR(n),
      Overlap(n), Status(n) {}

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

void CTemplateFittingResult::set_at_redshift(Int32 i,
                                             TFittingIsmIgmResult val) {
  ChiSquare[i] = val.chiSquare;
  ChiSquarePhot[i] = val.chiSquare_phot;
  FitAmplitude[i] = val.ampl;
  FitAmplitudeError[i] = val.ampl_err;
  FitAmplitudeSigma[i] = val.ampl_sigma;
  FitEbmvCoeff[i] = val.EbmvCoeff;
  FitMeiksinIdx[i] = val.MeiksinIdx;
  FitDtM[i] = val.sumCross;
  FitMtM[i] = val.sumT;
  LogPrior[i] = val.logprior;
  Overlap[i] = val.overlapRate;
  Status[i] = val.status;

  ChiSquareIntermediate[i] = std::move(val.ChiSquareInterm);
  IsmEbmvCoeffIntermediate[i] = std::move(val.IsmCalzettiCoeffInterm);
  IgmMeiksinIdxIntermediate[i] = std::move(val.IgmMeiksinIdxInterm);
}
