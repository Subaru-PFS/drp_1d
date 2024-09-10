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
#include "RedshiftLibrary/operator/powerlawresult.h"

using namespace NSEpic;

CPowerLawResult::CPowerLawResult(Int32 n)
    : COperatorResult("CPowerLawResult"), Redshifts(n), ChiSquare(n),
      fluxError(n), coefs(n), FitEbmvCoeff(n), FitMeiksinIdx(n),
      ChiSquareIntermediate(n), IsmEbmvCoeffIntermediate(n),
      IgmMeiksinIdxIntermediate(n), SNR(n) {}

void CPowerLawResult::set_at_redshift(Int32 zIdx, TPowerLawResult result) {
  ChiSquare[zIdx] = result.chiSquare;
  coefs[zIdx] = result.coefs;
  FitEbmvCoeff[zIdx] = result.ebmvCoef;
  FitMeiksinIdx[zIdx] = result.meiksinIdx;
  SNR[zIdx] = SNRCalculation(result.coefs);

  // [z][Calzetti][Meiksin])
  T2DList<Float64> ChiSquareIntermediate;
  T2DList<Float64> IgmMeiksinIdxIntermediate;
  T2DList<Float64> IsmEbmvCoeffIntermediate;
}


Float64 CPowerLawResult::SNRCalculation(TPowerLawCoefsPair const &coefs) const {
  return std::max(coefs.first.a / coefs.first.stda, coefs.second.a / coefs.second.stda);
}
