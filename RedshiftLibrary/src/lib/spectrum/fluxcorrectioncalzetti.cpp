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
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"
#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/log/log.h"

using namespace NSEpic;

CalzettiCorrection::CalzettiCorrection(TFloat64List _lbda,
                                       TFloat64List _fluxcorr)
    : lbda(_lbda), fluxcorr(_fluxcorr) {
  Float64 lbdaGridStep = lbda[1] - lbda[0];
  Float64 relativeLbdaGridStepTol = 1e-8;
  for (std::size_t t = 2; t < lbda.size(); t++) {
    Float64 lbda1 = lbda[t - 1];
    Float64 lbda2 = lbda[t];
    Float64 _lbdaGridStep;
    _lbdaGridStep = lbda2 - lbda1;
    Float64 relativeErrAbs =
        std::abs((_lbdaGridStep - lbdaGridStep) / lbdaGridStep);

    if (relativeErrAbs > relativeLbdaGridStepTol)
      THROWG(ErrorCode::BAD_CALZETTICORR, "lambdas are not regular sampled");
  }
  step = lbdaGridStep;
}

CSpectrumFluxCorrectionCalzetti::CSpectrumFluxCorrectionCalzetti(
    CalzettiCorrection _calzettiCorr, Float64 ebmv_start, Float64 ebmv_step,
    Float64 ebmv_n)
    : m_dataCalzetti(std::move(_calzettiCorr.fluxcorr)), m_nEbmvCoeff(ebmv_n),
      m_EbmvCoeffStep(ebmv_step), m_EbmvCoeffStart(ebmv_start),
      m_LambdaMin(_calzettiCorr.lbda.front()),
      m_LambdaMax(_calzettiCorr.lbda.back()), m_step(_calzettiCorr.step) {
  m_dataDustCoeff.resize(m_nEbmvCoeff * m_dataCalzetti.size());
  for (Int32 kDust = 0; kDust < m_nEbmvCoeff; kDust++) {
    Float64 coeffEBMV = GetEbmvValue(kDust);
    for (Int32 kCalzetti = 0; kCalzetti < ssize(m_dataCalzetti); kCalzetti++) {
      m_dataDustCoeff[kDust * m_dataCalzetti.size() + kCalzetti] =
          pow(10.0, -0.4 * m_dataCalzetti[kCalzetti] * coeffEBMV);
    }
  }
}

Float64 CSpectrumFluxCorrectionCalzetti::GetEbmvValue(Int32 k) const {
  Float64 coeffEBMV = m_EbmvCoeffStart + m_EbmvCoeffStep * (Float64)k;
  return coeffEBMV;
}

Int32 CSpectrumFluxCorrectionCalzetti::GetEbmvIndex(Float64 value) const {
  Int32 kEbmv = round((value - m_EbmvCoeffStart) / m_EbmvCoeffStep);
  return kEbmv;
}

Float64
CSpectrumFluxCorrectionCalzetti::GetDustCoeff(Int32 kDust,
                                              Float64 restLambda) const {
  if (kDust >= m_nEbmvCoeff)
    THROWG(ErrorCode::INTERNAL_ERROR, "ebmv index > nbEbmv ");
  Float64 coeffDust = 1.0;
  if (restLambda >= m_LambdaMin && restLambda < m_LambdaMax) {
    Int32 kCalzetti = Int32(round((restLambda - m_LambdaMin) / m_step));
    coeffDust = m_dataDustCoeff[kDust * m_dataCalzetti.size() + kCalzetti];
  }
  return coeffDust;
}
