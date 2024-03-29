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
#ifndef _REDSHIFT_SPECTRUM_FLUXCORRECTIONCALZETTI_
#define _REDSHIFT_SPECTRUM_FLUXCORRECTIONCALZETTI_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/range.h"

namespace NSEpic {
typedef struct CalzettiCorrection {
  CalzettiCorrection(TFloat64List _lbda, TFloat64List _fluxcorr);

  TFloat64List lbda;
  TFloat64List fluxcorr;
  Float64 step;
} CalzettiCorrection;

/**
 * \ingroup Redshift
 */
class CSpectrumFluxCorrectionCalzetti {

public:
  CSpectrumFluxCorrectionCalzetti(CalzettiCorrection _calzettiCorr,
                                  Float64 ebmv_start, Float64 ebmv_step,
                                  Float64 ebmv_n);

  Float64 getLambdaMin() const { return m_LambdaMin; };
  Float64 getLambdaMax() const { return m_LambdaMax; };
  Int32 GetNPrecomputedEbmvCoeffs() const { return m_nEbmvCoeff; };

  Float64 GetEbmvValue(Int32 k) const;

  Float64 GetDustCoeff(Int32 kDust, Float64 restLambda) const;
  Int32 GetEbmvIndex(Float64 ebmv) const;

  TFloat64List m_dataCalzetti;

  Int32 m_nEbmvCoeff = 0;
  Float64 m_EbmvCoeffStep;
  Float64 m_EbmvCoeffStart;
  TFloat64List m_dataDustCoeff;

private:
  Float64 m_LambdaMin;
  Float64 m_LambdaMax;
  Float64 m_step;
};

} // namespace NSEpic

#endif
