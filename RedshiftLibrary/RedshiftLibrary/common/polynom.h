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
#ifndef _REDSHIFT_COMMON_POLYNOM_
#define _REDSHIFT_COMMON_POLYNOM_

#include "RedshiftLibrary/common/datatypes.h"

namespace NSEpic {

struct TPolynomCoeffs {
  TPolynomCoeffs() = default;
  TPolynomCoeffs(Float64 a0_, Float64 a1_ = 0.0, Float64 a2_ = 0.0)
      : a0(a0_), a1(a1_), a2(a2_){};
  TPolynomCoeffs(const TFloat64List &coeffs);

  Float64 getValue(Float64 x) const;

  Float64 getValueAndGrad(Float64 x, TFloat64List &grad) const;

  static constexpr Int32 degree = 2;

  Float64 a0 = NAN;
  Float64 a1 = NAN;
  Float64 a2 = NAN;
};

class CPolynomCoeffsNormalized {
public:
  CPolynomCoeffsNormalized() = default;
  CPolynomCoeffsNormalized(Float64 x0_, Float64 scale_ = 1.0)
      : x0red(-x0_ / scale_), scale(scale_){};

  void getCoeffs(Float64 &a0_, Float64 &a1_, Float64 &a2_) const;
  void setCoeffs(Float64 a0_, Float64 a1_, Float64 a2_);
  Float64 getValue(Float64 x) const;
  Float64 getValueAndGrad(Float64 x, TFloat64List &grad) const;

  static constexpr Int32 degree = 2;

  Float64 a0 = NAN;
  Float64 a1 = NAN;
  Float64 a2 = NAN;

private:
  Float64 x0red = 0.0;
  Float64 scale = 1.0;
};

} // namespace NSEpic
#endif
