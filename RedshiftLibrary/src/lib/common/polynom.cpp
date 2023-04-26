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

#include "RedshiftLibrary/common/polynom.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"

using namespace NSEpic;

TPolynomCoeffs::TPolynomCoeffs(const TFloat64List &coeffs) {
  if (coeffs.size() <= degree)
    THROWG(INTERNAL_ERROR,
           Formatter()
               << "input array too small to initialize a polynomial of degree "
               << degree);
  a0 = coeffs[0];
  a1 = coeffs[1];
  a2 = coeffs[2];
}

Float64 TPolynomCoeffs::getValue(Float64 x) const {
  Float64 val = a0;
  val += a1 * x;
  val += a2 * x * x;
  return val;
}

Float64 TPolynomCoeffs::getValueAndGrad(Float64 x, TFloat64List &grad) const {
  grad.resize(degree + 1);
  grad[0] = 1.0;
  grad[1] = x;
  grad[2] = x * x;
  return grad[0] * a0 + grad[1] * a1 + grad[2] * a2;
}

constexpr Int32 TPolynomCoeffs::degree;
