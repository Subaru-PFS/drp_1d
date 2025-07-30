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
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include <tuple>

using namespace NSEpic;

CPolynomCoeffs::CPolynomCoeffs(const TFloat64List &coeffs) {
  if (coeffs.size() <= degree)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter()
               << "input array too small to initialize a polynomial of degree "
               << degree);
  m_a0 = coeffs[0];
  m_a1 = coeffs[1];
  m_a2 = coeffs[2];
}

Float64 CPolynomCoeffs::getValue(Float64 x) const {

  Float64 val = m_a2 * x;
  val += m_a1;
  val *= x;
  val += m_a0;

  return val;
}

TFloat64List CPolynomCoeffs::getPowers(Float64 x) const {
  TFloat64List grad(degree + 1);
  grad[0] = 1.0;
  grad[1] = x;
  grad[2] = x * x;
  return grad;
}

Float64 CPolynomCoeffs::getVariance(Float64 x) const {
  if ((m_covar.array() == 0.0).all())
    return NAN;
  Eigen::Vector3d powerx(getPowers(x).data());
  auto const var = powerx.transpose() * m_covar * powerx;
  return var;
}

CPolynomCoeffs CPolynomCoeffs::operator*(Float64 factor) const {
  CPolynomCoeffs poly;
  poly.m_a0 = m_a0 * factor;
  poly.m_a1 = m_a1 * factor;
  poly.m_a2 = m_a2 * factor;

  poly.m_covar = m_covar * (factor * factor);
  return poly;
}

constexpr Int32 CPolynomCoeffs::degree;

CPolynomCoeffsNormalized::CPolynomCoeffsNormalized(Float64 x0_, Float64 scale_)
    : m_x0red(-x0_ / scale_), m_scale(scale_),
      m_convCoeff({{1.0, m_x0red, m_x0red * m_x0red},
                   {0.0, 1 / m_scale, 2 * m_x0red / m_scale},
                   {0.0, 0.0, 1 / (m_scale * m_scale)}}),
      m_convCoeffInv(
          {{1.0, -m_scale * m_x0red, m_scale * m_scale * m_x0red * m_x0red},
           {0.0, m_scale, -2 * m_scale * m_scale * m_x0red},
           {0.0, 0.0, m_scale * m_scale}}){};

CPolynomCoeffs CPolynomCoeffsNormalized::getPolynomCoeffs() const {
  CPolynomCoeffs poly;
  auto const coeffs = Eigen::Vector3d{m_a0, m_a1, m_a2};
  auto const coeffs_out = m_convCoeff * coeffs;
  poly.m_a0 = coeffs_out(0);
  poly.m_a1 = coeffs_out(1);
  poly.m_a2 = coeffs_out(2);

  if ((m_covar.array() != 0.0).any())
    poly.m_covar = m_convCoeff * m_covar * m_convCoeff.transpose();
  return poly;
}

void CPolynomCoeffsNormalized::setFromPolynomCoeffs(
    const CPolynomCoeffs &poly) {
  auto const coeffs_in = Eigen::Vector3d{poly.m_a0, poly.m_a1, poly.m_a2};
  auto const coeffs = m_convCoeffInv * coeffs_in;
  m_a0 = coeffs(0);
  m_a1 = coeffs(1);
  m_a2 = coeffs(2);

  if ((poly.m_covar.array() != 0.0).any())
    m_covar = m_convCoeffInv * poly.m_covar * m_convCoeffInv.transpose();
}

std::pair<Float64, TFloat64List>
CPolynomCoeffsNormalized::getValueAndGradiant(Float64 x) const {
  Float64 xred = getXred(x);
  Float64 val = CPolynomCoeffs::getValue(xred);
  TFloat64List grad = CPolynomCoeffs::getCoeffGradiant(xred);
  return std::make_pair(val, grad);
}