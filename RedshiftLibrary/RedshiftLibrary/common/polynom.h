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

#include <cmath>

#include <Eigen/Core>
#include <tuple>
#include <vector>

#include "RedshiftLibrary/common/datatypes.h"

namespace NSEpic {

class CPolynomCoeffs {
public:
  CPolynomCoeffs() = default;
  CPolynomCoeffs(Float64 a0_, Float64 a1_ = 0.0, Float64 a2_ = 0.0,
                 Eigen::Matrix3d covar_ = Eigen::Matrix3d::Zero())
      : m_a0(a0_), m_a1(a1_), m_a2(a2_), m_covar(covar_){};
  CPolynomCoeffs(const TFloat64List &coeffs);

  virtual Float64 getValue(Float64 x) const;
  TFloat64List getPowers(Float64 x) const;
  virtual TFloat64List getCoeffGradiant(Float64 x) const {
    return getPowers(x);
  };
  virtual Float64 getVariance(Float64 x) const;
  CPolynomCoeffs operator*(Float64 factor) const;

  static constexpr Int32 degree = 2;

  Float64 m_a0 = NAN;
  Float64 m_a1 = NAN;
  Float64 m_a2 = NAN;

  Eigen::Matrix3d m_covar = Eigen::Matrix3d::Zero();
};

class CPolynomCoeffsNormalized : public CPolynomCoeffs {
public:
  CPolynomCoeffsNormalized() = default;
  CPolynomCoeffsNormalized(Float64 x0_, Float64 scale_ = 1.0);
  CPolynomCoeffs getPolynomCoeffs() const;
  void setFromPolynomCoeffs(const CPolynomCoeffs &);

  Float64 getValue(Float64 x) const override {
    return CPolynomCoeffs::getValue(getXred(x));
  };

  TFloat64List getCoeffGradiant(Float64 x) const override {
    return getPowers(getXred(x));
  };
  std::pair<Float64, TFloat64List> getValueAndGradiant(Float64 x) const;
  Float64 getVariance(Float64 x) const override {
    return CPolynomCoeffs::getVariance(getXred(x));
  };

private:
  Float64 getXred(Float64 x) const { return x / m_scale + m_x0red; };

  Float64 m_x0red = 0.0;
  Float64 m_scale = 1.0;
  Eigen::Matrix3d m_convCoeff = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d m_convCoeffInv = Eigen::Matrix3d::Zero();
};

} // namespace NSEpic
#endif
