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
#include "RedshiftLibrary/line/lineprofileASYM.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/line/lineprofileASYMFIT.h"
#include "RedshiftLibrary/log/log.h"
using namespace NSEpic;
using namespace std;

CLineProfileASYM::CLineProfileASYM(Float64 nsigmasupport,
                                   const TAsymParams &params,
                                   const std::string &centeringMethod)
    : CLineProfileASYM(ASYM, nsigmasupport, params, centeringMethod) {}

CLineProfileASYM::CLineProfileASYM(TProfile pltype, Float64 nsigmasupport,
                                   const TAsymParams &params,
                                   const std::string &centeringMethod)
    : CLineProfile(nsigmasupport, pltype), m_asym_sigma_coeff(params.sigma),
      m_asym_alpha(params.alpha), m_asym_delta(params.delta),
      m_centeringMethod(centeringMethod) {
  if (m_centeringMethod == "mean" || m_centeringMethod == "mode")
    m_constSigma = 2.5;
}

// asymfit -> asymfixed converter
CLineProfileASYM::CLineProfileASYM(const CLineProfileASYMFIT &other)
    : CLineProfileASYM(other.m_nsigmasupport, other.GetAsymParams(),
                       other.m_centeringMethod) {}

std::unique_ptr<CLineProfile> CLineProfileASYM::cloneToASYMFIT() const {
  return std::unique_ptr<CLineProfile>(new CLineProfileASYMFIT(*this));
}

Float64 CLineProfileASYM::GetXSurc(Float64 xc, Float64 &sigma,
                                   Float64 &xsurc) const {
  sigma = sigma * m_asym_sigma_coeff;
  Float64 delta = m_asym_alpha / std::sqrt(1. + m_asym_alpha * m_asym_alpha);
  Float64 muz = delta * sqrt(2. / M_PI);

  Float64 m0 = 0;
  if (m_centeringMethod == "none") { // ASYM
    m0 = 0;                          // i.e., no centering;
    delta = 0.;                      // reset it cause no centering here

    if (m_asym_delta) // temporary check
      THROWG(
          ErrorCode::INTERNAL_ERROR,
          "Problem in configuring ASYM lineprofile: asym_delta should be null");
  } else if (m_centeringMethod == "mean") // ASYMFIT, ASYMFIXED
  { // correction in order to have the line shifted on the mean: from
    // https://en.wikipedia.org/wiki/Skew_normal_distribution
    m0 = muz;
  } else if (m_centeringMethod ==
             "mode") { // correction in order to have the line shifted on the
                       // mode: from
                       // https://en.wikipedia.org/wiki/Skew_normal_distribution
    Float64 sigmaz = std::sqrt(1 - muz * muz);
    Float64 gamma1 = ((4 - M_PI) / 2.0) * pow(delta * std::sqrt(2 / M_PI), 3.) /
                     pow(1 - 2 * delta * delta / M_PI, 3. / 2.);
    m0 = muz - gamma1 * sigmaz / 2.0 - 0.5 * exp(-2 * M_PI / m_asym_alpha);
  }

  xc = xc + sigma * m0;

  Float64 xcd = xc + m_asym_delta;
  xsurc = xcd / sigma;
  return m0;
}

Float64 CLineProfileASYM::GetLineProfileVal(Float64 x, Float64 x0,
                                            Float64 sigma) const {
  if (!isValid()) {
    THROWG(ErrorCode::INTERNAL_ERROR, "LineProfile is not valid");
  }
  Float64 xsurc;
  GetXSurc(x - x0, sigma, xsurc);

  Float64 val =
      exp(-0.5 * xsurc * xsurc) * (1.0 + erf(m_asym_alpha / sqrt(2.0) * xsurc));
  return val;
}

Float64 CLineProfileASYM::GetNSigmaSupport() const {
  return m_nsigmasupport * m_asym_sigma_coeff * m_constSigma;
}

Float64 CLineProfileASYM::GetLineFlux(Float64 x0, Float64 sigma,
                                      Float64 A) const {
  return A * sigma * m_asym_sigma_coeff * sqrt(2 * M_PI);
}

Float64 CLineProfileASYM::GetLineProfileDerivX0(Float64 x, Float64 x0,
                                                Float64 sigma) const {
  if (!isValid()) {
    THROWG(ErrorCode::INTERNAL_ERROR, "LineProfile is not valid");
  }
  Float64 xsurc;
  Float64 alpha2 = m_asym_alpha * m_asym_alpha;
  Float64 xc = x - x0;

  GetXSurc(xc, sigma, xsurc);

  Float64 xsurc2 = xsurc * xsurc;
  Float64 val =
      xsurc / sigma * exp(-0.5 * xsurc2) *
          (1.0 + erf(m_asym_alpha / sqrt(2.0) * xsurc)) -
      m_asym_alpha / sqrt(2 * M_PI) / sigma * exp(-(1 + alpha2) / 2 * xsurc2);
  return val;
}

Float64 CLineProfileASYM::GetLineProfileDerivSigma(Float64 x, Float64 x0,
                                                   Float64 sigma) const {
  if (!isValid()) {
    THROWG(ErrorCode::INTERNAL_ERROR, "LineProfile is not valid");
  }
  Float64 alpha2 = m_asym_alpha * m_asym_alpha;

  Float64 xsurc;
  Float64 muz = GetXSurc(x - x0, sigma, xsurc);
  Float64 xsurc2 = xsurc * xsurc;

  Float64 valSym = exp(-0.5 * xsurc2);
  Float64 valAsym = (1.0 + erf(m_asym_alpha / sqrt(2.0) * xsurc));

  Float64 valSymD, valAsymD;
  if (m_centeringMethod == "none") {
    // temporary: double check muz == 0
    if (muz)
      THROWG(ErrorCode::INTERNAL_ERROR,
             "Problem: mu is not null for ASYM profile!");
  }
  valSymD = m_asym_sigma_coeff * (xsurc2 / sigma - muz * xsurc / sigma) *
            exp(-0.5 * xsurc2);
  valAsymD = m_asym_sigma_coeff * m_asym_alpha * sqrt(2) /
             (sigma * sqrt(M_PI)) * (muz - xsurc) * exp(-0.5 * xsurc2 * alpha2);

  Float64 val = valSym * valAsymD + valSymD * valAsym;
  return val;
}

TAsymParams CLineProfileASYM::GetAsymParams() const {
  return {m_asym_sigma_coeff, m_asym_alpha, m_asym_delta};
}

bool CLineProfileASYM::isValid() const {
  if (std::isnan(m_asym_sigma_coeff) || std::isnan(m_asym_alpha) ||
      std::isnan(m_asym_alpha)) {
    return 0;
  }
  return 1;
}
