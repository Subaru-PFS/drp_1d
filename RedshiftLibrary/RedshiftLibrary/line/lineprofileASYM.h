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
#ifndef _REDSHIFT_LINE_PROFILE_ASYM_
#define _REDSHIFT_LINE_PROFILE_ASYM_

#include <cmath>
#include <string>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/line/lineprofile.h"

namespace lineProfile_test {
class lineprofileASYM_test;
} // namespace lineProfile_test

namespace NSEpic {
/**
 * \ingroup Redshift
 */

class CLineProfileASYMFIT;
class CLineProfileASYM : public CLineProfile {
public:
  CLineProfileASYM(Float64 nsigmasupport = N_SIGMA_SUPPORT,
                   const TAsymParams &params = ASYM_DEFAULT_PARAMS,
                   const std::string &centeringMethod = "none");
  virtual ~CLineProfileASYM() = default;
  CLineProfileASYM(const CLineProfileASYM &other) = default;
  CLineProfileASYM(CLineProfileASYM &&other) = default;
  CLineProfileASYM &operator=(const CLineProfileASYM &other) = default;
  CLineProfileASYM &operator=(CLineProfileASYM &&other) = default;

  CLineProfileASYM(
      const CLineProfileASYMFIT &other); // ASYMFIT -> ASYM converter

  Float64 GetLineProfileVal(Float64 x, Float64 x0,
                            Float64 sigma) const override;
  Float64 GetLineFlux(Float64 x0, Float64 sigma,
                      Float64 A = 1.0) const override;
  Float64 GetLineProfileDerivX0(Float64 x, Float64 x0,
                                Float64 sigma) const override;
  Float64 GetLineProfileDerivSigma(Float64 x, Float64 x0,
                                   Float64 sigma) const override;
  Float64 GetNSigmaSupport() const override;
  TAsymParams GetAsymParams() const override;
  Float64 GetDelta() const override { return m_asym_delta; };
  virtual bool isAsym() const override { return true; };
  virtual bool isAsymFixed() const override { return true; };

  virtual const std::string &GetCenteringMethod() const {
    return m_centeringMethod;
  };

  std::unique_ptr<CLineProfile> cloneToASYMFIT() const;

private:
  friend class lineProfile_test::lineprofileASYM_test;
  virtual CLineProfile *CloneImplementation() const override {
    return new CLineProfileASYM(*this);
  }
  Float64 GetXSurc(Float64 xc, Float64 &sigma, Float64 &xsurc) const;

protected:
  CLineProfileASYM(
      TProfile pltype, Float64 nsigmasupport = N_SIGMA_SUPPORT,
      const TAsymParams &params = ASYMF_DEFAULT_PARAMS,
      const std::string &centeringMethod = "mean"); // used by asymfit ctor

  bool isValid() const;
  Float64 m_asym_sigma_coeff = ASYM_DEFAULT_PARAMS.sigma;
  Float64 m_asym_alpha = ASYM_DEFAULT_PARAMS.alpha;
  Float64 m_asym_delta = ASYM_DEFAULT_PARAMS.delta;
  std::string m_centeringMethod = "none";
  Float64 m_constSigma = 1; // vs 2.5 for AsymFit and AsymFixed
};
} // namespace NSEpic
#endif
