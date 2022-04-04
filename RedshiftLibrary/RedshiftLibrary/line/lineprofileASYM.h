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
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/line/lineprofile.h"
#include <math.h>
#include <string>
namespace NSEpic {
/**
 * \ingroup Redshift
 */
class CLineProfileASYM : public CLineProfile {
public:
  CLineProfileASYM(const Float64 nsigmasupport = 8.0,
                   const TAsymParams params = {1., 4.5, 0.},
                   const std::string centeringMethod = "none");

  Float64 GetLineProfile(Float64 x, Float64 x0,
                         const Float64 sigma) const override;
  Float64 GetLineFlux(Float64 A, const Float64 sigma) const override;
  Float64 GetLineProfileDerivZ(Float64 x, Float64 x0, Float64 redshift,
                               const Float64 sigma) const override;
  Float64 GetLineProfileDerivSigma(Float64 x, Float64 x0,
                                   const Float64 sigma) const override;
  Float64 GetNSigmaSupport() const override;

  Float64 GetAsymDelta() const override;
  const TAsymParams GetAsymParams() const override;
  virtual bool isAsymFixed() const override;
  virtual bool isAsymFit() const override;

  virtual ~CLineProfileASYM() = default;
  CLineProfileASYM(const CLineProfileASYM &other) = default;
  CLineProfileASYM(CLineProfileASYM &&other) = default;
  CLineProfileASYM &operator=(const CLineProfileASYM &other) = default;
  CLineProfileASYM &operator=(CLineProfileASYM &&other) = default;

private:
  virtual CLineProfile *CloneImplementation() const override {
    return new CLineProfileASYM(*this);
  }
  Float64 GetXSurc(Float64 xc, Float64 &sigma, Float64 &xsurc) const;

protected:
  bool isValid() const;
  Float64 m_asym_sigma_coeff = 1.0; // vs 2. for asymFit/Fixed
  Float64 m_asym_alpha = 4.5;
  Float64 m_asym_delta = 0.;
  std::string m_centeringMethod = "none";
  Float64 m_constSigma = 1; // vs 2.5 for AsymFit and AsymFixed
  CLineProfileASYM(
      const TProfile pltype, const Float64 nsigmasupport = 8.0,
      const TAsymParams params = {2., 2.5, 0.},
      const std::string centeringMethod = "mean"); // mainly called by asymfit
};
} // namespace NSEpic
#endif