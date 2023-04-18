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
#ifndef _REDSHIFT_LINE_PROFILE_
#define _REDSHIFT_LINE_PROFILE_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h"
#include <cmath>
#include <string>
namespace NSEpic {
class CSpectrumFluxCorrectionMeiksin;
/**
 * \ingroup Redshift
 * Abstract class for different line profiles
 */

/**
 * struct that holds ASYMFIXED profile parameters
 */
typedef struct TAsymParams {

  TAsymParams(Float64 sigma, Float64 alpha, Float64 delta)
      : sigma(sigma), alpha(alpha), delta(delta){};

  TAsymParams() = default;
  Float64 sigma = NAN, alpha = NAN, delta = NAN;
} TAsymParams;

// lineprofileAsym
static const TAsymParams ASYM_DEFAULT_PARAMS{1.0, 4.5, 0.};
static const Float64 ASYM_DEFAULT_CONSTSIGMA = 1.;
// lineprofileAsym Fit anf Fixed
static const TAsymParams ASYMF_DEFAULT_PARAMS{2.0, 2.0, 0.};
static const Float64 ASYMF_DEFAULT_CONSTSIGMA = 2.5;

typedef struct TSymIgmParams {
  TSymIgmParams() = default;
  TSymIgmParams(Int32 igmidx, Float64 redshift)
      : m_igmidx(igmidx), m_redshift(redshift){};

  Int32 m_igmidx = undefIdx;
  Float64 m_redshift = NAN;
} TSymIgmParams;

enum TProfile { NONE, SYM, LOR, ASYM, ASYMFIT, SYMIGM };

class CLineProfile {

public:
  CLineProfile(Float64 nsigmasupport = N_SIGMA_SUPPORT, TProfile name = NONE)
      : m_nsigmasupport(nsigmasupport), m_name(name){};

  virtual ~CLineProfile() = default;

  CLineProfile(const CLineProfile &other) = default;
  CLineProfile(CLineProfile &&other) = default;
  CLineProfile &operator=(const CLineProfile &other) = default;
  CLineProfile &operator=(CLineProfile &&other) = default;

  virtual Float64 GetLineProfileVal(Float64 x, Float64 x0,
                                    Float64 sigma) const = 0;

  virtual Float64 GetLineFlux(Float64 x0, Float64 sigma,
                              Float64 A = 1.0) const = 0;

  virtual Float64 GetLineProfileDerivX0(Float64 x, Float64 x0,
                                        Float64 sigma) const = 0;

  virtual Float64 GetLineProfileDerivSigma(Float64 x, Float64 x0,
                                           Float64 sigma) const = 0;

  virtual Float64 GetNSigmaSupport() const { return m_nsigmasupport; };

  const TProfile &GetName() const { return m_name; };
  virtual TAsymParams GetAsymParams() const { return TAsymParams(); };
  virtual Float64 GetDelta() const { return 0.; };
  virtual TSymIgmParams GetSymIgmParams() const { return TSymIgmParams(); };
  virtual bool isAsym() const { return false; };
  virtual bool isAsymFit() const { return false; };
  virtual bool isAsymFixed() const { return false; };
  virtual bool isSymIgm() const { return false; };
  virtual bool isSymIgmFit() const { return false; };
  virtual void SetSymIgmFit(bool val = true){};
  virtual void SetSymIgmFixed(){};
  virtual void SetAsymParams(const TAsymParams &params){};
  virtual void SetSymIgmParams(const TSymIgmParams &params){};
  virtual void resetParams(){};

  virtual Int32 getIGMIdxCount() const;

  std::unique_ptr<CLineProfile> Clone() const {
    return std::unique_ptr<CLineProfile>(CloneImplementation());
  }

private:
  virtual CLineProfile *CloneImplementation() const = 0;

protected:
  Float64 m_nsigmasupport;
  const TProfile m_name; // hack to avoid using dynamic casting
};

typedef std::unique_ptr<CLineProfile> CLineProfile_ptr;
typedef std::vector<CLineProfile_ptr> TProfileList;

inline Int32 CLineProfile::getIGMIdxCount() const {
  THROWG(INTERNAL_ERROR, "getIGMIdxCount is not defined "
                         "for non-SYMIGM lineprofile");
}
} // namespace NSEpic
#endif
