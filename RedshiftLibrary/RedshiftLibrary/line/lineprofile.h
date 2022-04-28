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
#include <cmath>
#include <string>
namespace NSEpic {
/**
 * struct that holds ASYMFIXED profile parameters
 */

enum TProfile {
  NONE,
  SYM,
  LOR,
  ASYM,
  ASYMFIT,
  // ASYMFIXED,//doesnt exist anymore since merged with ASYM
};
/**
 * \ingroup Redshift
 * Abstract class for different line profiles
 */
class CLineProfile {

public:
  CLineProfile(const Float64 nsigmasupport = N_SIGMA_SUPPORT);
  CLineProfile(const Float64 nsigmasupport = N_SIGMA_SUPPORT,
               const TProfile = NONE);
  virtual Float64 GetLineProfile(Float64 x, Float64 x0,
                                 Float64 sigma) const = 0;
  virtual Float64 GetLineFlux(Float64 A, Float64 sigma) const = 0;
  virtual Float64 GetLineProfileDerivZ(Float64 x, Float64 lambda0,
                                       Float64 redshift,
                                       Float64 sigma) const = 0;
  virtual Float64 GetLineProfileDerivSigma(Float64 x, Float64 x0,
                                           Float64 sigma) const = 0;
  virtual Float64 GetNSigmaSupport() const;

  const TProfile &GetName() const;
  virtual const TAsymParams GetAsymParams() const { return {NAN, NAN, NAN}; };
  virtual Float64 GetAsymDelta() const;
  virtual bool isAsymFit() const;
  virtual bool isAsymFixed() const;
  virtual void SetAsymParams(TAsymParams params){};
  virtual void resetAsymFitParams();
  CLineProfile(const CLineProfile &other) = default;
  CLineProfile(CLineProfile &&other) = default;
  CLineProfile &operator=(const CLineProfile &other) = default;
  CLineProfile &operator=(CLineProfile &&other) = default;

  virtual ~CLineProfile(){}; // to make sure derived objects are correctly
                             // deleted from a pointer to the base class
  // std::shared_ptr<CLineProfile> Clone () const {return
  // std::shared_ptr<CLineProfile>(CloneImplementation());}
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

inline const TProfile &CLineProfile::GetName() const { return m_name; }

// no need to define a constructor here
inline CLineProfile::CLineProfile(const Float64 nsigmasupport)
    : m_nsigmasupport(nsigmasupport), m_name(NONE) {}

inline CLineProfile::CLineProfile(const Float64 nsigmasupport,
                                  const TProfile name)
    : m_nsigmasupport(nsigmasupport), m_name(name) {}

inline Float64 CLineProfile::GetNSigmaSupport() const {
  return m_nsigmasupport;
}
inline Float64 CLineProfile::GetAsymDelta() const {
  return 0.; // default. Mainly used for asylfit/fixed
}
inline bool CLineProfile::isAsymFit() const {
  return 0; // default to no
}
inline bool CLineProfile::isAsymFixed() const {
  return 0; // default to no
}
inline void CLineProfile::resetAsymFitParams() { return; }

} // namespace NSEpic
#endif
