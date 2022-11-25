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
#ifndef _REDSHIFT_COMMON_WEIGHTS_
#define _REDSHIFT_COMMON_WEIGHTS_

#include "RedshiftLibrary/common/datatypes.h"

namespace NSEpic {

/**
 * \ingroup Redshift
 * Wavelength interval objects.
 */
class CMask {

public:
  CMask() = default;
  explicit CMask(Int32 weightsCount);
  CMask(Int32 weightsCount, Int32 defaultValue)
      : m_Mask(weightsCount, defaultValue){};
  const Mask *GetMasks() const;
  CMask &operator&=(const CMask &other);
  Int32 GetMasksCount() const;
  Mask operator[](const Int32 i) const;
  Mask &operator[](const Int32 i);
  Float64 CompouteOverlapRate(const CMask &other) const;
  Float64 IntersectAndComputeOverlapRate(const CMask &other) const;

  bool IntersectWith(const CMask &other);
  Int32 GetMaskedSampleCount() const;
  Int32 GetUnMaskedSampleCount() const;
  void SetSize(Int32 s);

private:
  TMaskList m_Mask;
};

inline Mask CMask::operator[](const Int32 i) const { return m_Mask[i]; }

inline Mask &CMask::operator[](const Int32 i) { return m_Mask[i]; }

inline Int32 CMask::GetMasksCount() const { return m_Mask.size(); }

inline const Mask *CMask::GetMasks() const { return m_Mask.data(); }

inline Int32 CMask::GetMaskedSampleCount() const {
  return m_Mask.size() - GetUnMaskedSampleCount();
}

inline void CMask::SetSize(Int32 s) { m_Mask.resize(s); }

inline Int32 CMask::GetUnMaskedSampleCount() const {
  Int32 n = 0;
  for (Int32 i = 0; i < (Int32)m_Mask.size(); i++) {
    n += m_Mask[i];
  }
  return n;
}

} // namespace NSEpic

#endif
