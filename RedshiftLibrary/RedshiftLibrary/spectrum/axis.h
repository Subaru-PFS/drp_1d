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
#ifndef _REDSHIFT_SPECTRUM_AXIS_
#define _REDSHIFT_SPECTRUM_AXIS_

#include "RedshiftLibrary/common/datatypes.h"

#include <vector>

namespace NSEpic {

/**
 * \ingroup Redshift
 */
class CSpectrumAxis {

public:
  CSpectrumAxis() = default;
  CSpectrumAxis(const CSpectrumAxis &other) = default;
  CSpectrumAxis(CSpectrumAxis &&other) = default;
  explicit CSpectrumAxis(Int32 n, Float64 value = 0.0) : m_Samples(n, value){};
  CSpectrumAxis(const Float64 *samples, Int32 n);
  CSpectrumAxis(const TFloat64List &samples) : m_Samples(samples){};
  CSpectrumAxis(TFloat64List &&samples) : m_Samples(std::move(samples)){};

  virtual ~CSpectrumAxis() = default;
  CSpectrumAxis &operator=(const CSpectrumAxis &other) = default;
  CSpectrumAxis &operator=(CSpectrumAxis &&other) = default;
  virtual CSpectrumAxis &operator*=(const Float64 op);
  virtual CSpectrumAxis &operator/=(const Float64 op);
  Float64 &operator[](const Int32 i);
  const Float64 &operator[](const Int32 i) const;
  CSpectrumAxis MaskAxis(const TFloat64List &mask) const;
  static const TFloat64List maskVector(const TFloat64List &mask,
                                       const TFloat64List &inputVector);

  const Float64 *GetSamples() const;
  Float64 *GetSamples();
  const TAxisSampleList &GetSamplesVector() const;
  TAxisSampleList &GetSamplesVector();
  Int32 GetSamplesCount() const;
  Int32 GetSamplesCount();
  virtual void SetSize(Int32 s);
  void clear();
  CSpectrumAxis extract(Int32 startIdx, Int32 endIdx) const;
  bool isEmpty() const;

protected:
  TAxisSampleList m_Samples;
  virtual void resetAxisProperties(){}; // by default it does nothing
};

inline Float64 &CSpectrumAxis::operator[](const Int32 i) {
  resetAxisProperties();
  return m_Samples[i];
}

inline const Float64 &CSpectrumAxis::operator[](const Int32 i) const {
  return m_Samples[i];
}

inline Int32 CSpectrumAxis::GetSamplesCount() const { return m_Samples.size(); }

inline Float64 *CSpectrumAxis::GetSamples() {
  resetAxisProperties();
  return m_Samples.data();
}

inline const Float64 *CSpectrumAxis::GetSamples() const {
  return m_Samples.data();
}

inline TAxisSampleList &CSpectrumAxis::GetSamplesVector() {
  resetAxisProperties();
  return m_Samples;
}

inline const TAxisSampleList &CSpectrumAxis::GetSamplesVector() const {
  return m_Samples;
}
inline bool CSpectrumAxis::isEmpty() const { return m_Samples.size() == 0; }

inline CSpectrumAxis CSpectrumAxis::extract(Int32 startIdx,
                                            Int32 endIdx) const {
  if (!m_Samples.size())
    return CSpectrumAxis();
  return CSpectrumAxis(TFloat64List(m_Samples.begin() + startIdx,
                                    m_Samples.begin() + endIdx + 1));
}
} // namespace NSEpic
#endif
