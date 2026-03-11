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
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/vectorOperations.h"
#include <algorithm>
#include <functional>

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
  CSpectrumAxis(const Float64 *samples, Int32 n)
      : m_Samples(samples, samples + n){};
  explicit CSpectrumAxis(const TFloat64List &samples) : m_Samples(samples){};
  explicit CSpectrumAxis(TFloat64List &&samples)
      : m_Samples(std::move(samples)){};

  virtual ~CSpectrumAxis() = default;
  CSpectrumAxis &operator=(const CSpectrumAxis &other) = default;
  CSpectrumAxis &operator=(CSpectrumAxis &&other) = default;
  virtual CSpectrumAxis &operator*=(Float64 op);
  virtual CSpectrumAxis &operator/=(Float64 op);
  virtual CSpectrumAxis &operator+=(CSpectrumAxis const &other);
  virtual CSpectrumAxis &operator-=(CSpectrumAxis const &other);
  Float64 operator[](Int32 i) const;
  CSpectrumAxis MaskAxis(const TMaskList &mask) const;
  // virtual CSpectrumAxis operator+(CSpectrumAxis other) const;
  // virtual CSpectrumAxis operator-(CSpectrumAxis other) const;
  // Hidden friend symmetric operators (ie non-member, but here for ADL )
  friend CSpectrumAxis operator*(CSpectrumAxis axis, Float64 op);
  friend CSpectrumAxis operator*(Float64 op, CSpectrumAxis axis);
  friend CSpectrumAxis operator+(CSpectrumAxis const &axis1,
                                 CSpectrumAxis const &axis2);
  friend CSpectrumAxis operator-(CSpectrumAxis const &axis1,
                                 CSpectrumAxis const &axis2);
  // member assymetric operator
  CSpectrumAxis operator/(Float64 op) const;

  const Float64 *GetSamples() const;
  const TAxisSampleList &GetSamplesVector() const &;
  TAxisSampleList &&GetSamplesVector() &&;
  virtual void setSamplesVector(TAxisSampleList axisList);
  Int32 GetSamplesCount() const;
  virtual void resize(Int32 s, Float64 valueDef = 0.0);
  virtual void clear();
  void Invert();
  void Negate();

  CSpectrumAxis extract(Int32 startIdx, Int32 endIdx) const;
  bool isEmpty() const;

protected:
  TAxisSampleList m_Samples;
  virtual void resetAxisProperties(){}; // by default it does nothing
};

inline Float64 CSpectrumAxis::operator[](Int32 i) const { return m_Samples[i]; }

inline CSpectrumAxis &CSpectrumAxis::operator*=(Float64 op) {
  resetAxisProperties();
  std::transform(m_Samples.cbegin(), m_Samples.cend(), m_Samples.begin(),
                 [op](Float64 sample) { return sample * op; });
  return *this;
}

inline CSpectrumAxis &CSpectrumAxis::operator/=(Float64 op) {
  operator*=(1 / op);
  return *this;
}

inline CSpectrumAxis &CSpectrumAxis::operator+=(CSpectrumAxis const &other) {
  if (other.GetSamplesCount() != GetSamplesCount())
    THROWG(ErrorCode::INTERNAL_ERROR, "Cannot sum axis of different sizes");
  resetAxisProperties();
  std::transform(m_Samples.cbegin(), m_Samples.cend(), other.m_Samples.cbegin(),
                 m_Samples.begin(), std::plus<>());
  return *this;
}

inline CSpectrumAxis &CSpectrumAxis::operator-=(CSpectrumAxis const &other) {
  if (other.GetSamplesCount() != GetSamplesCount())
    THROWG(ErrorCode::INTERNAL_ERROR,
           "Cannot subtract axis of different sizes");
  resetAxisProperties();
  std::transform(m_Samples.cbegin(), m_Samples.cend(), other.m_Samples.cbegin(),
                 m_Samples.begin(), std::minus<>());
  return *this;
}

inline CSpectrumAxis operator*(CSpectrumAxis axis, Float64 op) {
  CSpectrumAxis multipliedAxis(std::move(axis));
  multipliedAxis *= op;
  return multipliedAxis;
}

inline CSpectrumAxis operator*(Float64 op, CSpectrumAxis axis) {
  return std::move(axis) * op;
}

inline CSpectrumAxis operator+(CSpectrumAxis const &axis1,
                               CSpectrumAxis const &axis2) {
  CSpectrumAxis sumaxis(axis1);
  sumaxis += axis2;
  return sumaxis;
}

inline CSpectrumAxis operator-(CSpectrumAxis const &axis1,
                               CSpectrumAxis const &axis2) {
  CSpectrumAxis diffaxis(axis1);
  diffaxis -= axis2;
  return diffaxis;
}

// inline CSpectrumAxis CSpectrumAxis::operator+(CSpectrumAxis other) const {
//   CSpectrumAxis sumaxis(std::move(other));
//   sumaxis += *this;
//   return sumaxis;
// }

// inline CSpectrumAxis CSpectrumAxis::operator-(CSpectrumAxis other) const {
//   CSpectrumAxis diffaxis(std::move(other));
//   diffaxis -= *this;
//   return diffaxis;
// }

inline CSpectrumAxis CSpectrumAxis::operator/(Float64 op) const {
  CSpectrumAxis dividedAxis = *this;
  dividedAxis /= op;
  return dividedAxis;
}

inline Int32 CSpectrumAxis::GetSamplesCount() const { return m_Samples.size(); }

inline const Float64 *CSpectrumAxis::GetSamples() const {
  return m_Samples.data();
}

inline void CSpectrumAxis::setSamplesVector(TAxisSampleList axisList) {
  resetAxisProperties();
  m_Samples = std::move(axisList);
}

inline const TAxisSampleList &CSpectrumAxis::GetSamplesVector() const & {
  return m_Samples;
}

inline TAxisSampleList &&CSpectrumAxis::GetSamplesVector() && {
  return std::move(m_Samples);
}

inline void CSpectrumAxis::resize(Int32 s, Float64 valueDef) {
  m_Samples.resize(s, valueDef);
}

inline void CSpectrumAxis::clear() {
  resetAxisProperties();
  m_Samples.clear();
}

inline bool CSpectrumAxis::isEmpty() const { return m_Samples.size() == 0; }

inline void CSpectrumAxis::Invert() {
  std::transform(m_Samples.begin(), m_Samples.end(), m_Samples.begin(),
                 [](Float64 val) { return 1 / val; });
}

inline void CSpectrumAxis::Negate() {
  std::transform(m_Samples.begin(), m_Samples.end(), m_Samples.begin(),
                 std::negate<Float64>());
}

inline CSpectrumAxis CSpectrumAxis::extract(Int32 startIdx,
                                            Int32 endIdx) const {
  if (!m_Samples.size())
    return CSpectrumAxis();
  if (startIdx < 0 || startIdx >= GetSamplesCount())
    THROWG(ErrorCode::INTERNAL_ERROR, "startIdx out of bounds");
  if (endIdx < 0 || endIdx >= GetSamplesCount())
    THROWG(ErrorCode::INTERNAL_ERROR, "endIdx out of bounds");
  return CSpectrumAxis(TFloat64List(m_Samples.begin() + startIdx,
                                    m_Samples.begin() + endIdx + 1));
}

inline CSpectrumAxis
CSpectrumAxis::MaskAxis(const TMaskList &mask) const // mask is 0. or 1.
{
  return CSpectrumAxis(NSVectorOp::maskVector<Float64>(mask, m_Samples));
}

} // namespace NSEpic
#endif
