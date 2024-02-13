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
#ifndef _REDSHIFT_COMMON_RANGE_
#define _REDSHIFT_COMMON_RANGE_

#include <cmath>
#include <vector>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/log/log.h"

namespace NSEpic {

/**
 * \ingroup Redshift
 * Templated Range manipulation class
 */
template <typename T> class CRange {

public:
  CRange() : m_Begin(), m_End() {}

  CRange(const T begin, const T end) : m_Begin(begin), m_End(end) {}

  CRange(const std::vector<T> &v)
      : m_Begin(v.empty() ? T() : v.front()),
        m_End(v.empty() ? T() : v.back()) {}

  ~CRange() {}

  bool GetIsEmpty() const { return m_Begin == m_End; }

  void Set(const T &b, const T &e) {
    m_Begin = b;
    m_End = e;
  }

  void SetBegin(const T &v) { m_Begin = v; }

  void SetEnd(const T &v) { m_End = v; }

  CRange<T> operator-(const T &t) const {
    return CRange<T>(m_Begin - t, m_End - t);
  }
  CRange<T> operator+(const T &t) const {
    return CRange<T>(m_Begin + t, m_End + t);
  }
  CRange<T> operator*(const T &t) const {
    return CRange<T>(m_Begin * t, m_End * t);
  }
  CRange<T> operator/(const T &t) const {
    return CRange<T>(m_Begin / t, m_End / t);
  }

  CRange<T> operator+(const CRange<T> &rhs) const {
    return CRange<T>(m_Begin + rhs.m_Begin, m_End + rhs.m_End);
  }

  bool operator==(const CRange<T> &rhs) const {
    if (m_Begin == rhs.m_Begin && m_End == rhs.m_End)
      return true;
    return false;
  }

  bool operator!=(const CRange<T> &rhs) const {
    if (m_Begin != rhs.m_Begin || m_End != rhs.m_End)
      return true;
    return false;
  }

  bool operator<(const CRange<T> &rhs) const {
    if (m_Begin == rhs.m_Begin)
      return (m_End < rhs.m_End);
    return (m_Begin < rhs.m_Begin);
  }

  Float64 Interpolate(Float64 t) const {
    return m_Begin + (m_End - m_Begin) * t;
  }

  const T &GetBegin() const { return m_Begin; }

  const T &GetEnd() const { return m_End; }

  T GetLength() const { return m_End - m_Begin; }

  bool isSameSign(T offset) const {
    return (m_Begin + offset) * (m_End + offset) >= 0;
  }

  T Clamp(T v) const { return v < m_Begin ? m_Begin : (m_End < v ? m_End : v); }

  std::vector<T> SpreadOver_backward(T delta) const {
    if (GetIsEmpty() || delta == 0.0 || GetLength() < delta)
      return {m_End};

    Int32 count = GetLength() / delta + 1;

    std::vector<T> v(count);
    for (Int32 i = 0; i < count; i++)
      v[i] = m_End - delta * (count - 1 - i);

    return v;
  }

  std::vector<T> SpreadOver(T delta, bool backward = false) const {
    if (backward)
      return SpreadOver_backward(delta);

    if (GetIsEmpty() || delta == 0.0 || GetLength() < delta)
      return {m_Begin};

    Int32 count = GetLength() / delta + 1;

    std::vector<T> v(count);
    for (Int32 i = 0; i < count; i++)
      v[i] = m_Begin + delta * i;

    return v;
  }

  std::vector<T> SpreadOverEpsilon(T delta, bool backward = false,
                                   Float64 epsilon = 1e-8) const {
    static_assert(std::is_same<T, Float64>::value,
                  "not implemented"); // compile time check
    CRange<T> range_epsilon =
        backward ? CRange<T>{-epsilon, 0.0} : CRange<T>{0.0, epsilon};
    CRange<T> range_ = *this + range_epsilon;
    auto ret = range_.SpreadOver(delta, backward);
    ret.front() = std::max(ret.front(), m_Begin);
    ret.back() = std::min(ret.back(), m_End);
    return ret;
  }

  std::vector<T> SpreadOverLog_backward(T delta, T offset = 0.) const {
    static_assert(std::is_same<T, Float64>::value,
                  "not implemented"); // compile time check
    if (!isSameSign(offset))
      THROWG(INTERNAL_ERROR, "borders should be of same sign");
    if (GetIsEmpty() || delta == 0.0 ||
        GetLength() <
            (GetBegin() + offset) * exp(delta) - (GetBegin() + offset))
      return {m_End};

    T x = m_End + offset;
    T edelta = exp(-delta);
    Int32 count =
        (log(GetEnd() + offset) - log(GetBegin() + offset)) / delta + 1;

    std::vector<T> v(count);
    for (Int32 i = count - 1; i >= 0; i--) {
      v[i] = x - offset;
      x *= edelta;
    }
    return v;
  }

  std::vector<T> SpreadOverLog(T delta, T offset = 0.,
                               bool backward = false) const {
    static_assert(std::is_same<T, Float64>::value,
                  "not implemented"); // compile time check
    if (!isSameSign(offset))
      THROWG(INTERNAL_ERROR, "borders should be of same sign");

    if (backward)
      return SpreadOverLog_backward(delta, offset);

    if (GetIsEmpty() || delta == 0.0 ||
        GetLength() <
            (GetBegin() + offset) * exp(delta) - (GetBegin() + offset))
      return {m_Begin};

    T x = m_Begin + offset;
    T edelta = exp(delta);

    Int32 count =
        (log(GetEnd() + offset) - log(GetBegin() + offset)) / delta + 1;
    std::vector<T> v(count);
    for (Int32 i = 0; i < count; i++) {
      v[i] = x - offset;
      x *= edelta;
    }
    return v;
  }

  std::vector<T> SpreadOverLogEpsilon(T delta, T offset = 0.,
                                      bool backward = false,
                                      Float64 epsilon = 1e-8) const {
    static_assert(std::is_same<T, Float64>::value,
                  "not implemented"); // compile time check
    CRange<T> range_epsilon =
        backward ? CRange<T>{-epsilon, 0.0} : CRange<T>{0.0, epsilon};
    CRange<T> range_ = *this + range_epsilon;
    auto ret = range_.SpreadOverLog(delta, offset, backward);
    ret.front() = std::max(ret.front(), m_Begin);
    ret.back() = std::min(ret.back(), m_End);
    return ret;
  }

  // spread over log (z+1)
  std::vector<T> SpreadOverLogZplusOne(T delta, bool backward = false) const {
    static_assert(std::is_same<T, Float64>::value,
                  "not implemented"); // compile time check
    return SpreadOverLog(delta, 1., backward);
  }

  std::vector<T> SpreadOverLogZplusOneEpsilon(T delta, bool backward = false,
                                              Float64 epsilon = 1e-8) const {
    static_assert(std::is_same<T, Float64>::value,
                  "not implemented"); // compile time check
    CRange<T> range_epsilon =
        backward ? CRange<T>{-epsilon, 0.0} : CRange<T>{0.0, epsilon};
    CRange<T> range_ = *this + range_epsilon;
    auto ret = range_.SpreadOverLogZplusOne(delta, backward);
    ret.front() = std::max(ret.front(), m_Begin);
    ret.back() = std::min(ret.back(), m_End);
    return ret;
  }

  //  template<typename T>
  friend std::ostream &operator<<(std::ostream &out, const CRange<T> &range) {
    out << "[" << range.m_Begin << "," << range.m_End << "]";
    return out;
  }

  friend std::istream &operator>>(std::istream &in, CRange<T> &range) {
    in >> range.m_Begin;
    range.m_End = range.m_Begin;

    return in;
  }

  // create a centered vector
  const std::vector<T> spanCenteredWindow(T center, bool logsampling,
                                          T delta) const {
    CRange<T> range_left(GetBegin(), center);
    std::vector<T> vect =
        logsampling ? range_left.SpreadOverLogZplusOneEpsilon(delta, true)
                    : range_left.SpreadOverEpsilon(delta, true);

    CRange<T> range_right = {center, GetEnd()};
    std::vector<T> vect_part2 =
        logsampling ? range_right.SpreadOverLogZplusOneEpsilon(delta)
                    : range_right.SpreadOverEpsilon(delta);

    vect.insert(std::end(vect), std::begin(vect_part2) + 1,
                std::end(vect_part2));
    return vect;
  }

  // enclosed refers to having i_max referring to m_End or higher and i_min
  // referring to m_Begin or lower
  bool getEnclosingIntervalIndices(const std::vector<T> &ordered_values,
                                   const T &value, Int32 &i_min, Int32 &i_max,
                                   bool warning = true) const {
    if (value < m_Begin || value > m_End) {
      if (warning)
        Flag.warning(WarningCode::CRANGE_VALUE_OUTSIDERANGE,
                     Formatter()
                         << "CRange::" << __func__ << ": value " << value
                         << " not inside ]" << m_Begin << "," << m_End << "[");
      return false;
    } else if (m_Begin < ordered_values.front() ||
               m_End > ordered_values.back()) {
      if (warning)
        Flag.warning(WarningCode::CRANGE_VECTBORDERS_OUTSIDERANGE,
                     Formatter()
                         << "CRange::" << __func__ << ": ]" << m_Begin << ","
                         << m_End << "[ not inside ordered_values");
      return false;
    }
    typename std::vector<T>::const_iterator it =
        std::lower_bound(ordered_values.begin(), ordered_values.end(), value);
    typename std::vector<T>::const_iterator it_min =
        std::lower_bound(ordered_values.begin(), it, m_Begin);
    typename std::vector<T>::const_iterator it_max =
        std::lower_bound(it, ordered_values.end(), m_End);

    if (*it_min > m_Begin)
      it_min = it_min - 1;

    i_min = it_min - ordered_values.begin();
    i_max = it_max - ordered_values.begin();
    return true;
  }

  bool getEnclosingIntervalIndices(const std::vector<T> &ordered_values,
                                   Int32 &i_min, Int32 &i_max,
                                   bool warning = true) const {
    if (m_Begin < ordered_values.front() || m_End > ordered_values.back()) {
      if (warning)
        Flag.warning(WarningCode::CRANGE_VECTBORDERS_OUTSIDERANGE,
                     Formatter()
                         << "CRange::" << __func__ << ": ]" << m_Begin << ","
                         << m_End << "[ not inside ordered_values");
      return false;
    }

    typename std::vector<T>::const_iterator it_min =
        std::lower_bound(ordered_values.begin(), ordered_values.end(), m_Begin);
    typename std::vector<T>::const_iterator it_max =
        std::lower_bound(ordered_values.begin(), ordered_values.end(), m_End);

    if (*it_min > m_Begin)
      it_min = it_min - 1;

    i_min = it_min - ordered_values.begin();
    i_max = it_max - ordered_values.begin();
    return true;
  }

  // closed refers to having i_min referring to m_Begin index or higher and
  // i_max referring to m_End index or lower
  bool getClosedIntervalIndices(const std::vector<T> &ordered_values,
                                Int32 &i_min, Int32 &i_max,
                                bool warning = true) const {
    if (m_End < ordered_values.front() || m_Begin > ordered_values.back()) {
      if (warning)
        Flag.warning(WarningCode::CRANGE_VECTBORDERS_OUTSIDERANGE,
                     Formatter()
                         << "CRange::" << __func__ << ": ]" << m_Begin << ","
                         << m_End << "[ not inside ordered_values");
      return false;
    }

    typename std::vector<T>::const_iterator it_min =
        std::lower_bound(ordered_values.begin(), ordered_values.end(), m_Begin);
    typename std::vector<T>::const_iterator it_max =
        std::lower_bound(ordered_values.begin(), ordered_values.end(), m_End);

    if (it_max == ordered_values.end())
      --it_max;
    else if (*it_max > m_End)
      --it_max;

    i_min = it_min - ordered_values.begin();
    i_max = it_max - ordered_values.begin();
    if (i_min > i_max) {
      if (warning)
        Flag.warning(
            WarningCode::CRANGE_NO_INTERSECTION,
            Formatter()
                << "CRange::" << __func__
                << ": There is no sample inside range (min,max indices=["
                << i_min << "," << i_max << "]");
      return false;
    }
    return true;
  }

  static bool HasIntersection(const CRange<T> &a, const CRange<T> &b) {
    return std::max(a.GetBegin(), b.GetBegin()) <=
           std::min(a.GetEnd(), b.GetEnd());
  }
  bool HasIntersectionWith(const CRange<T> &r) {
    return HasIntersection(*this, r);
  }

  bool IntersectWith(const CRange<T> &r) { return Intersect(*this, r, *this); }

  static bool Intersect(const CRange<T> &a, const CRange<T> &b,
                        CRange<T> &intersect) {
    if (HasIntersection(a, b)) {
      intersect.Set(std::max(a.GetBegin(), b.GetBegin()),
                    std::min(a.GetEnd(), b.GetEnd()));
      return true;
    }

    return false;
  }

  static bool getUnion(const CRange<T> &a, const CRange<T> &b, CRange<T> &ab) {
    if (HasIntersection(a, b)) {
      ab.Set(std::min(a.GetBegin(), b.GetBegin()),
             std::max(a.GetEnd(), b.GetEnd()));
      return true;
    }
    return false;
  }
  bool unionWith(const CRange<T> &r) { return getUnion(*this, r, *this); }

  static std::vector<CRange<T>>
  joinIntersections(std::vector<CRange<T>> ranges) {

    if (ranges.size() == 1)
      return ranges;

    std::vector<CRange<T>> result;
    std::sort(ranges.begin(), ranges.end());
    auto current = ranges.front();
    for (auto it = ranges.cbegin() + 1, end = ranges.cend(); it != end; ++it) {
      if (current.GetEnd() >= it->GetBegin()) {
        current.SetEnd(std::max(current.GetEnd(), it->GetEnd()));
      } else {
        result.push_back(current);
        current = *it;
      }
    }
    result.push_back(current);
    return result;
  }

private:
  T m_Begin;
  T m_End;
};

typedef CRange<Int32> TInt32Range;
typedef CRange<Float64> TFloat64Range;
typedef TFloat64Range TLambdaRange;
typedef std::vector<TLambdaRange> TLambdaRangeList;
typedef std::vector<TInt32Range> TInt32RangeList;
typedef std::vector<TFloat64Range> TFloat64RangeList;
using CTLambdaRangePtrVector = std::vector<std::shared_ptr<const TLambdaRange>>;
} // namespace NSEpic

#endif
