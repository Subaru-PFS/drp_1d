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

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/log/log.h>

#include <cmath>
#include <vector>

namespace NSEpic {

/**
 * \ingroup Redshift
 * Templated Range manipulation class
 */
template <typename T> class CRange
{

  public:
    CRange(): m_Begin(), m_End() {}

    CRange(const T begin, const T end):
        m_Begin(begin), m_End(end) {}

    CRange(const std::vector<T> & v):
        m_Begin(v.empty() ? T() : v.front()), 
        m_End(v.empty() ?  T() : v.back()) 
        {}

    ~CRange() {}

    bool GetIsEmpty() const { return m_Begin == m_End; }

    void Set(const T &b, const T &e)
    {
        m_Begin = b;
        m_End = e;
    }

    void SetBegin(const T &v) { m_Begin = v; }

    void SetEnd(const T &v) { m_End = v; }

    CRange<T> operator-(T t) { return CRange<T>(m_Begin - t, m_End - t); }
    CRange<T> operator+(T t) { return CRange<T>(m_Begin + t, m_End + t); }
    CRange<T> operator*(T t) { return CRange<T>(m_Begin * t, m_End * t); }
    CRange<T> operator/(T t) { return CRange<T>(m_Begin / t, m_End / t); }

    T Interpolate(Float64 t) { return m_Begin + (m_End - m_Begin) * t; }

    const T &GetBegin() const { return m_Begin; }

    const T &GetEnd() const { return m_End; }

    T GetLength() const { return m_End - m_Begin; }

    bool IntersectWith(const CRange<T> r)
    {
        return Intersect(*this, r, *this);
    }

    static bool Intersect(const CRange<T> &a, const CRange<T> b,
                          CRange<T> &intersect)
    {
        if ((a.GetBegin() >= b.GetBegin() && a.GetBegin() <= b.GetEnd()) ||
            (a.GetEnd() >= b.GetBegin() && a.GetEnd() <= b.GetEnd()) ||
            (b.GetBegin() >= a.GetBegin() && b.GetBegin() <= a.GetEnd()) ||
            (b.GetEnd() >= a.GetBegin() && b.GetEnd() <= a.GetEnd()))
        {
            intersect.Set(std::max(a.GetBegin(), b.GetBegin()),
                          std::min(a.GetEnd(), b.GetEnd()));
            return true;
        }

        return false;
    }

    std::vector<T> SpreadOver(T delta) const
    {
        std::vector<T> v;

        if (GetIsEmpty() || delta == 0.0 || GetLength() < delta)
        {
            v.resize(1);
            v[0] = m_Begin;
            return v;
        }

        Int32 count = (GetLength() + epsilon) / delta;

        v.resize(count + 1);
        for (UInt32 i = 0; i < v.size(); i++)
        {
            v[i] = GetBegin() + delta * i;
        }

        return v;
    }

    std::vector<T> SpreadOverLog(T delta, T offset=0.) const
    {
        static_assert(std::is_same<T, Float64>::value,"not implemented");// compile time check 

        std::vector<T> v;
        if (GetIsEmpty() || delta == 0.0 || GetLength() < (GetBegin() + offset)*exp(delta) - (GetBegin() + offset) + epsilon)
        {
            v.resize(1);
            v[0] = m_Begin;
            return v;
        }

        T x = m_Begin + offset;
        T edelta = exp(delta);
        Int32 count = 0;
        Int32 maxCount = 1e8;
        while (x < (m_End + offset + epsilon) && count < maxCount)
        {
            v.push_back(x - offset);
            count++;
            x *= edelta;
        }
        return v;
    }  
    //spread over log (z+1)
    std::vector<T> SpreadOverLogZplusOne(T delta) const
    {
        static_assert(std::is_same<T, Float64>::value,"not implemented");// compile time check 
        return SpreadOverLog(delta, 1.);
    }
    //  template<typename T>
  friend std::ostream& operator<< (std::ostream &out, const CRange<T> &range)
  {
    out <<"["<< range.m_Begin<<","<< range.m_End<<"]";
    return out;
  }

  friend std::istream& operator>> (std::istream &in, CRange<T> &range)
  {
    in >> range.m_Begin;
    range.m_End = range.m_Begin;

    return in;
 }

  //enclosed refers to having i_max referring to m_End or higher and i_min referring to m_Begin or lower
  bool getEnclosingIntervalIndices(const std::vector<T>& ordered_values,const T& value, Int32& i_min,Int32& i_max, bool warning=true) const
  {
    if (value < m_Begin || value > m_End)
      {
        if (warning) Log.LogWarning("CRange::getEnclosingIntervalIndices: value %.5f not inside ]%.5f,%.5f[",value,m_Begin,m_End);
        return false;
      }
    else if(m_Begin < ordered_values.front() || m_End > ordered_values.back())
      {
        if (warning) Log.LogWarning("CRange::getEnclosingIntervalIndices: ]%.5f,%.5f[ not inside ordered_values",m_Begin,m_End);
        return false;
      }
    typename std::vector<T>::const_iterator it = std::lower_bound(ordered_values.begin(),ordered_values.end(),value);
    typename std::vector<T>::const_iterator it_min = std::lower_bound(ordered_values.begin(),it,m_Begin);
    typename std::vector<T>::const_iterator it_max = std::lower_bound(it,ordered_values.end(),m_End);

    if(*it_min > m_Begin) it_min = it_min -1;

    i_min = it_min - ordered_values.begin();
    i_max = it_max - ordered_values.begin();
    return true;
  }

  bool getEnclosingIntervalIndices(const std::vector<T>& ordered_values,Int32& i_min,Int32& i_max, bool warning=true) const
  {
    if(m_Begin < ordered_values.front() || m_End > ordered_values.back())
      {
        if (warning) Log.LogWarning("CRange::getEnclosingIntervalIndices: ]%.5f,%.5f[ not inside ordered_values",m_Begin,m_End);
        return false;
      }
    
    typename std::vector<T>::const_iterator it_min = std::lower_bound(ordered_values.begin(),ordered_values.end(),m_Begin);
    typename std::vector<T>::const_iterator it_max = std::lower_bound(ordered_values.begin(),ordered_values.end(),m_End);

    if(*it_min > m_Begin) it_min = it_min -1;

    i_min = it_min - ordered_values.begin();
    i_max = it_max - ordered_values.begin();
    return true;
  }

  //closed refers to having i_min referring to m_Begin index or higher and i_max referring to m_End index or lower  
  bool getClosedIntervalIndices(const std::vector<T>& ordered_values,Int32& i_min,Int32& i_max, bool warning=true) const
  {
    if(m_End < ordered_values.front() || m_Begin > ordered_values.back())
      {
        if (warning) Log.LogWarning("CRange::getClosedIntervalIndices: ]%.5f,%.5f[ not inside ordered_values",m_Begin,m_End);
        return false;
      }
    
    typename std::vector<T>::const_iterator it_min = std::lower_bound(ordered_values.begin(),ordered_values.end(),m_Begin);
    typename std::vector<T>::const_iterator it_max = std::lower_bound(ordered_values.begin(),ordered_values.end(),m_End);
    
    if (it_max==ordered_values.end()) --it_max;
    else if(*it_max > m_End) --it_max;

    i_min = it_min - ordered_values.begin();
    i_max = it_max - ordered_values.begin();
    if(i_min>i_max){
        if (warning) Log.LogWarning("CRange::getClosedIntervalIndices: There is no sample inside range (min,max indices=[%d,%d]",i_min, i_max);
        return false;
    }
    return true;
  } 

  private:
    T m_Begin;
    T m_End;
    Float64 epsilon = 1E-6;
};

typedef CRange<Int32> TInt32Range;
typedef CRange<Float64> TFloat64Range;
typedef TFloat64Range TLambdaRange;
typedef std::vector<TLambdaRange> TLambdaRangeList;
typedef std::vector<TInt32Range> TInt32RangeList;
typedef std::vector<TFloat64Range> TFloat64RangeList;

} // namespace NSEpic

#endif
