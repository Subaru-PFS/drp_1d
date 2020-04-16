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
    CRange() {}

    CRange(const T begin, const T end):
        m_Begin(begin), m_End(end) {}

    CRange(const std::vector<T> & v):
        m_Begin(v.front()), m_End(v.back()) {}

    ~CRange() {}

    Bool GetIsEmpty() const { return m_Begin == m_End; }

    void Set(const T &b, const T &e)
    {
        m_Begin = b;
        m_End = e;
    }

    void SetBegin(const T &v) { m_Begin = v; }

    void SetEnd(const T &v) { m_End = v; }

    CRange<T> operator-(Float64 t) { return CRange<T>(m_Begin - t, m_End - t); }

    T Interpolate(Float64 t) { return m_Begin + (m_End - m_Begin) * t; }

    const T &GetBegin() const { return m_Begin; }

    const T &GetEnd() const { return m_End; }

    T GetLength() const { return m_End - m_Begin; }

    Bool IntersectWith(const CRange<T> r)
    {
        return Intersect(*this, r, *this);
    }

    static Bool Intersect(const CRange<T> &a, const CRange<T> b,
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

    std::vector<T> SpreadOver(Float64 delta) const
    {
        std::vector<T> v;

        if (GetIsEmpty() || delta == 0.0 || GetLength() < delta)
        {
            v.resize(1);
            v[0] = m_Begin;
            return v;
        }

        Int32 count = GetLength() / delta;

        v.resize(count + 1);
        for (UInt32 i = 0; i < v.size(); i++)
        {
            v[i] = GetBegin() + delta * i;
        }

        return v;
    }

    std::vector<T> SpreadOverLog(Float64 delta) const
    {
	std::vector<T> v;
	if (GetIsEmpty() || delta == 0.0 || GetLength() < delta)
	{
	    v.resize(1);
	    v[0] = m_Begin;
	    return v;
	}

        Float64 x = m_Begin + 1.;
        Float64 edelta = exp(delta);
        Int32 count = 0;
        Int32 maxCount = 1e8;
        while (x < m_End+1. && count < maxCount)
        {
            v.push_back(x-1.);
            count++;
            x *= edelta;
        }

        return v;
    }
  //enclosed refers to having i_max referring to m_End or higher and i_min referring to m_Begin or lower
  bool getEnclosingIntervalIndices(std::vector<T>& ordered_values,const T& value, Int32& i_min,Int32& i_max) const
  {
    if (value < m_Begin || value > m_End)
      {
        Log.LogError("%.5f not inside ]%.5f,%.5f[",value,m_Begin,m_End);
        return false;
      }
    else if(m_Begin < ordered_values.front() || m_End > ordered_values.back())
      {
        Log.LogError("]%.5f,%.5f[ not inside ordered_values",m_Begin,m_End);
        return false;
      }
    typename std::vector<T>::iterator it = std::lower_bound(ordered_values.begin(),ordered_values.end(),value);
    typename std::vector<T>::iterator it_min = std::lower_bound(ordered_values.begin(),it,m_Begin);
    typename std::vector<T>::iterator it_max = std::lower_bound(it,ordered_values.end(),m_End);

    if(*it_min > m_Begin) it_min = it_min -1;

    i_min = it_min - ordered_values.begin();
    i_max = it_max - ordered_values.begin();
    return true;
  }

  bool getEnclosingIntervalIndices(std::vector<T>& ordered_values,Int32& i_min,Int32& i_max) const
  {
    if(m_Begin < ordered_values.front() || m_End > ordered_values.back())
      {
        Log.LogError("]%.5f,%.5f[ not inside ordered_values",m_Begin,m_End);
        return false;
      }
    
    typename std::vector<T>::iterator it_min = std::lower_bound(ordered_values.begin(),ordered_values.end(),m_Begin);
    typename std::vector<T>::iterator it_max = std::lower_bound(ordered_values.begin(),ordered_values.end(),m_End);

    if(*it_min > m_Begin) it_min = it_min -1;

    i_min = it_min - ordered_values.begin();
    i_max = it_max - ordered_values.begin();
    return true;
  }

    bool getExactEnclosingIntervalIndices(std::vector<T>& ordered_values,Int32& i_min,Int32& i_max) const
  {
    if(m_Begin < ordered_values.front() || m_End > ordered_values.back())
      {
        Log.LogError("]%.5f,%.5f[ not inside ordered_values",m_Begin,m_End);
        return false;
      }
    
    typename std::vector<T>::iterator it_min = std::lower_bound(ordered_values.begin(),ordered_values.end(),m_Begin);
    typename std::vector<T>::iterator it_max = std::lower_bound(ordered_values.begin(),ordered_values.end(),m_End);

    if(*it_min != m_Begin)
      {
        Log.LogError("Cannot find %.5f in ordered_values",m_Begin);
        return false;
      }
    if(*it_max != m_End)
     {
        Log.LogError("Cannot find %.5f in ordered_values",m_End);
        return false;
      }

    i_min = it_min - ordered_values.begin();
    i_max = it_max - ordered_values.begin();
    return true;
  }
  //closed refers to having i_min referring to m_Begin index or higher and i_max referring to m_End index or lower  
  bool getClosedIntervalIndices(std::vector<T>& ordered_values,Int32& i_min,Int32& i_max) const
  {
    if(m_Begin < ordered_values.front() || m_End > ordered_values.back())
      {
        Log.LogError("]%.5f,%.5f[ not inside ordered_values",m_Begin,m_End);
        return false;
      }
    
    typename std::vector<T>::iterator it_min = std::lower_bound(ordered_values.begin(),ordered_values.end(),m_Begin);
    typename std::vector<T>::iterator it_max = std::lower_bound(ordered_values.begin(),ordered_values.end(),m_End);

    if(*it_max > m_End) it_max = it_max -1;

    i_min = it_min - ordered_values.begin();
    i_max = it_max - ordered_values.begin();
    return true;
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

} // namespace NSEpic

#endif