#ifndef _CORE_COMMON_RANGE_
#define _CORE_COMMON_RANGE_

#include <RedshiftLibrary/common/datatypes.h>

#include <vector>

namespace NSEpic
{

/**
 * \ingroup Core
 * Templated Range manipulation class
 */
template <typename T> 
class CRange
{

public:

    CRange( )
    {
    }

    CRange( T begin, T end )
    {
        m_Begin = begin;
        m_End = end;
    }

    ~CRange( )
    {  
    }

    Bool GetIsEmpty() const
    {
        return m_Begin == m_End;
    }

    void Set( const T& b,  const T& e )
    {
        m_Begin = b;
        m_End = e;
    }

    void SetBegin( const T& v )
    {
        m_Begin = v;
    }

    void SetEnd( const T& v )
    {
        m_End = v;
    }

    CRange<T> operator - ( Float64 t )
    {
        return CRange<T>( m_Begin - t, m_End - t );
    }

    T Interpolate( Float64 t )
    {
        return m_Begin + ( m_End - m_Begin ) * t;
    }

    const T& GetBegin() const
    {
        return m_Begin;
    }

    const T& GetEnd() const
    {
        return m_End;
    }

    T GetLength() const
    {
        return m_End - m_Begin;
    }

    static Bool Intersect( const CRange<T>& a, const CRange<T> b, CRange<T>& intersect )
    {
        if( ( a.GetBegin() >= b.GetBegin() && a.GetBegin() <= b.GetEnd() ) ||
                ( a.GetEnd() >= b.GetBegin() && a.GetEnd() <= b.GetEnd() ) ||
                ( b.GetBegin() >= a.GetBegin() && b.GetBegin() <= a.GetEnd() ) ||
                            ( b.GetEnd() >= a.GetBegin() && b.GetEnd() <= a.GetEnd() ))
        {
            intersect.Set( std::max( a.GetBegin(), b.GetBegin() ), std::min( a.GetEnd(), b.GetEnd() ) );
            return true;
        }

        return false;
    }

    std::vector<T>   SpreadOver( Float64 delta ) const
    {
        std::vector<T> v;

        if( GetIsEmpty() || delta == 0.0  || GetLength() < delta )
        {
            v.resize( 1 );
            v[0] = m_Begin;
            return v;
        }

        Int32 count = GetLength() / delta;

        v.resize( count+1 );
        for( UInt32 i=0; i<v.size() ;i++ )
        {
            v[i] = GetBegin() + delta * i;
        }

        return v;
    }

    std::vector<T>   SpreadOverOnePlusX( Float64 delta ) const
    {
        std::vector<T> v;

        if( GetIsEmpty() || delta == 0.0  || GetLength() < delta )
        {
            v.resize( 1 );
            v[0] = m_Begin;
            return v;
        }


        v.push_back(m_Begin);
        Float64 x = m_Begin;
        Float64 step = delta/(1.+m_Begin);
        Int32 count = 0;
        Int32 maxCount = 1e12;
        while(x+step<m_End && count<maxCount)
        {
            x = x+step;
            v.push_back(x);
            step = delta/(1.+x);
            count ++;
        }

        return v;
    }

private:

    T     m_Begin;
    T     m_End;
};

typedef CRange<Int32>   TInt32Range;
typedef CRange<Float64> TFloat64Range;
typedef TFloat64Range   TLambdaRange;
typedef std::vector< TLambdaRange > TLambdaRangeList;
typedef std::vector< TInt32Range > TInt32RangeList;
typedef std::vector< TFloat64Range > TFloat64RangeList;

}

#endif
