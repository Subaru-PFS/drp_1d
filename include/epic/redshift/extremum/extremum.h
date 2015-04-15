#ifndef _REDSHIFT_EXTREMUM_EXTREMUM_
#define _REDSHIFT_EXTREMUM_EXTREMUM_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>

#include <vector>

namespace NSEpic
{

/**
 * Analyse a given input 2D array and find the n stronger extremum.
 *
 */
class CExtremum
{


public:

    CExtremum();
    CExtremum( const TFloat64Range& xRange, UInt32 maxPeakCount = 5, UInt32 refreshCount = 2 );
    ~CExtremum();

    void SetMaxPeakCount( UInt32 n );
    void SetRefreshCount( UInt32 n );
    void SetXRange( const TFloat64Range& r );

    Bool Find( const Float64* xAxis, const Float64* yAxis, UInt32 n, TPointList& maxPoint ) const;

private:

    Bool InternalFind2( const Float64* xAxis, const Float64* yAxis, UInt32 n, TPointList& maxPoint ) const;
    Bool InternalFind( const Float64* xAxis, const Float64* yAxis, UInt32 n, TPointList& maxPoint ) const;

    UInt32          m_MaxPeakCount;
    UInt32          m_RefreshCount;
    TFloat64Range   m_XRange;

};


}

#endif
