#ifndef _REDSHIFT_EXTREMUM_EXTREMUM_
#define _REDSHIFT_EXTREMUM_EXTREMUM_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>

#include <vector>

namespace NSEpic
{

/**
 * \ingroup Redshift
 * Analyse a given input 2D array and find the n stronger extremum.
 */
class CExtremum
{


public:

    CExtremum( Bool invertForMinSearch = false );
    CExtremum( const TFloat64Range& xRange, UInt32 maxPeakCount = 5, Float64 radius = 0.005, Bool invertForMinSearch=false, UInt32 refreshCount = 2);
    ~CExtremum();

    void SetMaxPeakCount( UInt32 n );
    void SetRefreshCount( UInt32 n );
    void SetXRange( const TFloat64Range& r );
    void SetSignSearch( Float64 val );
    Bool Find( const TFloat64List& xAxis, const TFloat64List& yAxis, TPointList& maxPoint ) const;
    Bool DefaultExtremum( const TFloat64List& xAxis, const TFloat64List& yAxis, TPointList& maxPoint ) ;
    Bool Cut_Threshold( TPointList& maxPoint, UInt32 meritCut, Float64 bestPDF, Int32 keepMinN) const;
private:

    Bool InternalFind2( const Float64* xAxis, const Float64* yAxis, UInt32 n, TPointList& maxPoint ) const;
    Bool InternalFind( const Float64* xAxis, const Float64* yAxis, UInt32 n, TPointList& maxPoint ) const;

    //Mira:
    Bool FindPeaks_extended(const TFloat64List& xAxis, const TFloat64List& yAxis, UInt32 n, TPointList& maxPoint ) const; 
    Bool InternalFind_refact_ext(const Float64* xAxis, const Float64* yAxis, UInt32 n, TPointList& maxPoint) const;
    
    UInt32          m_MaxPeakCount;
    UInt32          m_RefreshCount;
    TFloat64Range   m_XRange;

    Float64         m_SignSearch;
    Float64         m_Radius;
};


}

#endif
