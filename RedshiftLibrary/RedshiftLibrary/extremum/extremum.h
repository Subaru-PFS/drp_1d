#ifndef _REDSHIFT_EXTREMUM_EXTREMUM_
#define _REDSHIFT_EXTREMUM_EXTREMUM_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <iostream>
#include <vector>
using namespace std; 
namespace NSEpic
{

/**
 * \ingroup Redshift
 * Analyse a given input 2D array and find the n stronger extremum.
 */
class CExtremum
{


public:
    
    CExtremum( UInt32 maxPeakCount=10,
                Float64 peakSeparation= 0.005*2, 
                Float64 meritcut = -1, // <0 = no cut thresholding 
                Bool invertForMinSearch=false, 
                const TFloat64Range& xRange=TFloat64Range());

    void SetMaxPeakCount( UInt32 n );
    void SetXRange( const TFloat64Range& r );
    void SetSignSearch( Float64 val );
    void SetMeritCut( Float64 n );
    Bool Find( const TFloat64List& xAxis, const TFloat64List& yAxis, TPointList& maxPoint ) const;
    Bool DefaultExtremum( const TFloat64List& xAxis, const TFloat64List& yAxis, TPointList& maxPoint );

    Bool Cut_Threshold( TFloat64List& maxX, TFloat64List& maxY, Int32 keepMinN) const;
    //made public to do unit tests
    Bool Truncate( TFloat64List& xAxis, TFloat64List& yAxis, TPointList& maxPoint) const;
    Bool FilterOutNeighboringPeaksAndTruncate(TFloat64List& maxX, TFloat64List& maxY, UInt32 keepmin, TPointList& maxPoint)const;
    void SortIndexes(TFloat64List&  maxY) const;

private:
    Bool FindAllPeaks(const Float64* xAxis, const Float64* yAxis, UInt32 n, TFloat64List& maxX, TFloat64List& maxY) const;
    Bool FindAllPeaks(const Float64* xAxis, const Float64* yAxis, UInt32 n, TFloat64List& maxX, TFloat64List& maxY, Float64 SignSearch) const;
    TFloat64List Cut_Prominence_Merit( TFloat64List& maxX, TFloat64List& maxY, TFloat64List& minX, TFloat64List& minY) const;
 
    Bool verifyPeakSeparation( TFloat64List& maxX) const;
    Bool verifyPeakSeparation( TPointList& maxPoint) const;
    UInt32          m_MaxPeakCount;
    TFloat64Range   m_XRange;
    Float64         m_meritCut;
    Float64         m_SignSearch;
    Float64         m_extrema_separation;
    Bool            m_PeakSeparationActive = true;
    mutable TInt32List   m_sortedIndexes;
};

}

#endif
