#ifndef _REDSHIFT_EXTREMUM_EXTREMUM_
#define _REDSHIFT_EXTREMUM_EXTREMUM_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>

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

    CExtremum( Bool invertForMinSearch = false );
    CExtremum( const TFloat64Range& xRange, UInt32 maxPeakCount = 10, Float64 radius = 0.005, Bool invertForMinSearch=false);
    ~CExtremum();

    void SetMaxPeakCount( UInt32 n );
    void SetXRange( const TFloat64Range& r );
    void SetSignSearch( Float64 val );
    void SetMeritCut( UInt32 n );
    Bool Find( const TFloat64List& xAxis, const TFloat64List& yAxis, TPointList& maxPoint ) const;
    Bool DefaultExtremum( const TFloat64List& xAxis, const TFloat64List& yAxis, TPointList& maxPoint );
    //made it public to do unit tests
    Bool Cut_Threshold( vector <Float64>& maxX, vector <Float64>& maxY, Int32 keepMinN) const;
private:
    Bool FindAllPeaks(const Float64* xAxis, const Float64* yAxis, UInt32 n, vector <Float64>& maxX, vector <Float64>& maxY) const;
    Bool FilterOutNeighboringPeaks( vector <Float64>& maxX, vector <Float64>& maxY, UInt32 n) const;
    
    Bool Truncate( vector <Float64>& xAxis, vector <Float64>& yAxis, Int32 maxCount, TPointList& maxPoint) const;

    UInt32          m_MaxPeakCount;
    TFloat64Range   m_XRange;
    Float64         m_meritCut;
    Float64         m_SignSearch;
    Float64         m_Radius;

    //TO change
    UInt32 m_maxCount;
};


}

#endif
