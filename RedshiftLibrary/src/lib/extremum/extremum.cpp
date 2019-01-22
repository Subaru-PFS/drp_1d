#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/common/quicksort.h>
#include <RedshiftLibrary/common/datatypes.h>

#include <cmath>
#include <float.h>
#include <iostream>

#include <stdio.h>

#define PEAKS_MIN_THRESHOLD (3)
#define PEAKS_SMOOTH_LIMIT (20)

using namespace NSEpic;
using namespace std;

/**
 * Constructs CExtremum with default values, and SetSignSearch ( -1.0 ) if argument is true, SetSignSearch ( 1.0 ) otherwise.
 */
CExtremum::CExtremum( Bool invertForMinSearch ) :
    m_MaxPeakCount( 5 ),
    m_RefreshCount( 1 ),
    m_XRange( 0.0, 0.0 ),
    m_SignSearch( 1.0 )
{
    if( invertForMinSearch )
      {
        SetSignSearch( -1.0 );
      }
    else
      {
        SetSignSearch( 1.0 );
      }
}

/**
 * Member attribution constructor.
 */
CExtremum::CExtremum( const TFloat64Range& xRange, UInt32 maxPeakCount, Bool invertForMinSearch, UInt32 refreshCount ) :
    m_MaxPeakCount( maxPeakCount ),
    m_RefreshCount( refreshCount ),
    m_XRange( xRange ),
    m_SignSearch( 1.0 )
{
    if( invertForMinSearch )
      {
        SetSignSearch( -1.0 );
      }
    else
      {
        SetSignSearch( 1.0 );
      }
}

/**
 * Empty destructor.
 */
CExtremum::~CExtremum()
{

}

/**
 *  Sets m_XRange to r.
 */
void CExtremum::SetXRange( const TFloat64Range& r )
{
    m_XRange = r;
}

/**
 * Sets m_MaxPeakCount to n.
 */
void CExtremum::SetMaxPeakCount( UInt32 n )
{
    m_MaxPeakCount = n;
}

/**
 * Sets m_RefreshCount to n.
 */
void CExtremum::SetRefreshCount( UInt32 n )
{
    m_RefreshCount = n;
}

/**
 * Sets m_SignSearch to val.
 */
void CExtremum::SetSignSearch( Float64 val )
{
    m_SignSearch = val;
}

/**
 * Wrapper around InternalFind, this function validates and sets the search interval limits.
 */
Bool CExtremum::Find( const TFloat64List& xAxis, const TFloat64List& yAxis, TPointList& maxPoint ) const
{
    Int32 n = xAxis.size();
    const Float64* selectedXAxis = xAxis.data();
    const Float64* selectedYAxis = yAxis.data();

    Int32 rangeXBeginIndex = -1;
    Int32 rangeXEndIndex = -1;

    if( n == 0 )
        return false;

    // Find index in xAxis that correspond to the boundary specified by m_XRange
    if( m_XRange.GetIsEmpty() )
    {
        rangeXBeginIndex = 0;
        rangeXEndIndex = n-1;
    }
    else
    {
    	// Find index range for the given lambda range
        for( UInt32 i=0; i<n; i++ )
        {
            if( rangeXBeginIndex==-1 && xAxis[i]>=m_XRange.GetBegin() )
            {
                rangeXBeginIndex = i;
            }

            if( xAxis[i]<=m_XRange.GetEnd() )
            {
                rangeXEndIndex = i;
            }
        }
    }

    return InternalFind( selectedXAxis+rangeXBeginIndex, selectedYAxis+rangeXBeginIndex, (rangeXEndIndex-rangeXBeginIndex)+1, maxPoint );
}

/**
 * Attempts to find peaks, returning (when appropriate) the points where they reside in the maxPoint argument.
 */
Bool CExtremum::InternalFind( const Float64* xAxis, const Float64* yAxis, UInt32 n, TPointList& maxPoint ) const
{
    if( n == 0 )
        return false;

    //Method 1, use only 1 extremum
    /*
    maxPoint.resize( 1 );

    Float64 max = DBL_MIN ;
    Int32 maxIndex = 0;
    for( Int32 i=0; i<n; i++ )
    {
    	if( yAxis[i] > max ) {
    		max = yAxis[i];
    		maxIndex = i;
    	}
    }

    maxPoint[0].X = xAxis[maxIndex];
    maxPoint[0].Y = yAxis[maxIndex];

    return true;
    //*/

    vector<Float64> maxX( n );
    vector<Float64> maxY( n );

    // Tmp array can be considered as the "input" of each iteration.
    vector<Float64> tmpX( n );
    vector<Float64> tmpY( n );
    UInt32 tmpSize = n;
    for( UInt32 t=0; t<n; t++ )
    {
        tmpX[t] = xAxis[t];
        tmpY[t] = m_SignSearch*yAxis[t];
    }


    for( UInt32 count=0; count<m_RefreshCount; count++ )
    {
        Int32 maxCount = 0;

        // find first and last non Nan element
        Int32 firstNonNanInd = 0;
        for( Int32 iFirst=0; iFirst<tmpSize-1; iFirst++ )
	  {
            if( !std::isnan( (double) tmpY[iFirst] ) )
	      {
                firstNonNanInd = iFirst;
                break;
	      }
	  }
        Int32 lastNonNanInd = tmpSize-1;
        for( Int32 iLast=tmpSize-1; iLast>0; iLast-- )
	  {
            if( !std::isnan( tmpY[iLast] ) )
	      {
                lastNonNanInd = iLast;
                break;
	      }
	  }

        // First element
        if( tmpY[firstNonNanInd] > tmpY[firstNonNanInd+1] )
	  {
            maxX[maxCount] = tmpX[firstNonNanInd];
            maxY[maxCount] = tmpY[firstNonNanInd];
            maxCount++;
	  }

        for( Int32 i=firstNonNanInd+1; i<lastNonNanInd; i++ )
	  {
            if( ( tmpY[i]>tmpY[i-1] ) && ( tmpY[i]>tmpY[i+1] ) )
	      {
                maxX[maxCount] = tmpX[i];
                maxY[maxCount] = tmpY[i];
                maxCount++;
	      }
	  }

        // last elements
        if( tmpY[lastNonNanInd-1]<tmpY[lastNonNanInd] )
	  {
            maxX[maxCount] = tmpX[lastNonNanInd];
            maxY[maxCount] = tmpY[lastNonNanInd];
            maxCount++;
	  }

        //tmpX = vector<Float64>( maxCount );
        //tmpY = vector<Float64>( maxCount );

        tmpSize = maxCount;
        // Prepare for next iteration by storing every maximum found in this iteration in the tmp array used as input for the next iteration
        for( Int32 t=0; t<maxCount; t++ )
	  {
            tmpX[t]=maxX[t];
            tmpY[t]=maxY[t];
	  }

        /*//debug:
        // save median and xmad,  flux data
        FILE* f = fopen( "extremum_dbg.txt", "w+" );
        for( Int32 t=0;t<maxCount;t++)
        {
            fprintf( f, "%d %f %f\n", t, tmpX[t], tmpY[t]);
        }
        fclose( f );
        //*/

        if( maxCount == 0 )
	  {
            break;
	  }

        if( tmpSize <= PEAKS_SMOOTH_LIMIT )
	  {
            break;
	  }
    }

    Int32 nbPeaks = m_MaxPeakCount;
    if( tmpSize<m_MaxPeakCount )
      {
        nbPeaks = tmpSize;
      }

    // Extract best results
    CQuickSort<Float64> sort;

    vector<Int32> sortedIndexes( tmpSize );
    sort.SortIndexes( tmpY.data(), sortedIndexes.data(), sortedIndexes.size() );

    maxPoint.resize( nbPeaks );
    Int32 k = 0;

    for( Int32 i=0; i<nbPeaks; i++ )
      {
        Int32 j = sortedIndexes[(sortedIndexes.size()-1)-i];
        if( !std::isnan( tmpY[j] ) )
	  {
            maxPoint[k++] = SPoint( tmpX[j], m_SignSearch*tmpY[j] );
	  }
      }

    return true;
}

/**
 * Attempts to find peaks, returning (when appropriate) the points where they reside in the maxPoint argument.
 */
Bool CExtremum::InternalFind2( const Float64* xAxis, const Float64* yAxis, UInt32 n, TPointList& maxPoint ) const
{
    if( n == 0 )
        return false;

    vector<Float64> tmpX( n );
    vector<Float64> tmpY( n );
    Int32 tmpSize = n;
    for( UInt32 t=0; t<n; t++ )
      {
        if( !std::isnan( yAxis[t] ) )
	  {
            tmpX[t] = xAxis[t];
            tmpY[t] = yAxis[t];
	  }
      }

    // Extract best results
    CQuickSort<Float64> sort;

    vector<Int32> sortedIndexes( tmpSize );
    sort.SortIndexes( tmpY.data(), sortedIndexes.data(), tmpSize );

    maxPoint.resize( m_MaxPeakCount );
    Int32 k = 0;

    for( UInt32 i=0; i<m_MaxPeakCount; i++ )
      {
        UInt32 j = sortedIndexes[(sortedIndexes.size()-1)-i];

        maxPoint[k++] = SPoint( tmpX[j], tmpY[j] );
      }

    return true;
}
