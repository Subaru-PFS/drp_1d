#include <epic/redshift/extremum/extremum.h>

#include <epic/core/common/quicksort.h>
#include <epic/core/common/datatypes.h>

#include <math.h>
#include <float.h>
#include <iostream>

#define PEAKS_MIN_THRESHOLD (3)
#define PEAKS_SMOOTH_LIMIT (20)

using namespace NSEpic;
using namespace std;

CExtremum::CExtremum() :
    m_MaxPeakCount( 5 ),
    m_RefreshCount( 1 ),
    m_XRange( 0.0, 0.0 )
{

}

CExtremum::CExtremum( const TFloat64Range& xRange, UInt32 maxPeakCount, UInt32 refreshCount ) :
    m_MaxPeakCount( maxPeakCount ),
    m_RefreshCount( refreshCount ),
    m_XRange( xRange )
{

}

CExtremum::~CExtremum()
{

}

void CExtremum::SetXRange( const TFloat64Range& r )
{
    m_XRange = r;
}

void CExtremum::SetMaxPeakCount( UInt32 n )
{
    m_MaxPeakCount = n;
}

void CExtremum::SetRefreshCount( UInt32 n )
{
    m_RefreshCount = n;
}

Bool CExtremum::Find( const Float64* xAxis, const Float64* yAxis, UInt32 n, TPointList& maxPoint ) const
{
    const Float64* selectedXAxis = xAxis;
    const Float64* selectedYAxis = yAxis;

    Int32 rangeXBeginIndex = -1;
    Int32 rangeXEndIndex = -1;

    // Find index in xAxis that correspond to the boundary specified by m_XRange
    if( m_XRange.GetIsEmpty() )
    {
        rangeXBeginIndex = 0;
        rangeXEndIndex = n-1;
    }
    else
    {
    	// Find index range for the given lambda range
        for( Int32 i=0; i<n; i++ )
        {
            if( rangeXBeginIndex == -1 && xAxis[i] >= m_XRange.GetBegin() )
            {
                rangeXBeginIndex = i;
            }

            if( xAxis[i] <= m_XRange.GetEnd() )
            {
                rangeXEndIndex = i;
            }
        }
    }

    return InternalFind( selectedXAxis + rangeXBeginIndex, selectedYAxis + rangeXBeginIndex, ( rangeXEndIndex - rangeXBeginIndex ) + 1, maxPoint );
}

Bool CExtremum::InternalFind( const Float64* xAxis, const Float64* yAxis, UInt32 n, TPointList& maxPoint ) const
{
    if( n == 0 )
        return false;

    //Method 1, use only 1 extremum
    //*
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
    Int32 tmpSize = n;
    for( Int32 t=0;t<n;t++)
    {
        tmpX[t]=xAxis[t];
        tmpY[t]=yAxis[t];
    }

    for( Int32 count = 0; count < m_RefreshCount; count ++ )
    {
        Int32 maxCount = 0;

        // First element
        if( tmpY[0] > tmpY[1] )
        {
            maxX[maxCount] = tmpX[0];
            maxY[maxCount] = tmpY[0];
            maxCount++;
        }

        for( Int32 i=1; i<tmpSize-2; i++ )
        {
            if( ( tmpY[i] > tmpY[i-1] ) && ( tmpY[i] > tmpY[i+1] ) )
            {
                maxX[maxCount] = tmpX[i];
                maxY[maxCount] = tmpY[i];
                maxCount++;
            }
        }

        // last elements
        if( tmpY[tmpSize-2] < tmpY[tmpSize-1] )
        {
            maxX[maxCount] = tmpX[tmpSize-1];
            maxY[maxCount] = tmpY[tmpSize-1];
            maxCount++;
        }

        //tmpX = vector<Float64>( maxCount );
        //tmpY = vector<Float64>( maxCount );

        tmpSize = maxCount;
        // Prepare for next iteration by storing every maximum found in this iteration in the tmp array used as input for the next iteration
        for( Int32 t=0;t<maxCount;t++)
        {
            tmpX[t]=maxX[t];
            tmpY[t]=maxY[t];
        }

        if( maxCount == 0 )
            break;

        if( tmpSize <= PEAKS_SMOOTH_LIMIT )
            break;
    }

    Int32 nbPeaks = m_MaxPeakCount;
    if( tmpSize < m_MaxPeakCount )
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
        Int32 j = sortedIndexes[ (sortedIndexes.size()-1) - i];
        if( ! isnan( tmpY[j] ) )
        {
            maxPoint[ k++ ] = SPoint( tmpX[j], tmpY[j]);
        }
    }

    return true;
}

Bool CExtremum::InternalFind2( const Float64* xAxis, const Float64* yAxis, UInt32 n, TPointList& maxPoint ) const
{
    if( n == 0 )
        return false;

    vector<Float64> tmpX( n );
    vector<Float64> tmpY( n );
    Int32 tmpSize = n;
    for( Int32 t=0;t<n;t++)
    {
        if( ! isnan( yAxis[t] ) )
        {
            tmpX[t]=xAxis[t];
            tmpY[t]=yAxis[t];
        }
    }

    // Extract best results
    CQuickSort<Float64> sort;

    vector<Int32> sortedIndexes( tmpSize );
    sort.SortIndexes( tmpY.data(), sortedIndexes.data(), tmpSize );

    maxPoint.resize( m_MaxPeakCount );
    Int32 k = 0;

    for( Int32 i=0; i<m_MaxPeakCount; i++ )
    {
        Int32 j = sortedIndexes[ (sortedIndexes.size()-1) - i];

        maxPoint[ k++ ] = SPoint( tmpX[j], tmpY[j]);
    }

    return true;
}
