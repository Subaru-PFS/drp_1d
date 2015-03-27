#include <epic/redshift/peak/detection.h>

#include <epic/redshift/peak/peak.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/fluxaxis.h>
#include <epic/redshift/spectrum/spectralaxis.h>
#include <epic/redshift/common/median.h>

#include <math.h>

using namespace NSEpic;
IMPLEMENT_MANAGED_OBJECT(CPeakDetection)

CPeakDetection::CPeakDetection()
{

}

CPeakDetection::~CPeakDetection()
{

}

Bool CPeakDetection::Compute( const CSpectrum& spectrum, const TFloat64Range& lambdaRange )
{
    Bool mStarCheck = false;
    Float64 windowSize = 1.0f;
    Float64 cut = 1.0;

    if( mStarCheck )
    {

    }

    const CSpectrumFluxAxis& fluxAxis = spectrum.GetFluxAxis();
    const CSpectrumSpectralAxis& spectralAxis = spectrum.GetSpectralAxis();

    CSpectrumFluxAxis smoothedFluxAxis = fluxAxis;
    smoothedFluxAxis.ApplyMedianSmooth( 1 );

    UInt32 windowSampleCount = windowSize / spectrum.GetResolution();

    TInt32RangeList peaksBorders;
    FindPossiblePeaks( smoothedFluxAxis, spectralAxis, windowSampleCount, cut, peaksBorders );

    // No Peak detected, exit
    if( peaksBorders.size() == 0 )
    {
        return false;
    }

    RedefineBorders( peaksBorders, spectralAxis, smoothedFluxAxis, fluxAxis );

    for( UInt32 i=0; i<peaksBorders.size(); i++ )
    {
        TInt32Range fitRange = FindGaussianFitStartAndStop( i, peaksBorders, spectralAxis.GetSamplesCount() );
    }

    return true;
}

TInt32Range CPeakDetection::FindGaussianFitStartAndStop( Int32 i, const TInt32RangeList& peaksBorders, Int32 len )
{
    Int32 fitStart = peaksBorders[i].GetBegin();
    Int32 fitStop = peaksBorders[i].GetEnd();

    Int32 rate = fitStop - fitStart;
    fitStart = max( 0, fitStart - 2*rate );
    fitStop = min( len, fitStop + 2*rate );

    if( i>0 )
    {
        if( peaksBorders[i-1].GetEnd() > -1 )
            fitStart = max( peaksBorders[i-1].GetEnd(), fitStart );
    }

    if( i<peaksBorders.size() )
    {
        if( peaksBorders[i+1].GetBegin() > -1 )
            fitStop = max( peaksBorders[i+1].GetBegin(), fitStop );
    }

    return TInt32Range( fitStart, fitStop );
}

/**
 * Not yet implemented.
 * Does not seam to do any vital things (CF ez.function.lines.EZELfind: 697)
 */
Void CPeakDetection::RedefineBorders( TInt32RangeList& peakList, const CSpectrumAxis& waves, const CSpectrumAxis& smoothFluxAxis, const CSpectrumAxis& fluxAxis )
{

}


Void CPeakDetection::FindPossiblePeaks( const CSpectrumAxis& fluxAxis, const CSpectrumAxis& spectralAxis, UInt32 windowSampleCount, Float64 cut, TInt32RangeList& peakList )
{
    peakList.clear();

    std::vector<Float64> med;
    std::vector<Float64> xmad;

    med.reserve( fluxAxis.GetSamplesCount() );
    xmad.reserve( fluxAxis.GetSamplesCount() );

    // Compute median value for each sample over a window of size windowSampleCount
    const Float64* fluxData = fluxAxis.GetSamples();
    Int32 halfWindowSampleCount = windowSampleCount / 2;
    CMedian<Float64> medianFilter;
    for( Int32 i=0; i<fluxAxis.GetSamplesCount(); i++ )
    {
        UInt32 start = std::max( 0, i - halfWindowSampleCount );
        UInt32 stop = std::min( (Int32) fluxAxis.GetSamplesCount(), i + halfWindowSampleCount );

        med[i] = medianFilter.Find( fluxData + start, halfWindowSampleCount*2 + 1 );
        xmad[i] = XMad( fluxData+ start, halfWindowSampleCount*2 + 1, med[i] );
    }

    // Detect each point whose value is over the median precomputed median
    std::vector<Bool> points;
    points.resize( fluxAxis.GetSamplesCount() );
    Int32 j = 0;

    for( Int32 i=0; i<fluxAxis.GetSamplesCount(); i++ )
    {
        if( fluxData[i] > med[i]+0.5*cut*xmad[i] )
        {
            points[j++] = i;
        }
    }

    // No potential peak detected , exit
    if( j == 0 )
    {
        return;
    }

    // Find contiguous sample
    Int32 start = points[0];
    Int32 current = points[0];
    Int32 stop = points[0];

    for( Int32 i=1; i<j; i++ )
    {
        // Index stored in points[] are contiguous (i.e: 4,5,6)
        if( points[i]-1 == current )
        {
            current = points[i];
        }
        // if we hit a discontinuity, store the previous range of contiguous points representing a potential peak
        else
        {
            stop = current;
            peakList.push_back( TInt32Range(start, stop ) );
            start = points[i];
            current = start;
        }
    }
}

Float64 CPeakDetection::XMad( const Float64* x, Int32 n, Float64 median )
{
    std::vector<Float64> xdata;
    Float64 xmadm = 0.0;

    xdata.reserve( n );

    for( Int32 i=0;i<n; i++ )
    {
        xdata[i] = fabs( x[i]-median );
    }

    CQuickSort<Float64> sort;

    sort.Sort( xdata.data(), n);

    if( ((float)n)/2. - int(n/2.) == 0 )
    {
        UInt32 i1 = n/2;
        UInt32 i2 = n/2 + 1;
        xmadm = 0.5*(xdata[i1]+xdata[i2]);
    }
    else
    {
        UInt32 i1 = int(n/2) + 1;
        xmadm = xdata[i1];
    }

    return xmadm;
}

/*
def Xmad(data, xmed):
    """
    The XMAD subroutine calculates the Median Absolute Deviation from
    the sample median. The median, M , is subtracted from each
    ORDERED statistic and then the absolute value is taken. This new
    set of of statistics is then resorted so that they are ORDERED
    statistics. The MAD is then defined to be the median of this
    new set of statistics and is returned as XMADM. The MAD can
    be defined:

                   XMADM = median{ abs(x(i) - M) }

    where the x(i) are the values passed in the array XDATA, and
    the median, M, is passed in the array XLETTER. The set of stats
    in the brackets is assumed to be resorted. For more information
    see page 408 in UREDA.
    """
    ndat = len(data)
    xdata = []
    for item in data:
        xdata.append(abs(item-xmed))
    xdata.sort()
    xdata.insert(0, 0.)
    if (float(ndat)/2. - int(ndat/2.)) == 0:
        i1 = ndat/2
        i2 = ndat/2 + 1
        xmadm = 0.5*(xdata[i1]+xdata[i2])
    else:
        i1 = int(ndat/2) + 1
        xmadm = xdata[i1]
    return xmadm
*/

const TFloat64List& CPeakDetection::GetResults() const
{
    return m_Results;
}


