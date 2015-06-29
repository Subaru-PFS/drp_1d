#include <epic/redshift/operator/peakdetection.h>

#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/fluxaxis.h>
#include <epic/redshift/spectrum/spectralaxis.h>
#include <epic/redshift/common/median.h>
#include <epic/redshift/operator/peakdetectionresult.h>

#include <math.h>

#include <stdio.h>

using namespace NSEpic;
IMPLEMENT_MANAGED_OBJECT(CPeakDetection)

CPeakDetection::CPeakDetection(Float64 windowSize, Float64 cut, UInt32 medianSmoothHalfWidth, UInt32 enlargeRate)
{
    m_winsize = windowSize;
    m_cut = cut;
    m_medianSmoothHalfWidth = medianSmoothHalfWidth;
    m_enlargeRate = enlargeRate;
}

CPeakDetection::~CPeakDetection()
{

}

const CPeakDetectionResult* CPeakDetection::Compute( const CSpectrum& spectrum, const TLambdaRange& lambdaRange)
{
    CPeakDetectionResult* result = new CPeakDetectionResult();

    const CSpectrumFluxAxis& fluxAxis = spectrum.GetFluxAxis();
    const CSpectrumSpectralAxis& spectralAxis = spectrum.GetSpectralAxis();

    CSpectrumFluxAxis smoothedFluxAxis = fluxAxis;

    if( m_medianSmoothHalfWidth  )
    {
        smoothedFluxAxis.ApplyMedianSmooth( m_medianSmoothHalfWidth );
    }

    TInt32RangeList peaksBorders;
    FindPossiblePeaks( smoothedFluxAxis, spectralAxis, peaksBorders );

    // No Peak detected, exit
    if( peaksBorders.size() == 0 )
    {
        return NULL;
    }

    RedefineBorders( peaksBorders, spectralAxis, smoothedFluxAxis, fluxAxis );

    result->PeakList = peaksBorders;
    TInt32RangeList peaksBordersEnlarged= peaksBorders;
    if( m_enlargeRate )
    {
        for( UInt32 i=0; i<peaksBorders.size(); i++ )
        {
            TInt32Range fitRange = FindGaussianFitStartAndStop( i, peaksBorders, m_enlargeRate, spectralAxis.GetSamplesCount() );
            peaksBordersEnlarged[i] = fitRange;
        }
    }
    result->EnlargedPeakList = peaksBordersEnlarged;

    return result;
}

TInt32Range CPeakDetection::FindGaussianFitStartAndStop( Int32 i, const TInt32RangeList& peaksBorders, UInt32 enlargeRate, Int32 len )
{
    Int32 fitStart = peaksBorders[i].GetBegin();
    Int32 fitStop = peaksBorders[i].GetEnd()+1;

    Float64 width = fitStop - fitStart ;
    fitStart = max( 0, fitStart - (int)(enlargeRate*width) );
    fitStop = min( len, fitStop + (int)(enlargeRate*width) );

    if( i>0 )
    {
        if( peaksBorders[i-1].GetEnd() > -1 )
            fitStart = max( peaksBorders[i-1].GetEnd(), fitStart );
    }

    if( i<peaksBorders.size()-1 )
    {
        if( peaksBorders[i+1].GetBegin() > -1 )
            fitStop = min( peaksBorders[i+1].GetBegin(), fitStop );
    }

    return TInt32Range( fitStart, fitStop );
}

/**
 * Not fully implemented.
 * Does not seam to do any vital things (CF ez.function.lines.EZELfind: 697)
 *
 * _find_borders8, (CF ez.function.lines.EZELfind: 365)
 * 1. find center = max_position in the current range
 * 2. if center (=max position) is  on the left or on the right -> disable this peak
 * 3. todo: use waves and fluxAxis for proper logging
 */
Void CPeakDetection::RedefineBorders( TInt32RangeList& peakList, const CSpectrumAxis& waves, const CSpectrumAxis& smoothFluxAxis, const CSpectrumAxis& fluxAxis )
{
    const Float64* smoothFluxData = smoothFluxAxis.GetSamples();

    Int32 nPeaksInitial = peakList.size();
    for( Int32 iPeak=nPeaksInitial-1; iPeak>=0; iPeak-- )
    {
        int centerPos=-1;
        Float64 centerVal = -1e12;
        // find position of the maximum on the smoothed flux
        Int32 start = peakList[iPeak].GetBegin();
        Int32 stop = peakList[iPeak].GetEnd() + 1;

        for( Int32 i= start; i< stop; i++ )
        {
            if(centerVal<smoothFluxData[i]){
                centerVal = smoothFluxData[i];
                centerPos = i;
            }
        }

        int old_left= centerPos;
        int old_right= centerPos;
        int new_left= std::max(0, centerPos-1);
        int new_right= std::min((int)smoothFluxAxis.GetSamplesCount(), centerPos+1);

        if(new_left==0 || new_right==smoothFluxAxis.GetSamplesCount()){
            // if the max is found on the border, then erase this range
            peakList.erase(peakList.begin() + iPeak);
        }else{
            while(smoothFluxData[new_left]<=smoothFluxData[old_left]){
                if(new_left==0){
                    old_left=0;
                    break;
                }
                old_left=new_left;
                new_left=max(0, new_left-1);
            }
            while(smoothFluxData[new_right]<=smoothFluxData[old_right]){
                if(new_right==fluxAxis.GetSamplesCount()-1){
                    old_right=new_right;
                    break;
                }
                old_right=new_right;
                new_right=max(0, new_right+1);
            }
            if(new_right == new_left){
                peakList.erase(peakList.begin() + iPeak);
            }
        }

    }
}


Void CPeakDetection::FindPossiblePeaks( const CSpectrumAxis& fluxAxis, const CSpectrumSpectralAxis& spectralAxis, TInt32RangeList& peakList )
{
    peakList.clear();

    std::vector<Float64> med;
    std::vector<Float64> xmad;

    med.reserve( fluxAxis.GetSamplesCount() );
    xmad.reserve( fluxAxis.GetSamplesCount() );

    // Compute median value for each sample over a window of size windowSampleCount
    const Float64* fluxData = fluxAxis.GetSamples();
    CMedian<Float64> medianFilter;
    for( Int32 i=0; i<fluxAxis.GetSamplesCount(); i++ )
    {
        //old: regular sampling hypthesis
        //Int32 halfWindowSampleCount = windowSampleCount / 2;
        //UInt32 start = std::max( 0, i - halfWindowSampleCount );
        //UInt32 stop = std::min( (Int32) fluxAxis.GetSamplesCount(), i + halfWindowSampleCount );
        //irregular sampling compatible
        UInt32 start = std::max(0, spectralAxis.GetIndexAtWaveLength(spectralAxis[i]-m_winsize/2.0) );
        UInt32 stop = std::min( (Int32) fluxAxis.GetSamplesCount(), spectralAxis.GetIndexAtWaveLength(spectralAxis[i]+m_winsize/2.0)  );

        med[i] = medianFilter.Find( fluxData + start, stop - start );
        xmad[i] = XMad( fluxData+ start, stop - start , med[i] );
    }

    /*//debug:
    // save median and xmad,  flux data
    FILE* f = fopen( "peakdetection_dbg_median.txt", "w+" );
    for( Int32 i=0; i<fluxAxis.GetSamplesCount(); i++ )
    {
        if( med[i] < 0.0001 ){
            fprintf( f, "%e %e %e %e %e\n", spectralAxis[i], med[i], xmad[i], med[i]+0.5*m_cut*xmad[i], fluxData[i]);
        }else{
            fprintf( f, "%f %f %f %f %f\n", spectralAxis[i], med[i], xmad[i], med[i]+0.5*m_cut*xmad[i], fluxData[i]);
        }
    }
    fclose( f );
    //*/


    // Detect each point whose value is over the median precomputed median
    std::vector<Bool> points;
    points.resize( fluxAxis.GetSamplesCount() + 1 );
    Int32 j = 0;

    for( Int32 i=0; i<fluxAxis.GetSamplesCount(); i++ )
    {
        if( fluxData[i] > med[i]+0.5*m_cut*xmad[i] )
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

    for( Int32 i=1; i<j+1; i++ )
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
        UInt32 i1 = n/2 - 1;
        UInt32 i2 = n/2;
        xmadm = 0.5*(xdata[i1]+xdata[i2]);
    }
    else
    {
        UInt32 i1 = int(n/2);
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
