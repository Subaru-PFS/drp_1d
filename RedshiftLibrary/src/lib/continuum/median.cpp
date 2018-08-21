#include <RedshiftLibrary/continuum/median.h>

#include <RedshiftLibrary/common/quicksort.h>
#include <RedshiftLibrary/common/median.h>
#include <RedshiftLibrary/common/mean.h>
#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/spectrum/spectrum.h>

#include <math.h>

#include <algorithm>

using namespace NSEpic;
using namespace std;

/**
 * Attribution constructor.
 */
CContinuumMedian::CContinuumMedian()
{
    m_MeanSmoothAmplitude = 75;
    m_MedianSmoothCycles = 5;
    m_MedianSmoothAmplitude = 75;
    m_Even = false;
}

/**
 * Empty destructor.
 */
CContinuumMedian::~CContinuumMedian()
{

}

/**
 * Sets m_MeanSmoothAmplitude to width.
 */
void CContinuumMedian::SetMeanKernelWidth( Float32 width )
{
    m_MeanSmoothAmplitude = width;
}

/**
 * Sets m_MedianSmoothAmplitude to width.
 */
void CContinuumMedian::SetMedianKernelWidth( Float32 width )
{
    m_MedianSmoothAmplitude = width;
}

/**
 * Sets m_MedianSmoothCycles to count.
 */
void CContinuumMedian::SetMedianCycleCount( UInt32 count )
{
    m_MedianSmoothCycles = count;
}

/**
 * Take an array lenght N,
 * for each j element computes median value on interval j-n_range/2,j+n_range/2
 * an put this value in y_out array
 */
Int32 CContinuumMedian::MedianSmooth( const Float64 *y, Int32 n_points, Int32 n_range, Float64 *y_out)
{
    Int32 i;
    Int32 half, rest;
    Int32 start, stop;

    half = n_range/2;
    rest = n_range-2*half;

    CMedian<Float64> median;
    for( i=0; i<n_points; i++)
    {
        start = max( 0, i-half );
        stop = min( i+half+rest, n_points-1 );

        *(y_out+i) = median.Find( y+start, stop-start );
    }

    return 0;
}

/*
 * Input: y_input, which has N elemnts
 * Reflects on border this array
 * y_out extended array
 */
Int32 CContinuumMedian::OddMirror( const Float64* y_input, Int32 N, Int32 Nreflex, Float64* y_out )
{

    Int32 j;

    for( j=0; j<N; j++ )
    {
        *(y_out+j+Nreflex) = *(y_input+j);
    }
    for( j=0; j<Nreflex; j++ )
    {
        *(y_out+Nreflex-j-1) = 2**(y_input+0)-*(y_input+j);
        *(y_out+N+Nreflex+j) = 2**(y_input+N-1)-*(y_input+N-j-1);
    }
    return 0;
}

/**
 * Extend spectrum with even number of indexes using a "reflexion" of the spectrum itself.
 */
Int32 CContinuumMedian::EvenMirror( const Float64* y_input, Int32 N, Int32 Nreflex, Float64* y_out )
{
    int j;

    for( j=0; j<N; j++ )
    {
        *(y_out+j+Nreflex) = *(y_input+j);
    }
    for( j=0; j<Nreflex; j++ )
    {
        *(y_out+Nreflex-j-1) = *(y_input+j);
        *(y_out+N+Nreflex+j) = *(y_input+N-j-1);
    }
    return 0;
}

/**
 * This algorithm computes the spectrum continuum using the following procedure:
 * The spetrum is enlarged outside its borders, reflecting the first and last m_MeanSmoothAmplitude angstroms with respect to the wavelengths.
 * the added parts are reflected with respect to the fluxes.
 * On the enlarged spectrum, it is applied:
 * - m_MedianSmoothCycles times a median-smoothing with an amplitude of 2 * m_MedianSmoothAmplitude
 * - m_MedianSmoothCycles times a median-smoothing with an amplitude of m_MedianSmoothAmplitude
 * A mean smoothing with m_MeanSmoothingAmplitude / 4 is applied
 * The added parts are cut away
**/
Bool CContinuumMedian::RemoveContinuum( const CSpectrum& s, CSpectrumFluxAxis& noContinuumFluxAxis )
{
    Int32 k0 = 0;
    Int32 k1 = 0;
    Int32 nd;

    Int32 nreflex, nbig;
    Int32 j, k;

    const CSpectrumFluxAxis& fluxAxis = s.GetFluxAxis();

    Int32 norig = s.GetSampleCount();

    Float64 resolution = s.GetResolution();
    Float64 frac = m_MeanSmoothAmplitude / resolution - floor( m_MeanSmoothAmplitude/resolution );

    Int32 meanSmoothAmplitude = (Int32) m_MeanSmoothAmplitude/resolution;

    if( frac>=0.5 )
    {
        meanSmoothAmplitude += 1;
    }

    // set default
    if( meanSmoothAmplitude<=0 )
    {
        return false;
    }

    meanSmoothAmplitude = min( meanSmoothAmplitude,norig/2 );

    // set default
    if( m_MedianSmoothCycles<=0 )
    {
        return false;
    }

    frac = m_MedianSmoothAmplitude/resolution-floor( m_MedianSmoothAmplitude/resolution );

    m_MedianSmoothAmplitude = (Int32) m_MedianSmoothAmplitude/resolution;

    if( frac>=0.5 )
    {
        m_MedianSmoothAmplitude += 1;
    }

    // set default
    if( m_MedianSmoothAmplitude<=0 )
    {
        return false;
    }

    m_MedianSmoothAmplitude = max( meanSmoothAmplitude,m_MedianSmoothAmplitude );

    noContinuumFluxAxis.SetSize( norig );
    TFloat64List& noContinuumFluxAxisError = noContinuumFluxAxis.GetError();
    const TFloat64List& fluxAxisError = fluxAxis.GetError();
    for( j=0; j<norig; j++)
    {
        noContinuumFluxAxis[j] = 0;
        // Also copy error
        noContinuumFluxAxisError[j] = fluxAxisError[j];
    }

    // Find the first not null element, and put its index in k0
    k = 0;
    for( j=0; j<norig; j++ )
    {
        if( fluxAxis[j]!=0 )
        {
            k0 = j;
            break;
        }
    }

    // Find the last not null element, and put its index in k1
    for( j=norig-1; j>=0; j-- )
    {
        if( fluxAxis[j] != 0 )
        {
            k1 = j;
            break;
        }
    }

    // Strange test Here, ask MARCO why it's there
    k = 0;
    for( j=0; j<norig; j++ )
    {
        if( fluxAxis[j] != 0 )
        {
            k++;
        }
        if( k>10 )
            break;
    }

    if( k<=10 )
    {
        return false;
    }

    nd = k1 - k0 + 1;


    // 0<=k0<=k1<=n_orig [k0,k1]="effective spectrum"
    // lenght("effective spectrum")=nd
    {
        //set the reflection size
        Float64 tmp =  nd / 2.0;

        if( 5.0*meanSmoothAmplitude<tmp )
        {
            tmp = 5.0*meanSmoothAmplitude;
        }

        if( tmp<10.0 )
        {
            tmp = 10.;
        }

        nreflex = (Int32) tmp;
    }


    // spectrum reflected size
    nbig = nd+2*nreflex;

    // Allocate array
    vector<Float64> ysmoobig( nbig );

    // reflect original "effective spectrum" a set it in ysmoobig
    // if m_Even==1 reflects spectrum as an m_Even function
    // if m_Even==0 reflacts spectrum as an odd function
    if( m_Even )
        EvenMirror( fluxAxis.GetSamples()+k0, nd, nreflex, ysmoobig.data() );
    else
        OddMirror( fluxAxis.GetSamples()+k0, nd, nreflex, ysmoobig.data() );


    {
        // WARNING!!!
        // m_MedianSmoothAmplitude must be an odd number;
        // otherwise a sorted array will be sorted
        // by the median smooth iterations

        // median smoothing
        vector<Float64> temp( nbig );

        for( k=0; k<m_MedianSmoothCycles; k++ )
        {
            // median smooth size=2*m_MedianSmoothAmplitude+1
            MedianSmooth( ysmoobig.data(), nbig, 2*m_MedianSmoothAmplitude+1, temp.data() );
            ysmoobig = temp;
        }

        // m_MedianSmoothAmplitude must be odd
        m_MedianSmoothAmplitude = 2*(m_MedianSmoothAmplitude/2)+1;

        for( k=0; k<m_MedianSmoothCycles; k++ )
        {
            // median smooth size=m_MedianSmoothAmplitude
            MedianSmooth( ysmoobig.data(), nbig, m_MedianSmoothAmplitude, temp.data() );
            ysmoobig = temp;
        }
    }


    {
        // mean smoothing
        vector<Float64> temp( nbig );

        // mean smooth size=meanSmoothAmplitude/4
        MeanSmooth( ysmoobig.data(), nbig, (Int32) meanSmoothAmplitude/4, temp.data() );

        ysmoobig = temp;
    }

    // Copy spectrum before k0
    for( j=0; j<k0; j++ )
    {
        noContinuumFluxAxis[j] = fluxAxis[j];
    }

    // Set continuum inside "effective spectrum"
    for( j=k0; j<=k1; j++ )
    {
        noContinuumFluxAxis[j] = fluxAxis[j]-ysmoobig[j-k0+nreflex];
    }

    // Copy spectrum after k1
    if( k1+1<s.GetSampleCount() )
    {
        for( j=k1+1; j<s.GetSampleCount(); j++ )
        {
            noContinuumFluxAxis[j] = fluxAxis[j];
        }
    }

    return true;
}

/**
 * Compute an average smooth.
 */
Int32 CContinuumMedian::MeanSmooth( const Float64 *y, Int32 N, Int32 n, Float64 *y_out)
{
    Int32 i;
    Int32 start,end,half,rest;

    half = n/2;
    rest = n-2*half;

    CMean<Float64> mean;
    for( i=0; i<N; i++)
    {
        start = max( 0, i-half-rest );
        end = min( i+half, N-1 );

        *(y_out+i) = mean.Find( y+start, (end-start)+1);
    }
    return 0;
 }
