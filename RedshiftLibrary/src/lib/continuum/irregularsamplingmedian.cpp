#include <RedshiftLibrary/continuum/irregularsamplingmedian.h>

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
 * Default attributions constructor.
 */
CContinuumIrregularSamplingMedian::CContinuumIrregularSamplingMedian()
{
    m_MeanSmoothAmplitude = 75;     // Angstrom
    m_MedianSmoothCycles = 5;
    m_MedianSmoothAmplitude = 75;   // Angstrom
    m_Even = true;
}

/**
 * Empty destructor.
 */
CContinuumIrregularSamplingMedian::~CContinuumIrregularSamplingMedian()
{

}

/**
 * Sets the value for m_MeanSmoothAmplitude.
 */
Void CContinuumIrregularSamplingMedian::SetMeanKernelWidth( Float32 width )
{
    m_MeanSmoothAmplitude = width;
}

/**
 * Sets the value for m_MedianSmoothAmplitude.
 */
Void CContinuumIrregularSamplingMedian::SetMedianKernelWidth( Float32 width )
{
    m_MedianSmoothAmplitude = width;
}

/**
 * Sets m_MedianSmoothCycles to input.
 */
Void CContinuumIrregularSamplingMedian::SetMedianCycleCount( UInt32 count )
{
    m_MedianSmoothCycles = count;
}

/**
 * Take an array lenght N,
 * foreach j element computes median value on interval j-n_range/2,j+n_range/2
 * an put this value in y_out array
 */
Int32 CContinuumIrregularSamplingMedian::MedianSmooth( const Float64 *y, Int32 n_points, Int32 n_range, Float64 *y_out)
{
    Int32 i;

    Int32 half,rest;

    Int32 start, stop;

    half = n_range/2;
    rest = n_range-2*half;

    CMedian<Float64> median;
    for( i=0; i<n_points; i++ )
    {
        start = max( 0, i-half );
        stop = min( i+half+rest, n_points-1 );

        *( y_out+i ) = median.Find( y+start, stop-start );
    }

    return 0;
}

/**
 * Input: y_input, which has N elemnts
 * Reflects on border this array
 * y_out extended array
 */
Int32 CContinuumIrregularSamplingMedian::OddMirror( const Float64* y_input, Int32 N, Int32 Nreflex, Float64 y_input_begin_val, Float64 y_input_end_val, Float64* y_out )
{
    Int32 j;

    for( j=0; j<N; j++ )
    {
        *( y_out+j+Nreflex ) = *( y_input+j );
    }

    for( j=0; j<Nreflex; j++ )
    {
        *( y_out+Nreflex-j-1 ) = 2*y_input_begin_val-*( y_input+j );
        *( y_out+N+Nreflex+j ) = 2*y_input_end_val-*( y_input+N-j-1 );
    }
    return 0;
}

/**
 * Extends array with its own "reflection". Version for even-sized arrays.
 */
Int32 CContinuumIrregularSamplingMedian::EvenMirror( const Float64* y_input, Int32 N, Int32 Nreflex, Float64* y_out )
{
    int j;

    for( j=0; j<N; j++ )
    {
        *( y_out+j+Nreflex ) = *( y_input+j );
    }
    for( j=0; j<Nreflex; j++ )
    {
        *( y_out+Nreflex-j-1 ) = *( y_input+j );
        *( y_out+N+Nreflex+j ) = *( y_input+N-j-1 );
    }
    return 0;
}

/**
 * This algorithm computes the spectrum continuum using the following procedure:
 * estimate the 'middle' resolution (ex: for PFS the resolution doesn't vary more than ~ 25 % in a spectral axis')
 * compute the continuum using the median technique with that adapted resolution
**/
Bool CContinuumIrregularSamplingMedian::RemoveContinuum( const CSpectrum& s, CSpectrumFluxAxis& noContinuumFluxAxis )
{
    Float64 resolution = s.GetMeanResolution();

    ProcessRemoveContinuum( s, noContinuumFluxAxis, resolution );

    return true;
}

/**
 * Computes the continuum using the median technique, and the input resolution.
 */
Bool CContinuumIrregularSamplingMedian::ProcessRemoveContinuum( const CSpectrum& s, CSpectrumFluxAxis& noContinuumFluxAxis, Float64 resolution )
{
    Int32 k0 = 0;
    Int32 k1 = 0;
    Int32 nd;

    Int32 nreflex, nbig;
    Int32 j, k;

    const CSpectrumFluxAxis& fluxAxis = s.GetFluxAxis();

    Int32 norig = s.GetSampleCount();

    Float64 frac = m_MeanSmoothAmplitude / resolution - floor( m_MeanSmoothAmplitude/resolution );

    Int32 meanSmoothAmplitude = (Int32) m_MeanSmoothAmplitude/resolution;

    if( frac>=0.5 )
    {
        meanSmoothAmplitude+=1;
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

    m_MedianSmoothAmplitude = max( meanSmoothAmplitude, m_MedianSmoothAmplitude );

    noContinuumFluxAxis.SetSize( norig );
    TFloat64List& noContinuumFluxAxisError = noContinuumFluxAxis.GetError();
    const TFloat64List& fluxAxisError = fluxAxis.GetError();
    for( j=0; j<norig; j++ )
    {
        noContinuumFluxAxis[j] = 0;
        // Also copy error
        noContinuumFluxAxisError[j] = fluxAxisError[j];
    }

    // Find the first not null element, and put its index in k0
    k=0;
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
        if ( k>10 )
            break;
    }

    if ( k<=10 )
    {
        return false;
    }

    nd = k1 - k0 + 1;


    // 0<=k0<=k1<=n_orig [k0,k1]="effective spectrum"
    // lenght("effective spectrum")=nd
    {
        //set the reflection size
        Float64 tmp =  nd / 2.0;

        if( 5.0 * meanSmoothAmplitude < tmp )
        {
            tmp = 5.0*meanSmoothAmplitude;
        }

        if( tmp<10.0 )
        {
            tmp=10.;
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
    {
        EvenMirror( fluxAxis.GetSamples()+k0, nd, nreflex, ysmoobig.data() );
    }else{
        //estimate lin fit border values for odd reflex robustness
        //begin
        const CSpectrumSpectralAxis& spectralAxis = s.GetSpectralAxis();
        Int32 kBeginSup = k0+nreflex;
        if(kBeginSup>k1){
            kBeginSup = k1;
        }
        TFloat64Range rangeBegin = TFloat64Range( spectralAxis[k0], spectralAxis[kBeginSup] );
        Float64 a;
        Float64 b;
        bool retBegin = s.GetLinearRegInRange( rangeBegin,  a, b);
        Float64 FBegin = *(fluxAxis.GetSamples()+k0);
        if(retBegin){
            Float64 wlBegin = spectralAxis[k0];
            FBegin = wlBegin*a + b;
        }
        //end
        Int32 kEndInf = k1-nreflex-1;
        if(kEndInf<k0){
            kEndInf = k0;
        }
        TFloat64Range rangeEnd = TFloat64Range( spectralAxis[kEndInf], spectralAxis[k1] );
        bool retEnd = s.GetLinearRegInRange( rangeEnd,  a, b);
        Float64 FEnd = *(fluxAxis.GetSamples()+k1);
        if(retEnd){
            Float64 wlEnd = spectralAxis[k1];
            FEnd = wlEnd*a + b;
        }

        OddMirror( fluxAxis.GetSamples()+k0, nd, nreflex, FBegin, FEnd, ysmoobig.data() );
    }
    /*//debug:
    // save reflex data
    FILE* f = fopen( "median_even_dbg.txt", "w+" );
    for( Int32 t=0;t<ysmoobig.size();t++)
    {
        fprintf( f, "%d %f\n", t, ysmoobig[t]*1e17);
    }
    fclose( f );
    //*/


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
            MedianSmooth( ysmoobig.data(), nbig, 2*m_MedianSmoothAmplitude+1,temp.data() );
            ysmoobig = temp;
        }

        // m_MedianSmoothAmplitude must be odd
        m_MedianSmoothAmplitude= 2 * ( m_MedianSmoothAmplitude/2 ) + 1;

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
        noContinuumFluxAxis[j] = fluxAxis[j] - ysmoobig[j-k0+nreflex];
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
Int32 CContinuumIrregularSamplingMedian::MeanSmooth( const Float64 *y, Int32 N, Int32 n, Float64 *y_out )
{
    Int32 i;
    Int32 start,end,half,rest;

    half = n/2;
    rest = n-2*half;

    CMean<Float64> mean;
    for( i=0; i<N; i++ )
    {
        start = max( 0, i-half-rest );
        end = min( i+half, N-1 );

        *(y_out+i) = mean.Find( y+start, (end-start)+1);
    }
    return 0;
 }
