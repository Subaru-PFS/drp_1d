// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include "RedshiftLibrary/continuum/irregularsamplingmedian.h"

#include "RedshiftLibrary/common/median.h"
#include "RedshiftLibrary/common/mean.h"
#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

#include <math.h>
#include <algorithm>
#include <boost/algorithm/clamp.hpp>

using namespace NSEpic;
using namespace std;

/**
 * Sets the value for m_MeanSmoothAmplitude.
 */
void CContinuumIrregularSamplingMedian::SetMeanKernelWidth( Float32 width )
{
    m_MeanSmoothAmplitude = width;
}

/**
 * Sets the value for m_MedianSmoothAmplitude.
 */
void CContinuumIrregularSamplingMedian::SetMedianKernelWidth( Float32 width )
{
    m_MedianSmoothAmplitude = width;
}

/**
 * Sets m_MedianSmoothCycles to input.
 */
void CContinuumIrregularSamplingMedian::SetMedianCycleCount( Int32 count )
{
    m_MedianSmoothCycles = count;
}

/**
 * Sets m_even to input.
 */
void CContinuumIrregularSamplingMedian::SetMedianEvenReflection( bool evenReflection )
{
    m_MedianEvenReflection = evenReflection;
}

/**
 * Take an array lenght N,
 * foreach j element computes median value on interval j-n_range/2,j+n_range/2
 * an put this value in y_out array
 */
TFloat64List CContinuumIrregularSamplingMedian::MedianSmooth( const TFloat64List &y, Int32 n_range) const
{
    Int32 i;
    Int32 half,rest;
    Int32 start, stop;
    Int32 n_points = y.size();

    half = n_range/2;
    rest = n_range-2*half;    

    TFloat64List y_out(n_points);

    CMedian<Float64> median;
    for( i=0; i<n_points; i++ )
    {
        start = max( 0, i-half );
        stop = min( i+half+rest, n_points );

        y_out[i] = median.Find( y.begin()+start, y.begin()+stop );
    }

    return y_out;
}

/**
 * Compute an average smooth.
 */
TFloat64List CContinuumIrregularSamplingMedian::MeanSmooth( const TFloat64List &y, Int32 n) const
{
    Int32 i;
    Int32 start,end,half,rest;
    Int32 N = y.size();

    half = n/2;
    rest = n-2*half;

    TFloat64List y_out(N);

    CMean<Float64> mean;
    for( i=0; i<N; i++ )
    {
        start = max( 0, i-half );
        end = min( i+half+rest, N );
        
        y_out[i] = mean.Find(y.begin()+start, y.begin()+end);
    }
    return y_out;
 }

/**
 * Input: y_input, which has N elemnts
 * Reflects on border this array
 * y_out extended array
 */
TFloat64List CContinuumIrregularSamplingMedian::OddMirror(  const TFloat64List::const_iterator & begin, 
                                                            const TFloat64List::const_iterator & end,
                                                            Int32 Nreflex, Float64 y_input_begin_val, Float64 y_input_end_val) const
{
    TFloat64List y_out(std::distance(begin,end)+2*Nreflex);

    auto out_begin = y_out.begin()+Nreflex;
    std::copy(begin, end, out_begin);

    auto out_rbegin = std::reverse_iterator<TFloat64List::iterator>(out_begin);
    std::transform( begin, begin+Nreflex, out_rbegin, 
                    [y_input_begin_val](Float64 in){ 
                        return 2*y_input_begin_val - in;
                        });

    auto rbegin = std::reverse_iterator<TFloat64List::const_iterator>(end);
    out_begin = y_out.end()-Nreflex;
    std::transform(rbegin, rbegin+Nreflex, out_begin, [y_input_end_val](Float64 in){
        return 2*y_input_end_val - in;});

    return y_out;
}


/**
 * Extends array with its own "reflection". Version for even-sized arrays.
 */
TFloat64List CContinuumIrregularSamplingMedian::EvenMirror(  const TFloat64List::const_iterator & begin, 
                                                            const TFloat64List::const_iterator & end,
                                                            Int32 Nreflex) const
{
    TFloat64List y_out(std::distance(begin,end)+2*Nreflex);

    auto out_begin = y_out.begin()+Nreflex;
    std::copy(begin, end, out_begin);

    auto out_rbegin = std::reverse_iterator<TFloat64List::iterator>(out_begin);
    std::copy(begin, begin+Nreflex, out_rbegin );

    auto rbegin = std::reverse_iterator<TFloat64List::const_iterator>(end);
    out_begin = y_out.end()-Nreflex;
    std::copy(rbegin, rbegin+Nreflex, out_begin);

    return y_out;
}

//estimate lin fit border values for odd reflex robustness
Float64 CContinuumIrregularSamplingMedian::FitBorder(const CSpectrum& s, Int32 kstart, Int32 kend, bool isRightBorder) const
{
    const CSpectrumSpectralAxis& spectralAxis = s.GetSpectralAxis();
    const CSpectrumFluxAxis& fluxAxis = s.GetRawFluxAxis();

    TFloat64Range range = TFloat64Range(spectralAxis[kstart], spectralAxis[kend]);
    
    Float64 a;
    Float64 b;
    bool ret = s.GetLinearRegInRange( range,  a, b);
    Int32 k = isRightBorder ? kend : kstart;
    Float64 fitValue = ret ? spectralAxis[k]*a + b : fluxAxis[k];

    return fitValue;
}

/**
 * This algorithm computes the spectrum continuum using the following procedure:
 * estimate the 'middle' resolution (ex: for PFS the resolution doesn't vary more than ~ 25 % in a spectral axis')
 * compute the continuum using the median technique with that adapted resolution
**/
bool CContinuumIrregularSamplingMedian::RemoveContinuum( const CSpectrum& s, CSpectrumFluxAxis& noContinuumFluxAxis ) const
{
    Float64 resolution = s.GetMeanResolution();

    bool result = ProcessRemoveContinuum( s, noContinuumFluxAxis, resolution );

    return result;
}

/**
 * Computes the continuum using the median technique, and the input resolution.
 */
bool CContinuumIrregularSamplingMedian::ProcessRemoveContinuum( const CSpectrum& s, CSpectrumFluxAxis& noContinuumFluxAxis, Float64 resolution ) const
{
    Int32 k0 = 0;
    Int32 k1 = 0;
    Int32 nd;

    Int32 j, k;

    const CSpectrumFluxAxis& fluxAxis = s.GetRawFluxAxis();

    Int32 norig = s.GetSampleCount();

    Int32 meanSmoothAmplitude = round(m_MeanSmoothAmplitude/resolution);

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

    Int32 medianSmoothAmplitude = round(m_MedianSmoothAmplitude/resolution);

    // set default
    if( medianSmoothAmplitude<=0 )
    {
        return false;
    }

    medianSmoothAmplitude = max( meanSmoothAmplitude, medianSmoothAmplitude );

    noContinuumFluxAxis = CSpectrumFluxAxis(norig, 0.);
    noContinuumFluxAxis.GetError() = fluxAxis.GetError();
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

    //set the reflection size
    const Int32 nreflex = boost::algorithm::clamp(round(nd / 2.0), 10.0, 5.0 * meanSmoothAmplitude);

    // declare extend vector array
    TFloat64List ysmoobig;

    // reflect original "effective spectrum" and set it in ysmoobig
    // if m_Even==1 reflects spectrum as an m_Even function
    // if m_Even==0 reflacts spectrum as an odd function
    const TFloat64List & flux = fluxAxis.GetSamplesVector();
    if( m_MedianEvenReflection )
    {
        ysmoobig = EvenMirror( flux.begin()+k0, flux.begin()+k1+1, nreflex);
    }else{
        
        Float64 FBegin = FitBorder(s, k0, std::max(k0+nreflex,k1), false);

        Float64 FEnd = FitBorder(s, std::min(k0,k1-nreflex), k1, true);
        
        ysmoobig = OddMirror( flux.begin()+k0, flux.begin()+k1+1, nreflex, FBegin, FEnd);
    }

    {
        // WARNING!!!
        // medianSmoothAmplitude must be an odd number;
        // otherwise a sorted array will be sorted
        // by the median smooth iterations

        // median smoothing
        for( k=0; k<m_MedianSmoothCycles; k++ )
        {
            // median smooth size=2*medianSmoothAmplitude+1
            ysmoobig = MedianSmooth( ysmoobig, 2*medianSmoothAmplitude+1);
        }

        // medianSmoothAmplitude must be odd
        //medianSmoothAmplitude= 2 * ( medianSmoothAmplitude/2 ) + 1;
        medianSmoothAmplitude |= 1;

        for( k=0; k<m_MedianSmoothCycles; k++ )
        {
            // median smooth size=medianSmoothAmplitude
            ysmoobig = MedianSmooth( ysmoobig, medianSmoothAmplitude);
        }
    }


    {
        // mean smoothing
        // mean smooth size=meanSmoothAmplitude/4
        ysmoobig = MeanSmooth( ysmoobig, (Int32) meanSmoothAmplitude/4);
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
