#include <epic/redshift/operator/chisquare.h>

#include <epic/redshift/spectrum/axis.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/tools.h>
#include <epic/redshift/common/mask.h>
#include <epic/redshift/common/redshifts.h>

#include <epic/core/log/log.h>

#include <math.h>
#include <assert.h>

#define NOT_OVERLAP_VALUE NAN
#include <stdio.h>

using namespace NSEpic;
using namespace std;

COperatorChiSquare::COperatorChiSquare()
{

}

COperatorChiSquare::~COperatorChiSquare()
{

}

const TFloat64List& COperatorChiSquare::GetResults() const
{
    return m_Chisquare;
}

Float64 COperatorChiSquare::BasicFit( const CSpectrum& spectrum, const CTemplate& tpl,
                                const TFloat64Range& lambdaRange, Float64 redshift, Float64 overlapThreshold, Float64& overlapRate )
{

    if( spectrum.GetSpectralAxis().IsInLinearScale() == false || tpl.GetSpectralAxis().IsInLinearScale() == false )
    {
        Log.LogError("Failed to compute Chi Square, input spetrum or template are not in linear scale");
        return false;
    }

    Bool retVal;
    CSpectrumSpectralAxis shiftedTplSpectralAxis( tpl.GetSampleCount(), false );

    const CSpectrumSpectralAxis& spcSpectralAxis = spectrum.GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = spectrum.GetFluxAxis();

    const CSpectrumSpectralAxis& tplSpectralAxis = tpl.GetSpectralAxis();
    const CSpectrumFluxAxis& tplFluxAxis = tpl.GetFluxAxis();

    // Compute clamped lambda range over spectrum
    TFloat64Range spcLambdaRange;
    retVal = spcSpectralAxis.ClampLambdaRange( lambdaRange, spcLambdaRange );

    // Compute shifted template
    Float64 onePlusRedshift = 1.0 + redshift;
    shiftedTplSpectralAxis.ShiftByWaveLength( tplSpectralAxis, onePlusRedshift, CSpectrumSpectralAxis::nShiftForward );

    TFloat64Range intersectedLambdaRange( 0.0, 0.0 );

    // Compute clamped lambda range over template
    TFloat64Range tplLambdaRange;
    retVal = shiftedTplSpectralAxis.ClampLambdaRange( lambdaRange, tplLambdaRange );

    overlapRate = 0;

    // if there is any intersection between the lambda range of the spectrum and the lambda range of the template
    // Compute the intersected range
    if( TFloat64Range::Intersect( tplLambdaRange, spcLambdaRange, intersectedLambdaRange ) )
    {
        overlapRate = intersectedLambdaRange.GetLength() / spcLambdaRange.GetLength();
    }

    // Check for overlap rate
    if( overlapRate < overlapThreshold )
        return false;

    const Float64* Xtpl = shiftedTplSpectralAxis.GetSamples();
    const Float64* Ytpl = tplFluxAxis.GetSamples();
    const Float64* Xspc = spcSpectralAxis.GetSamples();
    const Float64* Yspc = spcFluxAxis.GetSamples();

    // Move cursors up to lambda range start
    Int32 j = 0;
    while( Xspc[j] < intersectedLambdaRange.GetBegin() )
        j++;

    Int32 k = 0;
    while( Xtpl[k] < intersectedLambdaRange.GetBegin() )
        k++;

    Int32 jStart = j;
    Int32 kStart = k;

    Float64 sumXDevs = 0.0;
    Float64 sumYDevs = 0.0;
    Float64 err2 = 0.0;
    Float64 fit = 0;
    Float64 tplInterpolatedFlux=-1;
    Float32 t = 0;
    Int32 numDevs = 0;
    const Float64* error = spcFluxAxis.GetError();

    while( k<tpl.GetSampleCount()-1 && Xtpl[k] <= intersectedLambdaRange.GetEnd() )
    {
        while( j<spectrum.GetSampleCount() && Xspc[j] <= Xtpl[k+1] )
        {
            t = ( Xspc[j] - Xtpl[k] ) / ( Xtpl[k+1] - Xtpl[k] );

            tplInterpolatedFlux = Ytpl[k] + ( Ytpl[k+1] - Ytpl[k] ) * t;

            numDevs++;
            err2 = 1.0 / error[j] * error[j];
            sumYDevs+=Yspc[j]*err2;
            sumXDevs+=tplInterpolatedFlux*err2;

            j++;
        }

        k++;
    }

    if ( numDevs==0 || sumYDevs==0 || sumXDevs==0 )
    {
        return false;
    }

    Float64 ampl = sumYDevs / sumXDevs;

    j = jStart;
    k = kStart;

    fit=0;

    while( k<tpl.GetSampleCount()-1 && Xtpl[k] <= intersectedLambdaRange.GetEnd() )
    {
        while( j<spectrum.GetSampleCount() && Xspc[j] <= Xtpl[k+1] )
        {
            t = ( Xspc[j] - Xtpl[k] ) / ( Xtpl[k+1] - Xtpl[k] );
            tplInterpolatedFlux = Ytpl[k] + ( Ytpl[k+1] - Ytpl[k] ) * t;

            fit += pow( Yspc[j] - ampl * tplInterpolatedFlux , 2.0 ) / pow( error[j], 2.0 );

            j++;
        }

        k++;
    }

    // Chi square reduct: it can introduces some problem?
    fit /= numDevs;

    return fit;
}


Bool COperatorChiSquare::Compute(const CSpectrum& spectrum, const CTemplate& tpl,
                          const TFloat64Range& lambdaRange, const CRedshifts& redshifts,
                          Float64 overlapThreshold )
{

    Int32 i;
    char fit_computation=1;

    Float64 Ampl;
    Int32 num_devs;

    m_Chisquare.resize( redshifts.GetRedshiftsCount() );
    m_Overlap.resize( redshifts.GetRedshiftsCount() );

    for (i=0;i<redshifts.GetRedshiftsCount();i++)
    {
        m_Chisquare[i] = BasicFit( spectrum, tpl, lambdaRange, redshifts[i], overlapThreshold, m_Overlap[i] );

    }


    return true;
}
