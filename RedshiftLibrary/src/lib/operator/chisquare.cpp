#include <RedshiftLibrary/operator/chisquare.h>

#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/tools.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/operator/chisquareresult.h>

#include <RedshiftLibrary/log/log.h>

#include <boost/numeric/conversion/bounds.hpp>

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


Void COperatorChiSquare::BasicFit( const CSpectrum& spectrum, const CTemplate& tpl,
                                const TFloat64Range& lambdaRange, Float64 redshift, Float64 overlapThreshold,
                                Float64& overlapRate, Float64& chiSquare, Float64& fitamplitude, EStatus& status )
{
    chiSquare = boost::numeric::bounds<float>::highest();
    fitamplitude = -1.0;
    overlapRate = 0.0;
    status = nStatus_DataError;

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

    // if there is any intersection between the lambda range of the spectrum and the lambda range of the template
    // Compute the intersected range
    TFloat64Range::Intersect( tplLambdaRange, spcLambdaRange, intersectedLambdaRange );

    CSpectrumFluxAxis itplTplFluxAxis;
    CSpectrumSpectralAxis itplTplSpectralAxis;
    CMask itplMask;
    CSpectrumFluxAxis::Rebin( intersectedLambdaRange, tplFluxAxis, shiftedTplSpectralAxis, spcSpectralAxis, itplTplFluxAxis, itplTplSpectralAxis, itplMask );

    CMask mask;
    spcSpectralAxis.GetMask( lambdaRange, mask );
    itplMask &= mask;
    overlapRate = mask.CompouteOverlapRate( itplMask );
    // Check for overlap rate
    if( overlapRate < overlapThreshold )
    {
        status = nStatus_NoOverlap;
        return ;
    }

    const Float64* Xtpl = itplTplSpectralAxis.GetSamples();
    const Float64* Ytpl = itplTplFluxAxis.GetSamples();
    const Float64* Xspc = spcSpectralAxis.GetSamples();
    const Float64* Yspc = spcFluxAxis.GetSamples();
    TFloat64Range logIntersectedLambdaRange( log( intersectedLambdaRange.GetBegin() ), log( intersectedLambdaRange.GetEnd() ) );
    //the spectral axis should be in the same scale
    TFloat64Range currentRange = logIntersectedLambdaRange;
    if( spcSpectralAxis.IsInLinearScale() != tplSpectralAxis.IsInLinearScale() )
        return;
    if(spcSpectralAxis.IsInLinearScale()){
        currentRange = intersectedLambdaRange;
    }

    // j cursor move over spectrum
    Int32 j = 0;
    while( j < spcSpectralAxis.GetSamplesCount() && Xspc[j] < currentRange.GetBegin() )
        j++;

    // k cursor move over template
    Int32 k = 0;
    while( k < itplTplSpectralAxis.GetSamplesCount() && Xtpl[k] < currentRange.GetBegin() )
        k++;

    Int32 jStart = j;
    Int32 kStart = k;

    Float64 sumXDevs = 0.0;
    Float64 sumYDevs = 0.0;
    Float64 err2 = 0.0;
    Float64 fit = 0;
    Int32 numDevs = 0;
    const Float64* error = spcFluxAxis.GetError();

    while( j<spcSpectralAxis.GetSamplesCount() && Xspc[j] <= currentRange.GetEnd() )
    {
        numDevs++;
        err2 = 1.0 / (error[j] * error[j]);
        sumYDevs+=Yspc[j]*err2;
        sumXDevs+=Ytpl[j]*err2;

        j++;
    }

    if ( numDevs==0 || sumYDevs==0 || sumXDevs==0 )
    {
        status = nStatus_DataError;
        return;
    }

    Float64 ampl = sumYDevs / sumXDevs;

    j = jStart;
    k = kStart;

    fit=0;

    Float64 s = 0;

    while( j<spcSpectralAxis.GetSamplesCount() && Xspc[j] <= currentRange.GetEnd() )
    {
        int k=j;
        {
            // fit
            fit += pow( Yspc[j] - ampl * Ytpl[k] , 2.0 ) / pow( error[j], 2.0 );
            s += Yspc[j];
        }
        j++;
    }

    // Chi square reduct: it can introduces some problem?
    fit /= numDevs;


    chiSquare = fit;
    fitamplitude = ampl;
    status = nStatus_OK;
}


 std::shared_ptr<COperatorResult>  COperatorChiSquare::Compute(const CSpectrum& spectrum, const CTemplate& tpl,
                          const TFloat64Range& lambdaRange, const TFloat64List& redshifts,
                          Float64 overlapThreshold , std::vector<CMask> additional_spcMasks_unused, string opt_interp_unused, Int32 opt_extinction_unused, Int32 opt_dustFitting_unused)
{

    if( spectrum.GetSpectralAxis().IsInLinearScale() == false || tpl.GetSpectralAxis().IsInLinearScale() == false )
    {
        Log.LogError("Failed to compute Cross correlation, input spectrum or template are not in log scale");
        return NULL;
    }


     std::shared_ptr<CChisquareResult>  result = std::shared_ptr<CChisquareResult>( new CChisquareResult() );

    result->ChiSquare.resize( redshifts.size() );
    result->FitAmplitude.resize( redshifts.size() );
    result->Redshifts.resize( redshifts.size() );
    result->Overlap.resize( redshifts.size() );
    result->Status.resize( redshifts.size() );

    result->Redshifts = redshifts;


    for (Int32 i=0;i<redshifts.size();i++)
    {
        BasicFit( spectrum, tpl, lambdaRange, result->Redshifts[i], overlapThreshold, result->Overlap[i], result->ChiSquare[i], result->FitAmplitude[i], result->Status[i] );
    }


    return result;

}
