#include <RedshiftLibrary/operator/correlation.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/spectrum/fluxaxis.h>
#include <RedshiftLibrary/spectrum/spectralaxis.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/tools.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/operator/correlationresult.h>

#include <math.h>
#include <assert.h>


#include <boost/date_time/posix_time/posix_time.hpp>

#define NOT_OVERLAP_VALUE NAN
#include <stdio.h>

using namespace NSEpic;
using namespace std;

COperatorCorrelation::COperatorCorrelation()
{
    m_TotalDuration = -1.0f;
}

COperatorCorrelation::~COperatorCorrelation()
{

}

Float64 COperatorCorrelation::GetComputationDuration() const
{
    return m_TotalDuration;
}

#include <fstream>

/**
 * Compute correlation factor between spectrum and tpl for each value specified in redshifts.
 * \note Spectrum AND template MUST be in log scale
 */
 std::shared_ptr<COperatorResult> COperatorCorrelation::Compute(const CSpectrum& spectrum,
                                                                   const CTemplate& tpl,
                                                                   const TFloat64Range& lambdaRange,
                                                                   const TFloat64List& redshifts,
                                                                   Float64 overlapThreshold,
                                                                   std::vector<CMask> additional_spcMasks_unused,
                                                                   std::string opt_interp_unused,
                                                                   Int32 opt_extinction_unused, Int32 opt_dustFitting_unused)
 {
    Bool retVal;

    // Input spectrum MUST be in log scales
    if( spectrum.GetSpectralAxis().IsInLogScale() == false || tpl.GetSpectralAxis().IsInLogScale() == false )
    {
        Log.LogError("Failed to compute Cross correlation, input spectrum or template are not in log scale");
        return NULL;
    }

    boost::posix_time::ptime  startTime = boost::posix_time::microsec_clock::local_time();

    DebugAssert( overlapThreshold > 0.0 && overlapThreshold <= 1.0 );

    CSpectrumSpectralAxis shiftedTplSpectralAxis( tpl.GetSampleCount(), true );

    const CSpectrumSpectralAxis& spcSpectralAxis = spectrum.GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = spectrum.GetFluxAxis();

    const CSpectrumSpectralAxis& tplSpectralAxis = tpl.GetSpectralAxis();
    const CSpectrumFluxAxis& tplFluxAxis = tpl.GetFluxAxis();


    // Compute clamped lambda range over spectrum
    TFloat64Range spcLambdaRange;
    retVal = spcSpectralAxis.ClampLambdaRange( lambdaRange, spcLambdaRange );
    DebugAssert( retVal );

    // Create Results bag
    std::shared_ptr<CCorrelationResult> result = std::shared_ptr<CCorrelationResult>( new CCorrelationResult() );
    result->Correlation.resize( redshifts.size() );
    result->Overlap.resize( redshifts.size() );
    result->Redshifts.resize( redshifts.size() );
    result->Status.resize( redshifts.size() );

    result->Redshifts = redshifts;

    for ( Int32 i=0; i<redshifts.size(); i++)
    {
        result->Correlation[i] = NAN;
        result->Status[i] = nStatus_DataError;
        result->Overlap[i] = 0;

        // Shift Template (Since template are created at Z=0)
        Float64 onePlusRedshift = 1.0 + redshifts[i];
        shiftedTplSpectralAxis.ShiftByWaveLength( tplSpectralAxis, onePlusRedshift, CSpectrumSpectralAxis::nShiftForward );

        TFloat64Range intersectedLambdaRange( 0.0, 0.0 );

        // Compute clamped lambda range over template
        TFloat64Range tplLambdaRange;
        retVal = shiftedTplSpectralAxis.ClampLambdaRange( lambdaRange, tplLambdaRange );
        DebugAssert( retVal );

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

        result->Overlap[i] = mask.CompouteOverlapRate( itplMask );
        if( result->Overlap[i] < overlapThreshold )
        {
            result->Status[i] = nStatus_NoOverlap;
            continue;
        }

        const Float64* error = spcFluxAxis.GetError();

        DebugAssert( error!= NULL );

        Float64 spcMean = 0.0;
        Float64 spcSDev = 0.0;

        if( !spcFluxAxis.ComputeMeanAndSDev( itplMask, spcMean, spcSDev, error ) )
        {
            result->Status[i] = nStatus_DataError;
            continue;
        }

        Float64 tplMean = 0.0;
        Float64 tplSDev = 0.0;
        if( !itplTplFluxAxis.ComputeMeanAndSDev( itplMask, tplMean, tplSDev, NULL ) )
        {
            result->Status[i] = nStatus_DataError;
            continue;
        }

        const Float64* Xtpl = itplTplSpectralAxis.GetSamples();
        const Float64* Ytpl = itplTplFluxAxis.GetSamples();
        const Float64* Xspc = spcSpectralAxis.GetSamples();
        const Float64* Yspc = spcFluxAxis.GetSamples();

        TFloat64Range logIntersectedLambdaRange( log( intersectedLambdaRange.GetBegin() ), log( intersectedLambdaRange.GetEnd() ) );

        // j cursor move over spectrum
        Int32 j = 0;
        while( j < spcSpectralAxis.GetSamplesCount() && Xspc[j] < logIntersectedLambdaRange.GetBegin() )
            j++;

        // j cursor move over template
        Int32 k = 0;
        while( k < itplTplSpectralAxis.GetSamplesCount() && Xtpl[k] < logIntersectedLambdaRange.GetBegin() )
            k++;


        Float64 sumCorr=0;
        Float64 sumWeight=0;
        Float64 tplInterpolatedFlux=-1;
        Float64 t = 0;
        // k index: move over template
        // j index: move over spectrum
        
        // Cross correlation function:
        // (t(x) - moy(t)) * (s(x) - moy(s))
        // ---------------------------------
        //          stddev(x) * stddev(s)
        // For each sample in the valid lambda range interval.
        while( j<spcSpectralAxis.GetSamplesCount() && Xspc[j] <= logIntersectedLambdaRange.GetEnd() )
        {
            sumWeight += 1.0 / error[j];
            sumCorr += ( Ytpl[k] - tplMean ) * ( Yspc[j] - spcMean ) / error[j];

            j++;
            k++;
        }

        if( sumWeight>0 )
        {
            result->Correlation[i] = sumCorr / ( tplSDev * spcSDev * sumWeight );
            result->Status[i] = nStatus_OK;
        }
    }

    boost::posix_time::time_duration diff = boost::posix_time::microsec_clock::local_time() - startTime;

    m_TotalDuration = diff.total_seconds();

    return result;
}


