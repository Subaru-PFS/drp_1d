#include <epic/redshift/operator/chicorr.h>

#include <epic/core/debug/assert.h>
#include <epic/core/log/log.h>
#include <epic/redshift/spectrum/fluxaxis.h>
#include <epic/redshift/spectrum/spectralaxis.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/tools.h>
#include <epic/redshift/common/mask.h>
#include <epic/redshift/operator/correlationresult.h>
#include <epic/redshift/operator/chisquareresult.h>
#include <epic/redshift/operator/operator.h>

#include <math.h>
#include <assert.h>


#include <boost/date_time/posix_time/posix_time.hpp>

#define NOT_OVERLAP_VALUE NAN
#include <stdio.h>

using namespace NSEpic;
using namespace std;

IMPLEMENT_MANAGED_OBJECT(COperatorChicorr)

COperatorChicorr::COperatorChicorr()
{
    m_TotalDuration = -1.0f;
}

COperatorChicorr::~COperatorChicorr()
{

}

Float64 COperatorChicorr::GetComputationDuration() const
{
    return m_TotalDuration;
}


/**

 */
int COperatorChicorr::Compute(const CSpectrum& spectrum, const CSpectrum& spectrumWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                                      const TFloat64Range& lambdaRange, const TFloat64List& redshifts,
                                      Float64 overlapThreshold, CCorrelationResult *result_corr, CChisquareResult *result_chi)
{
    Bool retVal;

    // Input spectrum MUST be in log scales
    if( spectrumWithoutCont.GetSpectralAxis().IsInLogScale() == false || tplWithoutCont.GetSpectralAxis().IsInLogScale() == false )
    {
        Log.LogError("Failed to compute Cross correlation / Fit, input spectrum or template are not in log scale");
        return -1;
    }

    boost::posix_time::ptime  startTime = boost::posix_time::microsec_clock::local_time();

    DebugAssert( overlapThreshold > 0.0 && overlapThreshold <= 1.0 );

    CSpectrumSpectralAxis shiftedTplSpectralAxis( tpl.GetSampleCount(), true );

    const CSpectrumSpectralAxis& spcSpectralAxis = spectrumWithoutCont.GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = spectrum.GetFluxAxis();
    const CSpectrumFluxAxis& spcWithoutContFluxAxis = spectrumWithoutCont.GetFluxAxis();

    const CSpectrumSpectralAxis& tplSpectralAxis = tplWithoutCont.GetSpectralAxis();
    const CSpectrumFluxAxis& tplFluxAxis = tpl.GetFluxAxis();
    const CSpectrumFluxAxis& tplWithoutContFluxAxis = tplWithoutCont.GetFluxAxis();

    CMask spcMask( spectrumWithoutCont.GetSampleCount() );
    CMask tplMask( tplWithoutCont.GetSampleCount() );

    // Compute clamped lambda range over spectrum
    TFloat64Range spcLambdaRange;
    retVal = spcSpectralAxis.ClampLambdaRange( lambdaRange, spcLambdaRange );
    DebugAssert( retVal );

    // Create Results bag
    result_corr->Correlation.resize( redshifts.size() );
    result_corr->Overlap.resize( redshifts.size() );
    result_corr->Redshifts.resize( redshifts.size() );
    result_corr->Status.resize( redshifts.size() );
    result_corr->Redshifts = redshifts;

    result_chi->ChiSquare.resize( redshifts.size() );
    result_chi->Redshifts.resize( redshifts.size() );
    result_chi->Overlap.resize( redshifts.size() );
    result_chi->Status.resize( redshifts.size() );
    result_chi->Redshifts = redshifts;

    for ( Int32 i=0; i<redshifts.size(); i++)
    {
        result_corr->Correlation[i] = NAN;
        result_corr->Status[i] = COperator::nStatus_DataError;
        result_corr->Overlap[i] = 0;

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
        if( TFloat64Range::Intersect( tplLambdaRange, spcLambdaRange, intersectedLambdaRange ) )
        {
            result_corr->Overlap[i] = intersectedLambdaRange.GetLength() / spcLambdaRange.GetLength();
            result_chi->Overlap[i] = result_corr->Overlap[i];
        }

        if( result_corr->Overlap[i] < overlapThreshold )
        {
            result_corr->Status[i] = COperator::nStatus_NoOverlap;
            continue;
        }

        const Float64* error = spcFluxAxis.GetError();

        DebugAssert( error!= NULL );

        Float64 spcMean = 0.0;
        Float64 spcSDev = 0.0;
        spcSpectralAxis.GetMask( intersectedLambdaRange, spcMask );
        if( !spcWithoutContFluxAxis.ComputeMeanAndSDev( spcMask, spcMean, spcSDev, error ) )
        {
            result_corr->Status[i] = COperator::nStatus_DataError;
            continue;
        }

        Float64 tplMean = 0.0;
        Float64 tplSDev = 0.0;
        shiftedTplSpectralAxis.GetMask( intersectedLambdaRange, tplMask );
        if( !tplWithoutContFluxAxis.ComputeMeanAndSDev( tplMask, tplMean, tplSDev, NULL ) )
        {
            result_corr->Status[i] = COperator::nStatus_DataError;
            continue;
        }

        const Float64* Xtpl = shiftedTplSpectralAxis.GetSamples();
        const Float64* YtplWithoutCont = tplWithoutContFluxAxis.GetSamples();
        const Float64* Ytpl = tplFluxAxis.GetSamples();
        const Float64* Xspc = spcSpectralAxis.GetSamples();
        const Float64* YspcWithoutCont = spcWithoutContFluxAxis.GetSamples();
        const Float64* Yspc = spcFluxAxis.GetSamples();

        TFloat64Range logIntersectedLambdaRange( log( intersectedLambdaRange.GetBegin() ), log( intersectedLambdaRange.GetEnd() ) );

        // Move cursors up to lambda range start
        Int32 j = 0;
        while( Xspc[j] < logIntersectedLambdaRange.GetBegin() )
            j++;

        Int32 k = 0;
        while( Xtpl[k] < logIntersectedLambdaRange.GetBegin() )
            k++;

        Int32 jStart = j;
        Int32 kStart = k;

        // xcorr vars
        Float64 sumCorr=0;
        Float64 sumWeight=0;
        // fit vars
        Float64 sumXDevs = 0.0;
        Float64 sumYDevs = 0.0;
        Float64 err2 = 0.0;
        Float64 fit = 0;
        Int32 numDevs = 0;


        Float64 t = 0;
        // k index: move over template
        // j index: move over spectrum
        vector<Float64> tplInterpolatedFlux;
        tplInterpolatedFlux.resize( spectrum.GetSampleCount() );
        vector<Float64> tplWithoutContInterpolatedFlux;
        tplWithoutContInterpolatedFlux.resize( spectrum.GetSampleCount() );
        Int64 iInterp = 0;

        // For each sample in the valid lambda range interval.
        while( k<tpl.GetSampleCount()-1 && Xtpl[k] <= logIntersectedLambdaRange.GetEnd() )
        {
            // For each sample in the spectrum that are in between two continous template sample
            while( j<spectrum.GetSampleCount() && Xspc[j] <= Xtpl[k+1] )
            {
                // perform linear interpolation of the flux
                t = ( Xspc[j] - Xtpl[k] ) / ( Xtpl[k+1] - Xtpl[k] );
                tplInterpolatedFlux[iInterp] = Ytpl[k] + ( Ytpl[k+1] - Ytpl[k] ) * t;
                tplWithoutContInterpolatedFlux[iInterp] = YtplWithoutCont[k] + ( YtplWithoutCont[k+1] - YtplWithoutCont[k] ) * t;

                numDevs++;
                err2 = 1.0 / error[j] * error[j];
                sumYDevs+=Yspc[j]*err2;
                sumXDevs+=tplInterpolatedFlux[iInterp]*err2;

                iInterp++;
                j++;
            }

            k++;
        }

        if ( numDevs==0 || sumYDevs==0 || sumXDevs==0 )
        {
            return -1;
        }

        Float64 ampl = sumYDevs / sumXDevs;

        fit=0;

        Float64 s = 0;

        for(j=jStart; j<spectrum.GetSampleCount(); j++){
            int k=j;
            {
                // fit
                fit += pow( Yspc[j] - ampl * tplInterpolatedFlux[k] , 2.0 ) / pow( error[j], 2.0 );
                s += Yspc[j];

                // xcorr
                sumWeight += 1.0 / error[j];
                sumCorr += ( tplWithoutContInterpolatedFlux[k] - tplMean ) * ( YspcWithoutCont[j] - spcMean ) / error[j];
            }
        }

        if( sumWeight>0 )
        {
            result_corr->Correlation[i] = sumCorr / ( tplSDev * spcSDev * sumWeight );
            result_corr->Status[i] = COperator::nStatus_OK;
        }


        // Chi square reduct: it can introduces some problem?
        fit /= numDevs;
        result_chi->ChiSquare[i] = fit;
        result_chi->Status[i] = COperator::nStatus_OK;
    }

    /*//debug:
    // save Correlation
    FILE* f = fopen( "full_dbg.txt", "w+" );
    for( Int32 i=0; i<redshifts.size(); i++ )
    {
        if( result_corr->Correlation[i] < 0.0001 ){
            fprintf( f, "%d %e %e\n", i, result_corr->Redshifts[i], result_corr->Correlation[i]);
        }else{
            fprintf( f, "%d %f %f\n", i, result_corr->Redshifts[i], result_corr->Correlation[i]);
        }
    }
    fclose( f );
    //*/


    //*//debug:
    // save Chisquare
    FILE* f = fopen( "full_chisquare_dbg.txt", "w+" );
    for( Int32 i=0; i<redshifts.size(); i++ )
    {
        if( result_chi->ChiSquare[i] < 0.0001 ){
            fprintf( f, "%d %e %e\n", i, result_chi->Redshifts[i], result_chi->ChiSquare[i]);
        }else{
            fprintf( f, "%d %f %f\n", i, result_chi->Redshifts[i], result_chi->ChiSquare[i]);
        }
    }
    fclose( f );
    //*/

    boost::posix_time::time_duration diff = boost::posix_time::microsec_clock::local_time() - startTime;

    m_TotalDuration = diff.total_seconds();

    return -1;
}



