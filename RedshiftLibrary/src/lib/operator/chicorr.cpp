#include <RedshiftLibrary/operator/chicorr.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/spectrum/fluxaxis.h>
#include <RedshiftLibrary/spectrum/spectralaxis.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/tools.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/operator/correlationresult.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <RedshiftLibrary/operator/operator.h>

#include <math.h>
#include <assert.h>


#include <boost/date_time/posix_time/posix_time.hpp>

#define NOT_OVERLAP_VALUE NAN
#include <stdio.h>

using namespace NSEpic;
using namespace std;


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

    const CSpectrumSpectralAxis& spcSpectralAxis = spectrumWithoutCont.GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = spectrum.GetFluxAxis();
    const CSpectrumFluxAxis& spcWithoutContFluxAxis = spectrumWithoutCont.GetFluxAxis();

    const CSpectrumSpectralAxis& tplSpectralAxis = tplWithoutCont.GetSpectralAxis();
    const CSpectrumFluxAxis tplFluxAxis = tpl.GetFluxAxis();
    const CSpectrumFluxAxis tplWithoutContFluxAxis = tplWithoutCont.GetFluxAxis();


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

    // Pre-Allocate the rebined template and mask with regard to the spectrum size
    m_templateRebined_bf.GetSpectralAxis().SetSize(spectrum.GetSampleCount());
    m_templateRebined_bf.GetFluxAxis().SetSize(spectrum.GetSampleCount());
    m_templateWOContRebined_bf.GetSpectralAxis().SetSize(spectrum.GetSampleCount());
    m_templateWOContRebined_bf.GetFluxAxis().SetSize(spectrum.GetSampleCount());
    m_mskRebined_bf.SetSize(spectrum.GetSampleCount());
    m_spcSpectralAxis_restframe.SetSize(spectrum.GetSampleCount());

    for ( Int32 i=0; i<redshifts.size(); i++)
    {
        result_corr->Correlation[i] = NAN;
        result_corr->Status[i] = COperator::nStatus_DataError;
        result_corr->Overlap[i] = 0;

        Float64 onePlusRedshift = 1.0 + redshifts[i];

        //shift lambdaRange backward to be in restframe
        TFloat64Range lambdaRange_restframe( lambdaRange.GetBegin() / onePlusRedshift, 
                                            lambdaRange.GetEnd() / onePlusRedshift );

        //get spectrum in restframe
        m_spcSpectralAxis_restframe.ShiftByWaveLength(spcSpectralAxis, onePlusRedshift, CSpectrumSpectralAxis::nShiftBackward);
        TFloat64Range spcLambdaRange_restframe(spcLambdaRange.GetBegin()/onePlusRedshift,
                                               spcLambdaRange.GetEnd()/onePlusRedshift);

        // Compute clamped lambda range over template in restframe
        TFloat64Range tplLambdaRange;
        retVal = tplSpectralAxis.ClampLambdaRange( lambdaRange_restframe, tplLambdaRange );
        DebugAssert( retVal );

        // if there is any intersection between the lambda range of the spectrum and the lambda range of the template
        // Compute the intersected range   
        TFloat64Range intersectedLambdaRange( 0.0, 0.0 );
        TFloat64Range::Intersect( tplLambdaRange, spcLambdaRange_restframe, intersectedLambdaRange );

        CSpectrum & itplTplWithoutContinuumSpectrum = m_templateWOContRebined_bf;
        CSpectrum & itplTplSpectrum = m_templateRebined_bf;
        CMask & itplMask = m_mskRebined_bf;

        tplWithoutCont.Rebin( intersectedLambdaRange, m_spcSpectralAxis_restframe, itplTplWithoutContinuumSpectrum, itplMask );
        
        const CSpectrumFluxAxis & itplTplWithoutContFluxAxis = itplTplWithoutContinuumSpectrum.GetFluxAxis();

        // same function for the spc with continuum, todo check ?
        tpl.Rebin( intersectedLambdaRange, m_spcSpectralAxis_restframe, itplTplSpectrum, itplMask );

        const CSpectrumFluxAxis & itplTplFluxAxis =  itplTplSpectrum.GetFluxAxis();
        const CSpectrumSpectralAxis& itplTplSpectralAxis = itplTplSpectrum.GetSpectralAxis();

        CMask mask;
        m_spcSpectralAxis_restframe.GetMask( lambdaRange_restframe, mask );
        itplMask &= mask;
        result_corr->Overlap[i] = mask.CompouteOverlapRate( itplMask );
        result_chi->Overlap[i] = result_corr->Overlap[i];

        if( result_corr->Overlap[i] < overlapThreshold )
        {
            result_corr->Status[i] = COperator::nStatus_NoOverlap;
            continue;
        }

        const TFloat64List& error = spcFluxAxis.GetError();

        DebugAssert( !error.empty() );

        Float64 spcMean = 0.0;
        Float64 spcSDev = 0.0;
        if( !spcWithoutContFluxAxis.ComputeMeanAndSDev( itplMask, spcMean, spcSDev, error ) )
        {
            result_corr->Status[i] = COperator::nStatus_DataError;
            continue;
        }

        Float64 tplMean = 0.0;
        Float64 tplSDev = 0.0;
        if( !itplTplWithoutContFluxAxis.ComputeMeanAndSDev( itplMask, tplMean, tplSDev, TFloat64List()) )
        {
            result_corr->Status[i] = COperator::nStatus_DataError;
            continue;
        }

        const TAxisSampleList & Xtpl = itplTplSpectralAxis.GetSamplesVector();
        const TAxisSampleList & YtplWithoutCont = itplTplWithoutContFluxAxis.GetSamplesVector();
        const TAxisSampleList & Ytpl = itplTplFluxAxis.GetSamplesVector();
        const TAxisSampleList & Xspc = m_spcSpectralAxis_restframe.GetSamplesVector();
        const TAxisSampleList & YspcWithoutCont = spcWithoutContFluxAxis.GetSamplesVector();
        const TAxisSampleList & Yspc = spcFluxAxis.GetSamplesVector();

        TFloat64Range logIntersectedLambdaRange( log( intersectedLambdaRange.GetBegin() ), log( intersectedLambdaRange.GetEnd() ) );

        // j cursor move over spectrum
        Int32 j = 0;
        while( j < m_spcSpectralAxis_restframe.GetSamplesCount() && Xspc[j] < logIntersectedLambdaRange.GetBegin() )
            j++;

        // k cursor move over template
        Int32 k = 0;
        while( k < itplTplSpectralAxis.GetSamplesCount() && Xtpl[k] < logIntersectedLambdaRange.GetBegin() )
            k++;

        Int32 jStart = j;

        // xcorr vars
        Float64 sumCorr=0;
        Float64 sumWeight=0;
        // fit vars
        Float64 sumXDevs = 0.0;
        Float64 sumYDevs = 0.0;
        Float64 err2 = 0.0;
        Float64 fit = 0;
        Int32 numDevs = 0;

        // Loop 1, compute the Chisquare coeff.
        // For each sample in the valid lambda range interval.
        while( j<m_spcSpectralAxis_restframe.GetSamplesCount() && Xspc[j] <= logIntersectedLambdaRange.GetEnd() )
        {
            numDevs++;
            err2 = 1.0 / (error[j] * error[j]);
            sumYDevs+=Yspc[j]*err2;
            sumXDevs+=Ytpl[j]*err2;

            j++;
        }

        if ( numDevs==0 || sumYDevs==0 || sumXDevs==0 )
        {
            return -1;
        }

        Float64 ampl = sumYDevs / sumXDevs;
        fit=0;
        Float64 s = 0;

        // Loop 2, compute the Chisquare sum and the Correlation sum.
        j=jStart;
        while( j<m_spcSpectralAxis_restframe.GetSamplesCount() && Xspc[j] <= logIntersectedLambdaRange.GetEnd() )
        {
            int k=j;
            {
                // fit
                fit += pow( Yspc[j] - ampl * Ytpl[k] , 2.0 ) / pow( error[j], 2.0 );
                s += Yspc[j];

                // xcorr
                sumWeight += 1.0 / error[j];
                sumCorr += ( YtplWithoutCont[k] - tplMean ) * ( YspcWithoutCont[j] - spcMean ) / error[j];
            }
            j++;
        }

        // Correlation output
        if( sumWeight>0 )
        {
            result_corr->Correlation[i] = sumCorr / ( tplSDev * spcSDev * sumWeight );
            result_corr->Status[i] = COperator::nStatus_OK;
        }
        // Chisquare output
        fit /= numDevs;
        result_chi->ChiSquare[i] = fit;
        result_chi->Status[i] = COperator::nStatus_OK;
    }

    boost::posix_time::time_duration diff = boost::posix_time::microsec_clock::local_time() - startTime;
    m_TotalDuration = diff.total_seconds();

    return -1;
}



