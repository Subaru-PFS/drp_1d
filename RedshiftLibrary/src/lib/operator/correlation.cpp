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
                                                                   Int32 opt_extinction_unused,
                                                                   Int32 opt_dustFitting_unused,
                                                                   CPriorHelper::TPriorZEList logpriorze_unused, 
                                                                   Bool keepigmism,
                                                                   Float64 FitDustCoeff,
                                                                   Float64 FitMeiksinIdx)
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

    // Pre-Allocate the rebined template and mask with regard to the spectrum size
    m_templateRebined_bf.GetSpectralAxis().SetSize(spectrum.GetSampleCount());
    m_templateRebined_bf.GetFluxAxis().SetSize(spectrum.GetSampleCount());
    m_mskRebined_bf.SetSize(spectrum.GetSampleCount());
    m_spcSpectralAxis_restframe.SetSize(spectrum.GetSampleCount());

    for ( Int32 i=0; i<redshifts.size(); i++)
    {
        result->Correlation[i] = NAN;
        result->Status[i] = nStatus_DataError;
        result->Overlap[i] = 0;

        Float64 onePlusRedshift = 1.0 + redshifts[i];

        //shift lambdaRange backward to be in restframe
        TFloat64Range lambdaRange_restframe( lambdaRange.GetBegin() / onePlusRedshift,
                                             lambdaRange.GetEnd() / onePlusRedshift );

        //redshift in restframe the tgtSpectralAxis
        m_spcSpectralAxis_restframe.ShiftByWaveLength(spcSpectralAxis, onePlusRedshift, CSpectrumSpectralAxis::nShiftBackward);
        TFloat64Range  spcLambdaRange_restframe(spcLambdaRange.GetBegin()/onePlusRedshift,
                                                spcLambdaRange.GetEnd()/onePlusRedshift);

        // Compute clamped lambda range over template in restframe
        TFloat64Range tplLambdaRange;
        retVal = tplSpectralAxis.ClampLambdaRange( lambdaRange_restframe, tplLambdaRange );
        DebugAssert( retVal );

        // if there is any intersection between the lambda range of the spectrum and the lambda range of the template
        // Compute the intersected range
        TFloat64Range intersectedLambdaRange( 0.0, 0.0 );
        TFloat64Range::Intersect( tplLambdaRange, spcLambdaRange_restframe, intersectedLambdaRange );

        CSpectrum& itplSpectrum = m_templateRebined_bf;
        CMask& itplMask = m_mskRebined_bf;

        tpl.Rebin( intersectedLambdaRange, m_spcSpectralAxis_restframe, itplSpectrum, itplMask );

        const CSpectrumFluxAxis& itplTplFluxAxis = itplSpectrum.GetFluxAxis();
        const CSpectrumSpectralAxis& itplTplSpectralAxis = itplSpectrum.GetSpectralAxis();

        CMask mask;
        m_spcSpectralAxis_restframe.GetMask( lambdaRange_restframe, mask );
        itplMask &= mask;

        result->Overlap[i] = mask.CompouteOverlapRate( itplMask );
        if( result->Overlap[i] < overlapThreshold )
        {
            result->Status[i] = nStatus_NoOverlap;
            continue;
        }

        const TFloat64List& error = spcFluxAxis.GetError();

        DebugAssert( !error.empty() );

        Float64 spcMean = 0.0;
        Float64 spcSDev = 0.0;

        if( !spcFluxAxis.ComputeMeanAndSDev( itplMask, spcMean, spcSDev, error ) )
        {
            result->Status[i] = nStatus_DataError;
            continue;
        }

        Float64 tplMean = 0.0;
        Float64 tplSDev = 0.0;
        if( !itplTplFluxAxis.ComputeMeanAndSDev( itplMask, tplMean, tplSDev, TFloat64List() ) )
        {
            result->Status[i] = nStatus_DataError;
            continue;
        }

        const TAxisSampleList & Xtpl = itplTplSpectralAxis.GetSamplesVector();
        const TAxisSampleList & Ytpl = itplTplFluxAxis.GetSamplesVector();
        const TAxisSampleList & Xspc = m_spcSpectralAxis_restframe.GetSamplesVector();
        const TAxisSampleList & Yspc = spcFluxAxis.GetSamplesVector();

        TFloat64Range logIntersectedLambdaRange( log( intersectedLambdaRange.GetBegin() ), log( intersectedLambdaRange.GetEnd() ) );

        // j cursor move over spectrum
        Int32 j = 0;
        while( j < m_spcSpectralAxis_restframe.GetSamplesCount() && Xspc[j] < logIntersectedLambdaRange.GetBegin() )
            j++;

        // j cursor move over template
        Int32 k = 0;
        while( k < itplTplSpectralAxis.GetSamplesCount() && Xtpl[k] < logIntersectedLambdaRange.GetBegin() )
            k++;


        Float64 sumCorr=0;
        Float64 sumWeight=0;

        // k index: move over template
        // j index: move over spectrum

        // Cross correlation function:
        // (t(x) - moy(t)) * (s(x) - moy(s))
        // ---------------------------------
        //          stddev(x) * stddev(s)
        // For each sample in the valid lambda range interval.
        while( j<m_spcSpectralAxis_restframe.GetSamplesCount() && Xspc[j] <= logIntersectedLambdaRange.GetEnd() )
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
