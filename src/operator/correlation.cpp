#include <epic/redshift/operator/correlation.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/fluxaxis.h>
#include <epic/redshift/spectrum/spectralaxis.h>
#include <epic/redshift/continuum/median.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/tools.h>
#include <epic/redshift/common/mask.h>

#include <math.h>
#include <assert.h>

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


const TFloat64List& COperatorCorrelation::GetResults() const
{
    return m_Correlation;
}

const COperatorCorrelation::TStatusList& COperatorCorrelation::GetStatus() const
{
    return m_Status;
}

Float64 COperatorCorrelation::GetComputationDuration() const
{
    return m_TotalDuration;
}

/**
 * Compute correlation factor between spectrum and tpl for each value specified in redshifts.
 * \note Spectrum AND template MUST be in log scale
 */
Bool COperatorCorrelation::Compute(   const CSpectrum& spectrum, const CTemplate& tpl,
                                      const TFloat64Range& lambdaRange, const CRedshifts& redshifts,
                                      Float64 overlapThreshold )
{
    Int32 i,j;
    Bool retVal;

    // Input spectrum MUST be in log scales
    if( spectrum.GetSpectralAxis().IsInLogScale() == false || tpl.GetSpectralAxis().IsInLogScale() == false )
        return false;

    boost::posix_time::ptime  startTime = boost::posix_time::microsec_clock::local_time();

    DebugAssert( overlapThreshold > 0.0 && overlapThreshold <= 1.0 );

    CContinuumMedian continuum;

    CMask spcMask( spectrum.GetSampleCount() );
    CMask itplTplMask( spectrum.GetSampleCount() );
    CSpectrumSpectralAxis shiftedTplSpectralAxis( tpl.GetSampleCount() );
    CSpectrumFluxAxis itplTplFluxAxis( spectrum.GetSampleCount() );

    const CSpectrumSpectralAxis& spcSpectralAxis = spectrum.GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = spectrum.GetFluxAxis();

    const CSpectrumSpectralAxis& tplSpectralAxis = tpl.GetSpectralAxis();
    const CSpectrumFluxAxis& tplFluxAxis = tpl.GetFluxAxis();

    m_Correlation.resize( redshifts.GetRedshiftsCount() );
    m_Overlap.resize( redshifts.GetRedshiftsCount() );
    m_Status.resize( redshifts.GetRedshiftsCount() );

    spectrum.GetSpectralAxis().GetMask( lambdaRange, spcMask );

    UInt32 t1 = spcMask.GetUnMaskedSampleCount();

    // Compute total lambda range, i.e the total lambda range which overlap the spectrum
    TFloat64Range clampedLambdaRange;
    retVal = spectrum.GetSpectralAxis().ClampLambdaRange( lambdaRange, clampedLambdaRange );
    DebugAssert( retVal );


    for ( i=0; i<redshifts.GetRedshiftsCount(); i++)
    {
        m_Correlation[i] = NAN;
        m_Status[i] = nStatus_InvalidInputData;
        m_Overlap[i] = 0;

        Float64 onePlusRedshift = 1.0 + redshifts[i];
        Float64 logOnePlusRedshift = log( onePlusRedshift );

        Int32 shiftedTplStartIndex = tplSpectralAxis.GetIndexAtWaveLength( log( clampedLambdaRange.GetBegin() ) - logOnePlusRedshift );
        Int32 shiftedTplEndIndex = tplSpectralAxis.GetIndexAtWaveLength( log( clampedLambdaRange.GetEnd() ) - logOnePlusRedshift );

        if( shiftedTplStartIndex == -1 || shiftedTplEndIndex == -1 )
        {
            m_Status[i] = nStatus_InvalidInputData;
            continue;
        }

        if ( shiftedTplStartIndex>0 )
            shiftedTplStartIndex--;

        if ( shiftedTplEndIndex < tpl.GetSampleCount() - 1 )
            shiftedTplEndIndex++;

        // Because we are on a log scale basis, shifting goes
        // From: ShiftedLambda = (UnshiftedLambda) * ( 1 + redshift )
        // To:   log( ShiftedLambda ) = log( (UnshiftedLambda) * ( 1 + redshift ) ) = log( UnshiftedLambda  ) + log( 1 + redshift )
        DebugAssert( shiftedTplStartIndex <= shiftedTplEndIndex );
        DebugAssert( shiftedTplEndIndex < shiftedTplSpectralAxis.GetSamplesCount() );
        DebugAssert( shiftedTplEndIndex < tplSpectralAxis.GetSamplesCount() );

        // NOte that everything here happen in log scale
        for ( j=shiftedTplStartIndex; j<=shiftedTplEndIndex; j++ )
        {
            shiftedTplSpectralAxis[j] = tplSpectralAxis[j] + logOnePlusRedshift;
        }

        // template rebin over spectrum lambda
        CSpectrumTools::Interpolate(
                    shiftedTplSpectralAxis,
                    tplFluxAxis,

                    shiftedTplStartIndex,
                    ( shiftedTplEndIndex - shiftedTplStartIndex )+ 1,

                    spcSpectralAxis,
                    itplTplFluxAxis,

                    itplTplMask
                    );

        // Compute the intersection set bewteen interpolated template mask and spectrum mask
        // and stoer it in itplTplMask
        retVal = itplTplMask.IntersectWith( spcMask );
        DebugAssert( retVal );

        // Compute overlap rate between intersection set and template
        m_Overlap[i] = spcMask.CompouteOverlapRate( itplTplMask );
        DebugAssert( m_Overlap[i] != -1.0 );

        if( m_Overlap[i] < overlapThreshold )
        {
            m_Status[i] = nStatus_AboveOverlapThreshold;
            continue;
        }

        const Float64* error = spcFluxAxis.GetError();

        DebugAssert( error!= NULL );

        Float64 spcMean = 0.0;
        Float64 spcSDev = 0.0;
        if( !spcFluxAxis.ComputeMeanAndSDev( spcMask, spcMean, spcSDev, error ) )
        {
            m_Status[i] = nStatus_InvalidInputData;
            continue;
        }

        Float64 tplMean = 0.0;
        Float64 tplSDev = 0.0;
        if( !tplFluxAxis.ComputeMeanAndSDev( spcMask, tplMean, tplSDev, NULL ) )
        {
            m_Status[i] = nStatus_InvalidInputData;
            continue;
        }

        Float64 sumCorr=0;
        Float64 sumWeight=0;
        for (j=0;j<spectrum.GetSampleCount();j++)
        {
            if (itplTplMask[j])
            {
                sumWeight += 1.0 / error[j];
                sumCorr += ( itplTplFluxAxis[j] - tplMean ) * ( spcFluxAxis[j] - spcMean ) / error[j];
            }
        }

        if( sumWeight>0 )
        {
            m_Correlation[i] = sumCorr / ( tplSDev * spcSDev * sumWeight );
            m_Status[i] = nStatus_OK;
        }
    }

    boost::posix_time::time_duration diff = boost::posix_time::microsec_clock::local_time() - startTime;

    m_TotalDuration = diff.total_seconds();

    return true;
}


