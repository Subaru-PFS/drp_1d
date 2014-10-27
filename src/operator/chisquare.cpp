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
                                const TFloat64Range& lambdaRange, Float64 redshift, Float64 overlapThreshold,
                                Float64& ampl, Float64& overlapRate, Int32& numDevs )
{

    if( spectrum.GetSpectralAxis().IsInLinearScale() == false || tpl.GetSpectralAxis().IsInLinearScale() == false )
    {
        Log.LogError("Failed to compute Chi Square, input spetrum or template are not in linear scale");
        return false;
    }

    CMask mask( spectrum.GetSampleCount() );
    CMask itplMask( spectrum.GetSampleCount() );
    CSpectrumSpectralAxis shiftedTplSpectralAxis( tpl.GetSampleCount() );
    CSpectrumFluxAxis itplTplFluxAxis( spectrum.GetSampleCount() );

    const CSpectrumSpectralAxis& spcSpectralAxis = spectrum.GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = spectrum.GetFluxAxis();

    const CSpectrumAxis& tplSpectralAxis = tpl.GetSpectralAxis();
    const CSpectrumAxis& tplFluxAxis = tpl.GetFluxAxis();

    spectrum.GetSpectralAxis().GetMask( lambdaRange, mask );

    Float64 onePlusRedshift = 1.0 + redshift;
    Float64 logOnePlusRedshift = log( onePlusRedshift );

    // Compute shifted template
    for( Int32 i = 0; i<tplSpectralAxis.GetSamplesCount(); i++ )
    {
        shiftedTplSpectralAxis[i] = tplSpectralAxis[i] * onePlusRedshift;
    }

    // template rebin over spectrum lambda
    CSpectrumTools::Interpolate(
                shiftedTplSpectralAxis,
                tplFluxAxis,
                0,
                tpl.GetSampleCount(),

                spcSpectralAxis,
                itplTplFluxAxis,
                itplMask
                );

    if( itplMask.IntersectWith( mask ) == false )
        return false;

    overlapRate = mask.CompouteOverlapRate( itplMask );
    if( overlapRate < overlapThreshold )
        return false;

    Float64 sum_xdevs = 0.0;
    Float64 sum_ydevs = 0.0;
    Float64 err2 = 0.0;
    Float64 fit = 0;

    const Float64* spcFluxError = spcFluxAxis.GetError();
    for( Int32 i = 0; i<spectrum.GetSampleCount(); i++ )
    {
        if( itplMask[i] )
        {
            numDevs++;
            err2 = 1.0 / pow( spcFluxError[i], 2.0 );
            sum_ydevs+=spcFluxAxis[i]*err2;
            sum_xdevs+=itplTplFluxAxis[i]*err2;
        }
    }

    if( numDevs > 0 )
    {
        ampl = sum_ydevs / sum_xdevs;
    }
    else
    {
        return false;
    }

    if ( sum_ydevs==0 )
    {
        return false;
    }

    if ( sum_xdevs==0 )
    {
        return false;
    }

    fit=0;
    for ( Int32 j=0;j<spectrum.GetSampleCount();j++)
    {
        if (itplMask[j])
        {
            fit += pow( spcFluxAxis[j] - ampl * itplTplFluxAxis[j] , 2.0 ) / pow( spcFluxError[j], 2.0 );
        }
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
        Float64 ampl = 0.0;
        Float64 overlapRate = 0.0;
        Int32   numDevs = 0;

        m_Chisquare[i] = BasicFit( spectrum, tpl, lambdaRange, redshifts[i], overlapThreshold, ampl, m_Overlap[i], numDevs );

    }


    return true;
}
