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


void COperatorChiSquare::BasicFit( const CSpectrum& spectrum, const CTemplate& tpl,
                                const TFloat64Range& lambdaRange, Float64 redshift, Float64 overlapThreshold,
                                Float64& overlapRate, Float64& chiSquare, Float64& fitamplitude, EStatus& status )
{
    chiSquare = boost::numeric::bounds<float>::highest();
    fitamplitude = -1.0;
    overlapRate = 0.0;
    status = nStatus_DataError;

    const CSpectrumSpectralAxis& spcSpectralAxis = spectrum.GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = spectrum.GetFluxAxis();

    const CSpectrumSpectralAxis& tplSpectralAxis = tpl.GetSpectralAxis();
    const CSpectrumFluxAxis& tplFluxAxis = tpl.GetFluxAxis();


    Float64 onePlusRedshift = 1.0 + redshift;

    //shift lambdaRange backward to be in restframe
    TFloat64Range spcLambdaRange_restframe;
    TFloat64Range lambdaRange_restframe( lambdaRange.GetBegin() / onePlusRedshift,
                                         lambdaRange.GetEnd() / onePlusRedshift );

    //redshift in restframe the spectrum spectralaxis
    m_spcSpectralAxis_restframe.ShiftByWaveLength(spcSpectralAxis, onePlusRedshift, CSpectrumSpectralAxis::nShiftBackward);
    m_spcSpectralAxis_restframe.ClampLambdaRange( lambdaRange_restframe, spcLambdaRange_restframe );

    // Compute clamped lambda range over template
    TFloat64Range tplLambdaRange;
    tplSpectralAxis.ClampLambdaRange( lambdaRange_restframe, tplLambdaRange );
    // Compute the intersected range
    TFloat64Range intersectedLambdaRange( 0.0, 0.0 );
    TFloat64Range::Intersect( tplLambdaRange, spcLambdaRange_restframe, intersectedLambdaRange );

    CSpectrum& itplTplSpectrum = m_templateRebined_bf;
    CMask&  itplMask = m_mskRebined_bf;

    tpl.Rebin( intersectedLambdaRange, m_spcSpectralAxis_restframe, itplTplSpectrum, itplMask );

    CSpectrumFluxAxis& itplTplFluxAxis = itplTplSpectrum.GetFluxAxis();
    const CSpectrumSpectralAxis& itplTplSpectralAxis = itplTplSpectrum.GetSpectralAxis();

    CMask mask;
    m_spcSpectralAxis_restframe.GetMask( lambdaRange_restframe, mask );
    itplMask &= mask;
    overlapRate = mask.CompouteOverlapRate( itplMask );
    // Check for overlap rate
    if( overlapRate < overlapThreshold )
    {
        status = nStatus_NoOverlap;
        return ;
    }

    const TAxisSampleList & Xtpl = itplTplSpectralAxis.GetSamplesVector();
    const TAxisSampleList & Ytpl = itplTplFluxAxis.GetSamplesVector();
    const TAxisSampleList & Xspc = m_spcSpectralAxis_restframe.GetSamplesVector();
    const TAxisSampleList & Yspc = spcFluxAxis.GetSamplesVector();
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
    while( j < m_spcSpectralAxis_restframe.GetSamplesCount() && Xspc[j] < currentRange.GetBegin() )
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
    const TFloat64List& error = spcFluxAxis.GetError();

    while( j<m_spcSpectralAxis_restframe.GetSamplesCount() && Xspc[j] <= currentRange.GetEnd() )
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

    while( j<m_spcSpectralAxis_restframe.GetSamplesCount() && Xspc[j] <= currentRange.GetEnd() )
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
                          Float64 overlapThreshold , std::vector<CMask> additional_spcMasks_unused, string opt_interp_unused,
                                                               Int32 opt_extinction_unused,
                                                               Int32 opt_dustFitting_unused,
                                                               CPriorHelper::TPriorZEList logpriorze_unused,
                                                               Bool keepigmism,
                                                               Float64 FitDustCoeff,
                                                               Float64 FitMeiksinIdx)
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

    // Pre-Allocate the rebined template and mask with regard to the spectrum size
    m_templateRebined_bf.GetSpectralAxis().SetSize(spectrum.GetSampleCount());
    m_templateRebined_bf.GetFluxAxis().SetSize(spectrum.GetSampleCount());
    m_mskRebined_bf.SetSize(spectrum.GetSampleCount());
    m_spcSpectralAxis_restframe.SetSize(spectrum.GetSampleCount());


    for (Int32 i=0;i<redshifts.size();i++)
    {
        BasicFit( spectrum, tpl, lambdaRange, result->Redshifts[i], overlapThreshold, result->Overlap[i], result->ChiSquare[i], result->FitAmplitude[i], result->Status[i] );
    }


    return result;

}
