#include <epic/redshift/operator/chisquare2.h>

#include <epic/redshift/spectrum/axis.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/tools.h>
#include <epic/redshift/common/mask.h>
#include <epic/redshift/operator/chisquareresult.h>

#include <epic/core/log/log.h>

#include <boost/numeric/conversion/bounds.hpp>

#include <math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <algorithm>    // std::sort

#include <assert.h>

#define NOT_OVERLAP_VALUE NAN
#include <stdio.h>

using namespace NSEpic;
using namespace std;

COperatorChiSquare2::COperatorChiSquare2()
{

}

COperatorChiSquare2::~COperatorChiSquare2()
{

}


Void COperatorChiSquare2::BasicFit( const CSpectrum& spectrum, const CTemplate& tpl, CTemplate& tplFineBuffer,
                                const TFloat64Range& lambdaRange, Float64 redshift, Float64 overlapThreshold,
                                Float64& overlapRate, Float64& chiSquare, EStatus& status )
{
    chiSquare = boost::numeric::bounds<float>::highest();
    overlapRate = 0.0;
    status = nStatus_DataError;

    Bool retVal;

    //CSpectrumSpectralAxis shiftedTplSpectralAxis( tpl.GetSampleCount(), false );

    const CSpectrumSpectralAxis& spcSpectralAxis = spectrum.GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = spectrum.GetFluxAxis();

    const CSpectrumSpectralAxis& tplSpectralAxis = tpl.GetSpectralAxis();
    const CSpectrumFluxAxis& tplFluxAxis = tpl.GetFluxAxis();

    // Compute clamped lambda range over spectrum
    TFloat64Range spcLambdaRange;
    retVal = spcSpectralAxis.ClampLambdaRange( lambdaRange, spcLambdaRange );

    // Compute shifted template
    Float64 onePlusRedshift = 1.0 + redshift;
    m_shiftedTplSpectralAxis_bf.ShiftByWaveLength( tplSpectralAxis, onePlusRedshift, CSpectrumSpectralAxis::nShiftForward );
    TFloat64Range intersectedLambdaRange( 0.0, 0.0 );

    // Compute clamped lambda range over template
    TFloat64Range tplLambdaRange;
    retVal = m_shiftedTplSpectralAxis_bf.ClampLambdaRange( lambdaRange, tplLambdaRange );

    // if there is any intersection between the lambda range of the spectrum and the lambda range of the template
    // Compute the intersected range
    TFloat64Range::Intersect( tplLambdaRange, spcLambdaRange, intersectedLambdaRange );

    //UInt32 tgtn = spcSpectralAxis.GetSamplesCount() ;
    CSpectrumFluxAxis& itplTplFluxAxis = m_templateRebined_bf.GetFluxAxis();
    CSpectrumSpectralAxis& itplTplSpectralAxis = m_templateRebined_bf.GetSpectralAxis();
    CMask& itplMask = m_mskRebined_bf;

    CSpectrumFluxAxis& tplFineFluxAxis = tplFineBuffer.GetFluxAxis();

    //CSpectrumFluxAxis::Rebin( intersectedLambdaRange, tplFluxAxis, shiftedTplSpectralAxis, spcSpectralAxis, itplTplFluxAxis, itplTplSpectralAxis, itplMask );
    CSpectrumFluxAxis::Rebin2( intersectedLambdaRange, tplFluxAxis, tplFineFluxAxis, redshift, m_shiftedTplSpectralAxis_bf, spcSpectralAxis, itplTplFluxAxis, itplTplSpectralAxis, itplMask );

    /*//overlapRate, Method 1
    CMask mask;
    spcSpectralAxis.GetMask( lambdaRange, mask );
    itplMask &= mask;
    overlapRate = mask.CompouteOverlapRate( itplMask );
    //*/

    //overlapRate, Method 2
    //CMask mask;
    //spcSpectralAxis.GetMask( lambdaRange, mask );
    //overlapRate = mask.IntersectAndComputeOverlapRate( itplMask );

    //overlapRate, Method 3
    overlapRate = spcSpectralAxis.IntersectMaskAndComputeOverlapRate( lambdaRange, itplMask );

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

    // Tonry&Davis formulation
    Float64 sumCross = 0.0;
    Float64 sumT = 0.0;

    //if(0)
    while( j<spcSpectralAxis.GetSamplesCount() && Xspc[j] <= currentRange.GetEnd() )
    {
        numDevs++;
        err2 = 1.0 / (error[j] * error[j]);
        sumYDevs+=Yspc[j]*err2;
        sumXDevs+=Ytpl[j]*err2;

        // Tonry&Davis formulation
        //sumCross+=Yspc[j]*Ytpl[j]*err2;
        //sumT+=Ytpl[j]*Ytpl[j]*err2;

        j++;
    }

    if ( numDevs==0 || sumYDevs==0 || sumXDevs==0 )
    {
        status = nStatus_DataError;
        return;
    }

    //Float64 ampl = 1.0;
    Float64 ampl = sumYDevs / sumXDevs; //EZ formulation
    //Float64 ampl = sumCross / sumT; // Tonry&Davis formulation

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
    status = nStatus_OK;
}


const COperatorResult* COperatorChiSquare2::Compute(const CSpectrum& spectrum, const CTemplate& tpl,
                          const TFloat64Range& lambdaRange, const TFloat64List& redshifts,
                          Float64 overlapThreshold )
{

    if( spectrum.GetSpectralAxis().IsInLinearScale() == false || tpl.GetSpectralAxis().IsInLinearScale() == false )
    {
        Log.LogError("Failed to compute Cross correlation, input spectrum or template are not in log scale");
        return NULL;
    }

    // Pre-Allocate the rebined template and mask with regard to the spectrum size
    m_templateRebined_bf.GetSpectralAxis().SetSize(spectrum.GetSampleCount());
    m_templateRebined_bf.GetFluxAxis().SetSize(spectrum.GetSampleCount());
    m_mskRebined_bf.SetSize(spectrum.GetSampleCount());
    m_shiftedTplSpectralAxis_bf.SetSize( tpl.GetSampleCount());


    TFloat64List sortedRedshifts = redshifts;
    std::sort(sortedRedshifts.begin(), sortedRedshifts.end());
    //*/
    // Precalculate a fine grid template to be used for the 'closest value' rebin method
    Int32 n = tpl.GetSampleCount();
    CSpectrumFluxAxis tplFluxAxis = tpl.GetFluxAxis();
    CSpectrumSpectralAxis tplSpectralAxis = tpl.GetSpectralAxis();
    //Float64 dLambdaTgt =  1.0 * ( spectrum.GetMeanResolution()*0.9 )/( 1+sortedRedshifts[sortedRedshifts.size()-1] );
    Float64 dLambdaTgt =  0.1;
    //Float64 lmin = tplSpectralAxis[0];
    Float64 lmin = 0;
    Float64 lmax = tplSpectralAxis[n-1];
    Int32 nTgt = (lmax-lmin)/dLambdaTgt + 2.0/dLambdaTgt;
    CTemplate       templateFine;
    templateFine.GetSpectralAxis().SetSize(nTgt);
    templateFine.GetFluxAxis().SetSize(nTgt);
    Float64* Yfine = templateFine.GetFluxAxis().GetSamples();
    Float64* Xfine = templateFine.GetSpectralAxis().GetSamples();

    //inialise and allocate the gsl objects
    Float64* Ysrc = tplFluxAxis.GetSamples();
    Float64* Xsrc = tplSpectralAxis.GetSamples();
    // linear
    //gsl_interp *interpolation = gsl_interp_alloc (gsl_interp_linear,n);
    //gsl_interp_init(interpolation, Xsrc, Ysrc, n);
    //gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();

    //spline
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
    gsl_spline_init (spline, Xsrc, Ysrc, n);
    gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();

    Int32 k = 0;
    Float64 x = 0.0;
    for(k=0; k<nTgt; k++){
        x = lmin + k*dLambdaTgt;
        Xfine[k] = x;
        if(x < tplSpectralAxis[0] || x > tplSpectralAxis[n-1]){
            Yfine[k] = 0.0;
        }else{
            //Yfine[k] = gsl_interp_eval(interpolation, Xsrc, Ysrc, x, accelerator);
            Yfine[k] = gsl_spline_eval (spline, x, accelerator);
        }
    }
    //*/

    //*//debug:
    // save templateFine
    FILE* f = fopen( "template_fine.txt", "w+" );
    for(Int32 m=0; m<nTgt; m++){
        if( Yfine[m] < 0.0001 ){
            fprintf( f, "%e %e\n", Xfine[m], Yfine[m]);
        }else{
            fprintf( f, "%f %f\n", Xfine[m], Yfine[m]);
        }
    }
    fclose( f );
    //*/

    // use original tpl
    //CTemplate _tpl = tpl;
    // use fine tpl
    //CTemplate _tpl = templateFine;

    CChisquareResult* result = new CChisquareResult();
    result->ChiSquare.resize( sortedRedshifts.size() );
    result->Redshifts.resize( sortedRedshifts.size() );
    result->Overlap.resize( sortedRedshifts.size() );
    result->Status.resize( sortedRedshifts.size() );

    result->Redshifts = sortedRedshifts;


    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        BasicFit( spectrum, tpl, templateFine, lambdaRange, result->Redshifts[i], overlapThreshold, result->Overlap[i], result->ChiSquare[i], result->Status[i] );
    }


    return result;

}
