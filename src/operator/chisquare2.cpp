#include <epic/redshift/operator/chisquare2.h>

#include <epic/redshift/spectrum/axis.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/tools.h>
#include <epic/redshift/common/mask.h>
#include <epic/redshift/operator/chisquareresult.h>
#include <epic/redshift/extremum/extremum.h>

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


Void COperatorChiSquare2::BasicFit( const CSpectrum& spectrum, const CTemplate& tpl, Float64* pfgTplBuffer,
                                const TFloat64Range& lambdaRange, Float64 redshift, Float64 overlapThreshold,
                                Float64& overlapRate, Float64& chiSquare, EStatus& status , Float64 amplitude)
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

    //CSpectrumFluxAxis::Rebin( intersectedLambdaRange, tplFluxAxis, shiftedTplSpectralAxis, spcSpectralAxis, itplTplFluxAxis, itplTplSpectralAxis, itplMask );
    CSpectrumFluxAxis::Rebin2( intersectedLambdaRange, tplFluxAxis, pfgTplBuffer, redshift, m_shiftedTplSpectralAxis_bf, spcSpectralAxis, itplTplFluxAxis, itplTplSpectralAxis, itplMask );

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

    //EZ formulation
    //Float64 sumXDevs = 0.0;
    //Float64 sumYDevs = 0.0;
    // Tonry&Davis formulation
    Float64 sumCross = 0.0;
    Float64 sumT = 0.0;

    Float64 err2 = 0.0;
    Float64 fit = 0;
    Int32 numDevs = 0;
    const Float64* error = spcFluxAxis.GetError();



    //if(0)
    while( j<spcSpectralAxis.GetSamplesCount() && Xspc[j] <= currentRange.GetEnd() )
    {
        numDevs++;
        err2 = 1.0 / (error[j] * error[j]);
        //EZ formulation
        //sumYDevs+=Yspc[j]*err2;
        //sumXDevs+=Ytpl[j]*err2;
        //sumYDevs+=Yspc[j];
        //sumXDevs+=Ytpl[j];

        // Tonry&Davis formulation
        sumCross+=Yspc[j]*Ytpl[j]*err2;
        sumT+=Ytpl[j]*Ytpl[j]*err2;
        //sumCross+=Yspc[j]*Ytpl[j];
        //sumT+=Ytpl[j]*Ytpl[j];

        j++;
    }

    //if ( numDevs==0 || sumYDevs==0 || sumXDevs==0 ) //EZ formulation
    if ( numDevs==0 || sumCross==0 || sumT==0 ) // Tonry&Davis formulation
    {
        status = nStatus_DataError;
        return;
    }

    //Float64 ampl = 1.0;
    //Float64 ampl = sumYDevs / sumXDevs; //EZ formulation
    //Float64 ampl = sumT; //EZ formulation
    Float64 ampl = max(0.0, sumCross / sumT); // Tonry&Davis formulation
    if(amplitude !=-1){
        ampl = amplitude;
    }

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
            //fit += pow( Yspc[j] - ampl * Ytpl[k] , 2.0 );
            s += Yspc[j];
        }
        j++;
    }

    // Chi square reduct: it can introduces some problem?
    fit /= numDevs;


    chiSquare = fit;
    //chiSquare = ampl;
    status = nStatus_OK;
}


const COperatorResult* COperatorChiSquare2::Compute(const CSpectrum& spectrum, const CTemplate& tpl,
                          const TFloat64Range& lambdaRange, const TFloat64List& redshifts,
                          Float64 overlapThreshold )
{

    if( spectrum.GetSpectralAxis().IsInLinearScale() == false || tpl.GetSpectralAxis().IsInLinearScale() == false )
    {
        Log.LogError("Chisquare2, input spectrum or template are not in log scale (ignored)");
        //return NULL;
    }

    // Pre-Allocate the rebined template and mask with regard to the spectrum size
    m_templateRebined_bf.GetSpectralAxis().SetSize(spectrum.GetSampleCount());
    m_templateRebined_bf.GetFluxAxis().SetSize(spectrum.GetSampleCount());
    m_mskRebined_bf.SetSize(spectrum.GetSampleCount());
    m_shiftedTplSpectralAxis_bf.SetSize( tpl.GetSampleCount());

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

    // pfg with std::vector
    //CTemplate       templateFine;
    //templateFine.GetSpectralAxis().SetSize(nTgt);
    //templateFine.GetFluxAxis().SetSize(nTgt);
    //Float64* precomputedFineGridTplFlux = templateFine.GetFluxAxis().GetSamples();
    // pfg with malloc
    Float64* precomputedFineGridTplFlux = (Float64*)malloc(nTgt*sizeof(Float64));
    // pfg with static array => doesn't work
    //nTgt = 999999;
    //Float64 precomputedFineGridTplFlux[999999];
    //Log.LogInfo( "nTgt: %d samples", nTgt);

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
        if(x < tplSpectralAxis[0] || x > tplSpectralAxis[n-1]){
            precomputedFineGridTplFlux[k] = 0.0;
        }else{
            //precomputedFineGridTplFlux[k] = gsl_interp_eval(interpolation, Xsrc, Ysrc, x, accelerator);
            precomputedFineGridTplFlux[k] = gsl_spline_eval (spline, x, accelerator);
        }
    }
    //*/

    /*//debug:
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

    TFloat64List sortedRedshifts = redshifts;
    std::sort(sortedRedshifts.begin(), sortedRedshifts.end());

    CChisquareResult* result = new CChisquareResult();
    result->ChiSquare.resize( sortedRedshifts.size() );
    result->Redshifts.resize( sortedRedshifts.size() );
    result->Overlap.resize( sortedRedshifts.size() );
    result->Status.resize( sortedRedshifts.size() );

    result->Redshifts = sortedRedshifts;


    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        BasicFit( spectrum, tpl, precomputedFineGridTplFlux, lambdaRange, result->Redshifts[i], overlapThreshold, result->Overlap[i], result->ChiSquare[i], result->Status[i] );
    }

    // extrema
    Int32 extremumCount = 10;
    TPointList extremumList;
    TFloat64Range redshiftsRange(result->Redshifts[0], result->Redshifts[result->Redshifts.size()-1]);
    CExtremum extremum( redshiftsRange, extremumCount, true);
    extremum.Find( result->Redshifts, result->ChiSquare, extremumList );
    // Refine Extremum with a second maximum search around the z candidates:
    // This corresponds to the finer xcorrelation in EZ Pandora (in standard_DP fctn in SolveKernel.py)
    Float64 radius = 0.001;
    for( Int32 i=0; i<extremumList.size(); i++ )
    {
        Float64 x = extremumList[i].X;
        Float64 left_border = max(redshiftsRange.GetBegin(), x-radius);
        Float64 right_border=min(redshiftsRange.GetEnd(), x+radius);

        TPointList extremumListFine;
        TFloat64Range rangeFine = TFloat64Range( left_border, right_border );
        CExtremum extremumFine( rangeFine , 1, true);
        extremumFine.Find( result->Redshifts, result->ChiSquare, extremumListFine );
        extremumList[i] = extremumListFine[0];
    }
    // store extrema results
    result->Extrema.resize( extremumCount );
    for( Int32 i=0; i<extremumList.size(); i++ )
    {

        result->Extrema[i] = extremumList[i].X;
    }

    free(precomputedFineGridTplFlux);
    return result;

}

const COperatorResult* COperatorChiSquare2::ExportChi2versusAZ(const CSpectrum& spectrum, const CTemplate& tpl,
                          const TFloat64Range& lambdaRange, const TFloat64List& redshifts,
                          Float64 overlapThreshold )
{
    if( spectrum.GetSpectralAxis().IsInLinearScale() == false || tpl.GetSpectralAxis().IsInLinearScale() == false )
    {
        Log.LogError("Chisquare2, input spectrum or template are not in log scale (ignored)");
        //return NULL;
    }

    // Pre-Allocate the rebined template and mask with regard to the spectrum size
    m_templateRebined_bf.GetSpectralAxis().SetSize(spectrum.GetSampleCount());
    m_templateRebined_bf.GetFluxAxis().SetSize(spectrum.GetSampleCount());
    m_mskRebined_bf.SetSize(spectrum.GetSampleCount());
    m_shiftedTplSpectralAxis_bf.SetSize( tpl.GetSampleCount());


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

    // pfg with std::vector
    //CTemplate       templateFine;
    //templateFine.GetSpectralAxis().SetSize(nTgt);
    //templateFine.GetFluxAxis().SetSize(nTgt);
    //Float64* precomputedFineGridTplFlux = templateFine.GetFluxAxis().GetSamples();
    // pfg with malloc
    Float64* precomputedFineGridTplFlux = (Float64*)malloc(nTgt*sizeof(Float64));
    // pfg with static array => doesn't work
    //nTgt = 999999;
    //Float64 precomputedFineGridTplFlux[999999];
    //Log.LogInfo( "nTgt: %d samples", nTgt);

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
        if(x < tplSpectralAxis[0] || x > tplSpectralAxis[n-1]){
            precomputedFineGridTplFlux[k] = 0.0;
        }else{
            //precomputedFineGridTplFlux[k] = gsl_interp_eval(interpolation, Xsrc, Ysrc, x, accelerator);
            precomputedFineGridTplFlux[k] = gsl_spline_eval (spline, x, accelerator);
        }
    }
    //*/

    /*//debug:
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

    TFloat64List sortedRedshifts;// = redshifts;
    // for sc_020086471_F02P016_vmM1_red_107_1_atm_clean : zref = 1.3455, aref=4.6772031621836956e-17, zcalc = 2.1634
    Float64 zcenter = 1.3455;
    Float64 acenter = 4.6772031621836956e-17;
    Float64 arange = acenter/2.0/2.0; //
    //Float64 zcenter = 2.1634;
    //Float64 acenter = 3.3228859582551038e-17;
    //Float64 arange = acenter/4.0/2.0; //
    // for
    //Float64 zcenter = 1.134;
    //Float64 acenter = 4.6772031621836956e-17;
    //Float64 arange = acenter/2.0/2.0; //


    // fill the redshift grid
    //sortedRedshifts.push_back(zcenter);
    Float64 zmin = zcenter - 0.025;
    Float64 zmax = zcenter + 0.025;
    Float64 zstep = 0.0001;
    Int32 nz = (Int32)((zmax-zmin)/zstep);
    for (Int32 i=0;i<nz;i++)
    {
        Float64 z = zmin + zstep*i;
        sortedRedshifts.push_back(z);
    }
    std::sort(sortedRedshifts.begin(), sortedRedshifts.end());

    // fill the amplitude grid
    TFloat64List sortedAmplitudes;
    //sortedAmplitudes.push_back(-1); //auto amplitude
    Float64 amin = acenter - arange;
    Float64 amax = acenter + arange;
    Int32 na = 200;
    Float64 astep = ((amax-amin)/na);
    for (Int32 i=0;i<na;i++)
    {
        Float64 a = amin + astep*i;
        sortedAmplitudes.push_back(a);
    }
    std::sort(sortedAmplitudes.begin(), sortedAmplitudes.end());

    CChisquareResult* result = new CChisquareResult();
    result->ChiSquare.resize( sortedRedshifts.size() );
    result->Redshifts.resize( sortedRedshifts.size() );
    result->Overlap.resize( sortedRedshifts.size() );
    result->Status.resize( sortedRedshifts.size() );

    result->Redshifts = sortedRedshifts;


    FILE* fa = fopen( "chi2_versus_ampl_z_axes_ampl.txt", "w+" );
    for (Int32 j=0;j<sortedAmplitudes.size();j++)
    {
        Float64 ampl = sortedAmplitudes[j];
        fprintf( fa, "%.15e\n", ampl);
    }
    fclose( fa );

    FILE* fz = fopen( "chi2_versus_ampl_z_axes_z.txt", "w+" );
    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        fprintf( fz, "%.15e\n", sortedRedshifts[i]);
    }
    fclose( fz );

    FILE* f = fopen( "chi2_versus_ampl_z.txt", "w+" );
    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        for (Int32 j=0;j<sortedAmplitudes.size();j++)
        {
            Float64 ampl = sortedAmplitudes[j];
            BasicFit( spectrum, tpl, precomputedFineGridTplFlux, lambdaRange, result->Redshifts[i], overlapThreshold, result->Overlap[i], result->ChiSquare[i], result->Status[i], ampl );
            fprintf( f, "%.15e", result->ChiSquare[i]);
            if(j<sortedAmplitudes.size()-1){
                fprintf( f, "\t");
            }
        }
        fprintf( f, "\n");
    }
    fclose( f );



    free(precomputedFineGridTplFlux);
    return result;

}