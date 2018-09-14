#include <RedshiftLibrary/spectrum/fluxaxis.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/common/median.h>
#include <RedshiftLibrary/common/mean.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/spectralaxis.h>

#include <math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <RedshiftLibrary/log/log.h>

using namespace NSEpic;
using namespace std;

CSpectrumFluxAxis::CSpectrumFluxAxis()
{

}

CSpectrumFluxAxis::CSpectrumFluxAxis( UInt32 n ) :
    CSpectrumAxis( n ),
    m_StatError( n )
{
    for( UInt32 i=0; i<n; i++ )
    {
        m_StatError[i] = 1.0;
    }
}

CSpectrumFluxAxis::CSpectrumFluxAxis( const Float64* samples, UInt32 n ) :
    CSpectrumAxis( samples, n ),
    m_StatError( n )
{
    for( UInt32 i=0; i<n; i++ )
    {
        m_StatError[i] = 1.0;
    }
}

CSpectrumFluxAxis::CSpectrumFluxAxis( const Float64* samples, UInt32 n,
				      const Float64* error, UInt32 m ) :
    CSpectrumAxis( samples, n ),
    m_StatError( n )
{
    for( UInt32 i=0; i<n; i++ )
    {
        m_StatError[i] = error[i];
    }
}

CSpectrumFluxAxis::~CSpectrumFluxAxis()
{

}

CSpectrumFluxAxis& CSpectrumFluxAxis::operator=(const CSpectrumFluxAxis& other)
{
    m_StatError = other.m_StatError;
    CSpectrumAxis::operator=( other );
    return *this;
}

Bool CSpectrumFluxAxis::Rebin( const TFloat64Range& range, const CSpectrumFluxAxis& sourceFluxAxis, const CSpectrumSpectralAxis& sourceSpectralAxis, const CSpectrumSpectralAxis& targetSpectralAxis,
                               CSpectrumFluxAxis& rebinedFluxAxis, CSpectrumSpectralAxis& rebinedSpectralAxis, CMask& rebinedMask  )
{
    if( sourceFluxAxis.GetSamplesCount() != sourceSpectralAxis.GetSamplesCount() )
        return false;

    //the spectral axis should be in the same scale
    TFloat64Range logIntersectedLambdaRange( log( range.GetBegin() ), log( range.GetEnd() ) );
    TFloat64Range currentRange = logIntersectedLambdaRange;
    if( sourceSpectralAxis.IsInLinearScale() != targetSpectralAxis.IsInLinearScale() )
        return false;
    if(sourceSpectralAxis.IsInLinearScale()){
        currentRange = range;
    }

    rebinedFluxAxis.SetSize( targetSpectralAxis.GetSamplesCount() );
    rebinedSpectralAxis.SetSize( targetSpectralAxis.GetSamplesCount() );
    rebinedMask.SetSize( targetSpectralAxis.GetSamplesCount() );

    const Float64* Xsrc = sourceSpectralAxis.GetSamples();
    const Float64* Ysrc = sourceFluxAxis.GetSamples();
    const Float64* Xtgt = targetSpectralAxis.GetSamples();
    Float64* Yrebin = rebinedFluxAxis.GetSamples();
    Float64* Xrebin = rebinedSpectralAxis.GetSamples();

    // Move cursors up to lambda range start
    Int32 j = 0;
    while( j<targetSpectralAxis.GetSamplesCount() && Xtgt[j] < currentRange.GetBegin() )
    {
        rebinedMask[j] = 0;
        Yrebin[j] = 0.0;
        j++;
    }
    // Include ALL lambda range
    //if( j > 0 )
    //    j--;

    Int32 k = 0;

    // For each sample in the valid lambda range interval.
    while( k<sourceSpectralAxis.GetSamplesCount()-1 && Xsrc[k] <= currentRange.GetEnd() )
    {
        // For each sample in the target spectrum that are in between two continous source sample
        while( j<targetSpectralAxis.GetSamplesCount() && Xtgt[j] <= Xsrc[k+1] )
        {
            // perform linear interpolation of the flux
            Float64 t = ( Xtgt[j] - Xsrc[k] ) / ( Xsrc[k+1] - Xsrc[k] );
            Xrebin[j] = Xsrc[k] + ( Xsrc[k+1] - Xsrc[k] ) * t;
            Yrebin[j] = Ysrc[k] + ( Ysrc[k+1] - Ysrc[k] ) * t;

            rebinedMask[j] = 1;

            j++;
        }

        k++;
    }

    while( j < targetSpectralAxis.GetSamplesCount() )
    {
        rebinedMask[j] = 0;
        Yrebin[j] = 0.0;
        j++;
    }

    return true;
}

///
/// * This rebin method targets processing speed:
/// - it uses already allocated rebinedFluxAxis, rebinedSpectralAxis and rebinedMask
/// - opt_interp = 'lin' : linear interpolation is performed
/// - opt_interp = 'precomputedfinegrid' : nearest grid point interpolation is performed using pfgTplBuffer which is the precomputed fine grid
/// - opt_interp = 'spline' : GSL/spline interpolation is performed (TODO - not tested)
/// - opt_interp = 'ngp' : nearest grid point is performed (TODO - not tested)
///
Bool CSpectrumFluxAxis::Rebin2( const TFloat64Range& range, const CSpectrumFluxAxis& sourceFluxAxis, const Float64* pfgTplBuffer, Float64 sourcez, const CSpectrumSpectralAxis& sourceSpectralAxis, const CSpectrumSpectralAxis& targetSpectralAxis,
                               CSpectrumFluxAxis& rebinedFluxAxis, CSpectrumSpectralAxis& rebinedSpectralAxis, CMask& rebinedMask , const std::string opt_interp )
{
    if( sourceFluxAxis.GetSamplesCount() != sourceSpectralAxis.GetSamplesCount() )
    {
        return false;
    }

    if(pfgTplBuffer==0 && opt_interp=="precomputedfinegrid")
    {
        return false;
    }

    //the spectral axis should be in the same scale
    TFloat64Range logIntersectedLambdaRange( log( range.GetBegin() ), log( range.GetEnd() ) );
    TFloat64Range currentRange = logIntersectedLambdaRange;
    if( sourceSpectralAxis.IsInLinearScale() != targetSpectralAxis.IsInLinearScale() )
        return false;
    if(sourceSpectralAxis.IsInLinearScale()){
        currentRange = range;
    }

    //rebinedFluxAxis.SetSize( targetSpectralAxis.GetSamplesCount() );
    //rebinedSpectralAxis.SetSize( targetSpectralAxis.GetSamplesCount() );
    //rebinedMask.SetSize( targetSpectralAxis.GetSamplesCount() );

    const Float64* Xsrc = sourceSpectralAxis.GetSamples();
    const Float64* Ysrc = sourceFluxAxis.GetSamples();
    const Float64* Xtgt = targetSpectralAxis.GetSamples();
    Float64* Yrebin = rebinedFluxAxis.GetSamples();
    Float64* Xrebin = rebinedSpectralAxis.GetSamples();

    // Move cursors up to lambda range start
    Int32 j = 0;
    while( j<targetSpectralAxis.GetSamplesCount() && Xtgt[j] < currentRange.GetBegin() )
    {
        rebinedMask[j] = 0;
        Yrebin[j] = 0.0;
        j++;
    }
    // Include ALL lambda range
    //if( j > 0 )
    //    j--;

    if(opt_interp=="lin"){
        //* //Original interp.
        Int32 k = 0;
        Float64 t = 0.0;
        // For each sample in the valid lambda range interval.
        while( k<sourceSpectralAxis.GetSamplesCount()-1 && Xsrc[k] <= currentRange.GetEnd() )
        {
            // For each sample in the target spectrum that are in between two continous source sample
            while( j<targetSpectralAxis.GetSamplesCount() && Xtgt[j] <= Xsrc[k+1] )
            {
                // perform linear interpolation of the flux
                t = ( Xtgt[j] - Xsrc[k] ) / ( Xsrc[k+1] - Xsrc[k] );
                Xrebin[j] = Xsrc[k] + ( Xsrc[k+1] - Xsrc[k] ) * t;
                Yrebin[j] = Ysrc[k] + ( Ysrc[k+1] - Ysrc[k] ) * t;

                // closest value
                //Xrebin[j] = Xsrc[k];
                //Yrebin[j] = Ysrc[k];

                rebinedMask[j] = 1;

                j++;
            }

            k++;
        }
        //*/
    }else if(opt_interp=="precomputedfinegrid"){
        //* // Precomputed FINE GRID nearest sample, 20150801
        Int32 k = 0;
        Float64 dl = 0.1;
        Float64 Coeffk = 1.0/dl/(1+sourcez);
        // For each sample in the target spectrum
        while( j<targetSpectralAxis.GetSamplesCount() && Xtgt[j] <= currentRange.GetEnd() )
        {
            k = (int)(Xtgt[j]*Coeffk+0.5);
            // k = 10;
            Xrebin[j] = Xtgt[j];
            Yrebin[j] = pfgTplBuffer[k];
            rebinedMask[j] = 1;

            j++;

        }
        //*/
    }else if(opt_interp=="spline"){
        //* // GSL method spline
        //initialise and allocate the gsl objects
        Int32 n = sourceSpectralAxis.GetSamplesCount();
        //lin
        //gsl_interp *interpolation = gsl_interp_alloc (gsl_interp_linear,n);
        //gsl_interp_init(interpolation, Xsrc, Ysrc, n);
        //gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();
        //spline
        gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
        gsl_spline_init (spline, Xsrc, Ysrc, n);
        gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();

        // For each sample in the valid lambda range interval.
        while( j<targetSpectralAxis.GetSamplesCount() && Xtgt[j] <= currentRange.GetEnd() )
        {
            Xrebin[j] = Xtgt[j];
            //Yrebin[j] = gsl_interp_eval(interpolation, Xsrc, Ysrc, Xrebin[j], accelerator);
            Yrebin[j] = gsl_spline_eval (spline, Xrebin[j], accelerator);
            rebinedMask[j] = 1;

            j++;
        }
        //*/
    }else if(opt_interp=="ngp"){
        //* //nearest sample, lookup
        Int32 k = 0;
        Int32 kprev = 0;
        Int32 n = sourceSpectralAxis.GetSamplesCount();
        //    Int32 jmax = gsl_interp_bsearch (Xtgt, currentRange.GetEnd(), 0, targetSpectralAxis.GetSamplesCount());
        //    for(k=j; k<=jmax; k++){

        //        Xrebin[k] = Xtgt[k];
        //        Yrebin[k] = 0;
        //        rebinedMask[k] = 1;
        //    }
        // For each sample in the valid lambda range interval.
        while( j<targetSpectralAxis.GetSamplesCount() && Xtgt[j] <= currentRange.GetEnd() )
        {
            k = gsl_interp_bsearch (Xsrc, Xtgt[j], kprev, n);
            kprev = k;
            // closest value
            Xrebin[j] = Xtgt[j];
            Yrebin[j] = Ysrc[k];

            rebinedMask[j] = 1;
            j++;
        }
        //return false;
        //*/
    }

    while( j < targetSpectralAxis.GetSamplesCount() )
    {
        rebinedMask[j] = 0;
        Yrebin[j] = 0.0;
        j++;
    }

    return true;
}

Bool CSpectrumFluxAxis::RebinVarianceWeighted( const CSpectrumFluxAxis& sourceFluxAxis, const CSpectrumSpectralAxis& sourceSpectralAxis, const CSpectrumFluxAxis& sourceError,
                                               const CSpectrumSpectralAxis& targetSpectralAxis,
                                               CSpectrumFluxAxis& rebinedFluxAxis, CSpectrumSpectralAxis& rebinedSpectralAxis, CSpectrumFluxAxis& rebinedError,
                                               const std::string opt_interp )
{
    if( sourceFluxAxis.GetSamplesCount() != sourceSpectralAxis.GetSamplesCount() )
    {
        return false;
    }

    //only linear interpolation supported for now
    if(opt_interp!="lin")
    {
        return false;
    }

    TFloat64Range currentRange = targetSpectralAxis.GetLambdaRange();

    //rebinedFluxAxis.SetSize( targetSpectralAxis.GetSamplesCount() );
    //rebinedSpectralAxis.SetSize( targetSpectralAxis.GetSamplesCount() );
    //rebinedMask.SetSize( targetSpectralAxis.GetSamplesCount() );

    const Float64* Xsrc = sourceSpectralAxis.GetSamples();
    const Float64* Ysrc = sourceFluxAxis.GetSamples();
    const Float64* Errsrc = sourceError.GetSamples();
    const Float64* Xtgt = targetSpectralAxis.GetSamples();
    Float64* Yrebin = rebinedFluxAxis.GetSamples();
    Float64* Errrebin = rebinedError.GetSamples();
    Float64* Xrebin = rebinedSpectralAxis.GetSamples();

    if(opt_interp=="lin"){
        Int32 j = 0;
        Int32 k = 0;
        Float64 t = 0.0;
        //Float64 varianceCompensation=1.;
        Float64 xStepCompensation=1.;
        Float64 xSrcStep=1.;
        Float64 xDestStep=1.;
        // For each sample in the valid lambda range interval.
        while( k<sourceSpectralAxis.GetSamplesCount()-1 && Xsrc[k] <= currentRange.GetEnd() )
        {
            // For each sample in the target spectrum that are in between two continous source sample
            while( j<targetSpectralAxis.GetSamplesCount() && Xtgt[j] <= Xsrc[k+1] )
            {
                // perform linear interpolation of the flux
                xSrcStep = ( Xsrc[k+1] - Xsrc[k] );
                t = ( Xtgt[j] - Xsrc[k] ) / xSrcStep;
                //varianceCompensation = (Errsrc[k]*Errsrc[k])/(Errsrc[k+1]*Errsrc[k+1]+Errsrc[k]*Errsrc[k]); //to be verified and tested!
                //t*=varianceCompensation;

                Xrebin[j] = Xsrc[k] + xSrcStep * t;
                Yrebin[j] = Ysrc[k] + ( Ysrc[k+1] - Ysrc[k] ) * t;

                Errrebin[j] = sqrt(Errsrc[k]*Errsrc[k] + ( Errsrc[k+1]*Errsrc[k+1] - Errsrc[k]*Errsrc[k] ) * t*t);
                //*
                if(j<targetSpectralAxis.GetSamplesCount()-1)
                {
                    xDestStep = Xtgt[j+1]-Xtgt[j];
                    xStepCompensation = xSrcStep/xDestStep;
                }else if(j>0){
                    xDestStep = Xtgt[j]-Xtgt[j-1];
                    xStepCompensation = xSrcStep/xDestStep;
                }else{
                    xStepCompensation = 1.0;
                }
                Errrebin[j] *= sqrt(xStepCompensation);
                //*/

                // closest value
                //Xrebin[j] = Xsrc[k];
                //Yrebin[j] = Ysrc[k];
                j++;
            }

            k++;
        }
    }

    return true;
}

void CSpectrumFluxAxis::SetSize( UInt32 s )
{
    CSpectrumAxis::SetSize( s );
    m_StatError.resize( s );
    m_StatError.assign(s, 1.0);
}

Bool CSpectrumFluxAxis::ApplyMedianSmooth( UInt32 kernelHalfWidth )
{
    if( kernelHalfWidth == 0 )
        return false;

    if( GetSamplesCount() < ( kernelHalfWidth ) + 1 )
        return false;

    CSpectrumAxis tmp = *this;

    Int32 left = 0;
    Int32 right = 0;
    CMedian<Float64> median;
    Int32 N = GetSamplesCount();

    for( Int32 i=0; i<N; i++ )
    {
        left = max( 0, i - (Int32)kernelHalfWidth );
        right = min( (Int32)N - 1, i + (Int32)kernelHalfWidth );

        (*this)[i] = median.Find( tmp.GetSamples()+left, ( right - left ) +1 );
    }

    return true;
}

Bool CSpectrumFluxAxis::ApplyMeanSmooth( UInt32 kernelHalfWidth )
{
    if( kernelHalfWidth == 0 )
        return false;

    if( GetSamplesCount() < ( kernelHalfWidth ) + 1 )
        return false;

    CSpectrumAxis tmp = *this;

    Int32 left = 0;
    Int32 right = 0;
    CMean<Float64> mean;

    for( Int32 i=0; i<GetSamplesCount(); i++ )
    {
        left = max( 0, i - (Int32)kernelHalfWidth );
        right = min( (Int32)GetSamplesCount() - 1, i + (Int32)kernelHalfWidth );

        (*this)[i] = mean.Find( tmp.GetSamples()+left, ( right - left ) +1 );
    }

    return true;
}


Bool CSpectrumFluxAxis::ComputeMeanAndSDev( const CMask& mask, Float64& mean, Float64& sdev, const TFloat64List error ) const
{
  if( !error.empty() )
    {
        return ComputeMeanAndSDevWithError( mask, mean, sdev, error );
    }
    else
    {
        return ComputeMeanAndSDevWithoutError( mask, mean, sdev );
    }
}

Bool CSpectrumFluxAxis::ComputeMeanAndSDevWithoutError( const CMask& mask, Float64& mean,  Float64& sdev) const
{
    DebugAssert( mask.GetMasksCount() == GetSamplesCount() );

    Int32 j;
    Float64 sum,var,ep;
    Float64 ndOfSampleUsed;

    sum=0.0;
    ndOfSampleUsed=0;
    for ( j=0; j < GetSamplesCount(); j++)
    {
        DebugAssert( mask[j] == 1 || mask[j] == 0 );

        sum += mask[j] * m_Samples[j];
        ndOfSampleUsed += mask[j];
    }

    if( ndOfSampleUsed > 1 )
    {
        mean = sum / ndOfSampleUsed;

        var=0.0;
        ep = 0.0;
        for( j=0; j < GetSamplesCount() ;j++ )
        {
            sum = mask[j] * ( m_Samples[j] - mean );
            ep+=sum;
            var += sum * sum;
        }

        sdev = sqrt( (var - ep*ep / ndOfSampleUsed) / (ndOfSampleUsed-1) );
    }
    else
    {
        mean=NAN;
        sdev=NAN;
        return false;
    }

    return true;
}

Bool CSpectrumFluxAxis::ComputeMeanAndSDevWithError( const CMask& mask, Float64& mean, Float64& sdev, const TFloat64List error ) const
{
    DebugAssert( mask.GetMasksCount() == GetSamplesCount() );

    Int32 j;

    Float64 sum, var, errorSum, err;
    Int32 ndOfSampleUsed=0;

    sum=0.0;
    errorSum= 0.0;
    ndOfSampleUsed=0;
    for (j=0;j< GetSamplesCount();j++)
    {
        err = 1.0 / ( error[j] * error[j] );

        sum += mask[j] * m_Samples[j] * err;
        errorSum += mask[j] * err;

        ndOfSampleUsed += mask[j];
    }

    if (ndOfSampleUsed>1)
    {
        mean = sum / errorSum;

        var=0.0;
        for (j=0;j < GetSamplesCount();j++)
        {
            sum = mask[j] * ( m_Samples[j] - mean );

            var += sum * sum;
        }

        sdev = sqrt( var / ( ndOfSampleUsed - 1 ) );
    }
    else
    {
        mean=NAN;
        sdev=NAN;
        return false;
    }

    return true;
}


Float64  CSpectrumFluxAxis::ComputeRMSDiff( const CSpectrumFluxAxis& other )
{
    Float64 er2 = 0.f;
    Float64 er = 0.f;

    int n = GetSamplesCount();
    Float64 weight = (Float64)n;
    for (int j=0;j<n;j++)
    {
        er2 += (m_Samples[j]-other[j])*(m_Samples[j]-other[j]) / weight;

    }

    er = sqrt(er2);
    return er;
}

Bool CSpectrumFluxAxis::Subtract(const CSpectrumFluxAxis& other)
{
    Int32 N = GetSamplesCount();
    for( UInt32 i=0; i<N; i++ )
    {
        m_Samples[i] = m_Samples[i]-other[i];
    }
    return true;
}

Bool CSpectrumFluxAxis::Invert()
{
    Int32 N = GetSamplesCount();
    for( UInt32 i=0; i<N; i++ )
    {
        m_Samples[i] = -m_Samples[i];
    }
    return true;
}


