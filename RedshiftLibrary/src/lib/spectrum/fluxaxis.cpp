#include <RedshiftLibrary/spectrum/fluxaxis.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/common/median.h>
#include <RedshiftLibrary/common/mean.h>
#include <RedshiftLibrary/common/mask.h>

#include <math.h>
#include <RedshiftLibrary/log/log.h>

using namespace NSEpic;
using namespace std;



CSpectrumFluxAxis::CSpectrumFluxAxis( UInt32 n, Float64 value ) :
    CSpectrumAxis( n, value),
    m_StdError(n)
{

}

CSpectrumFluxAxis::CSpectrumFluxAxis(const CSpectrumAxis & otherFlux, const CSpectrumNoiseAxis & otherError ):
   CSpectrumAxis(otherFlux),
   m_StdError(otherError)
{}


CSpectrumFluxAxis::CSpectrumFluxAxis( const Float64* samples, UInt32 n ) :
    CSpectrumAxis( samples, n ),
    m_StdError(n)
{

}
CSpectrumFluxAxis::CSpectrumFluxAxis( const TFloat64List samples, UInt32 n ) :
    CSpectrumAxis( samples, n ),
    m_StatError( n, 1.0 )
{

}

CSpectrumFluxAxis::CSpectrumFluxAxis( const Float64* samples, UInt32 n,
				      const Float64* error,  const UInt32 m) :
    CSpectrumAxis( samples, n ),
    m_StdError( error, m )
{
    if(m!=n){
        Log.LogError("FluxAxis and NoiseAxis do not have equal size.");
        throw runtime_error("FluxAxis and NoiseAxis do not have equal size.");
    }
}

void CSpectrumFluxAxis::SetSize( UInt32 s )
{
    CSpectrumAxis::SetSize( s );
    m_StdError.SetSize(s);
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


Bool CSpectrumFluxAxis::ComputeMeanAndSDev( const CMask& mask, Float64& mean, Float64& sdev, const CSpectrumNoiseAxis error ) const
{
  if( !error.isEmpty() )
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

Bool CSpectrumFluxAxis::ComputeMeanAndSDevWithError( const CMask& mask, Float64& mean, Float64& sdev, const CSpectrumNoiseAxis error ) const
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
