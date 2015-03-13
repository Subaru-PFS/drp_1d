#include <epic/redshift/spectrum/fluxaxis.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/common/median.h>
#include <epic/redshift/common/mean.h>
#include <epic/redshift/common/mask.h>

#include <math.h>

#include <epic/core/log/log.h>

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

CSpectrumFluxAxis::~CSpectrumFluxAxis()
{

}

CSpectrumFluxAxis& CSpectrumFluxAxis::operator=(const CSpectrumFluxAxis& other)
{
    m_StatError = other.m_StatError;
    CSpectrumAxis::operator=( other );
}

Void CSpectrumFluxAxis::SetSize( UInt32 s )
{
    CSpectrumAxis::SetSize( s );
    m_StatError.resize( s );

    for( UInt32 i=0; i<s; i++ )
    {
        m_StatError[i] = 1.0;
    }
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


Bool CSpectrumFluxAxis::ComputeMeanAndSDev( const CMask& mask, Float64& mean, Float64& sdev, const Float64* error ) const
{
    if( error )
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
    Float64 sum,var;
    Int32 ndOfSampleUsed;

    sum=0.0;
    ndOfSampleUsed=0;
    for ( j=0; j < GetSamplesCount(); j++)
    {
        DebugAssert( mask[j] == 1 || mask[j] == 0 );

        sum += mask[j] * m_Samples[j];
        ndOfSampleUsed += mask[j];
    }


    if( ndOfSampleUsed > 0 )
    {
        mean = sum / ndOfSampleUsed;

        var=0.0;
        for( j=0; j < GetSamplesCount() ;j++ )
        {
            sum = mask[j] * ( m_Samples[j] - mean );
            var += sum * sum;
        }

        sdev = sqrt( ( 1.0 / ndOfSampleUsed ) * var );
    }
    else
    {
        mean=NAN;
        sdev=NAN;
        return false;
    }

    return true;
}

Bool CSpectrumFluxAxis::ComputeMeanAndSDevWithError( const CMask& mask, Float64& mean, Float64& sdev, const Float64* error ) const
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
        err = error[j] * error[j];
        sum += mask[j] * m_Samples[j] / err;
        //errorSum += mask[j] / err;

        errorSum += mask[j] * err;
        ndOfSampleUsed += mask[j];
    }

    if (ndOfSampleUsed>0)
    {
        mean = sum / errorSum;

        var=0.0;
        for (j=0;j < GetSamplesCount();j++)
        {
            sum= mask[j] * ( m_Samples[j] - mean );
            var += sum * sum;
        }

        sdev = sqrt( ( 1.0 / ndOfSampleUsed ) * var );
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
}


