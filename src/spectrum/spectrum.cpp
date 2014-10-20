#include <epic/redshift/spectrum/spectrum.h>

#include <epic/core/serializer/serializer.h>
#include <epic/core/debug/assert.h>

#include <math.h>
#include <stdio.h>

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CSpectrum )

CSpectrum::CSpectrum()
{

}

CSpectrum::~CSpectrum()
{

}

/**
 * Convert the spectral axis to a neperian logarithm scale
 */
Bool CSpectrum::ConvertToLogScale()
{
    return m_SpectralAxis.ConvertToLogScale();
}

/**
 * Convert the spectral axis to a linear scale
 */
Bool CSpectrum::ConvertToLinearScale()
{
    return m_SpectralAxis.ConvertToLinearScale();
}

Float64 CSpectrum::GetResolution() const
{
    return m_SpectralAxis.GetResolution();
}

/**
 * Return the lambda range of the entire spectrum.
 * Range is always expressed in linear scale NOT in log scale even if the underlying spcetrum is in log scale
 */
TLambdaRange CSpectrum::GetLambdaRange() const
{
    return m_SpectralAxis.GetLambdaRange();
}


Bool CSpectrum::Serialize( CSerializer& ar )
{
    Int16 version = 1;

    if( ar.BeginScope( "Spectrum", version ) == version )
    {
        ar.Serialize( m_FluxAxis, "FluxAxis" );
        ar.Serialize( m_SpectralAxis, "SpectralAxis" );
        ar.EndScope();
        return true;
    }

    return false;
}

