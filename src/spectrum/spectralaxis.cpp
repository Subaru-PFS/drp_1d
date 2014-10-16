#include <epic/redshift/spectrum/spectralaxis.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/common/mask.h>

#include <math.h>

using namespace __NS__;
using namespace std;

CSpectrumSpectralAxis::CSpectrumSpectralAxis() :
    m_SpectralFlags( 0 )
{

}

CSpectrumSpectralAxis::CSpectrumSpectralAxis( UInt32 n ) :
    CSpectrumAxis( n ),
    m_SpectralFlags( 0 )
{
}

CSpectrumSpectralAxis::CSpectrumSpectralAxis( const Float64* samples, UInt32 n ) :
    CSpectrumAxis( samples, n ),
    m_SpectralFlags( 0 )
{
}

CSpectrumSpectralAxis::~CSpectrumSpectralAxis()
{

}

Void CSpectrumSpectralAxis::ShiftByWaveLength( Float64 wavelengthOffset )
{
    if( m_SpectralFlags & nFLags_LogScale )
        wavelengthOffset = log( wavelengthOffset );

    for( Int32 i=0; i< GetSamplesCount(); i++ )
    {
        m_Samples[i] = m_Samples[i] + wavelengthOffset;
    }
}

Float64 CSpectrumSpectralAxis::GetResolution( Float64 atWavelength ) const
{
    if( GetSamplesCount() < 2 )
        return 0.0;

    Int32 i = 0;

    if( atWavelength > -1.0 )
    {
        i = GetIndexAtWaveLength( atWavelength );

        if( i >= m_Samples.size()-1 )
            i = m_Samples.size()-2;

        return m_Samples[i+1] - m_Samples[i];
    }
    else
    {
        return GetLambdaRange().GetLength() / m_Samples.size();
    }

    return 0;
}

Bool CSpectrumSpectralAxis::IsInLogScale() const
{
    return m_SpectralFlags & nFLags_LogScale;
}

Bool CSpectrumSpectralAxis::IsInLinearScale() const
{
    return !(m_SpectralFlags & nFLags_LogScale);
}

Bool CSpectrumSpectralAxis::HasUniformResolution() const
{
    Float64 r1 = GetResolution();
    Float64 r2 = 0.0;
    Float64 d = 0.0;

    for( Int32 i=0; i<GetSamplesCount()-1; i++ )
    {
        r2 = m_Samples[i+1] - m_Samples[i];
        d = fabs( r2 - r1 );
        if( d > 0.01 )
            return false;
    }

    return true;
}

TLambdaRange CSpectrumSpectralAxis::GetLambdaRange() const
{
    if( m_Samples.size() < 2 )
        return TLambdaRange( 0.0, 0.0 );

    if( m_SpectralFlags & nFLags_LogScale )
    {
        return TLambdaRange( exp( m_Samples[0] ), exp( m_Samples[m_Samples.size()-1] ) );
    }

    return TLambdaRange( m_Samples[0], m_Samples[m_Samples.size()-1] );
}

Void CSpectrumSpectralAxis::GetMask( const TFloat64Range& lambdaRange,  CMask& mask ) const
{
    TFloat64Range range = lambdaRange;

    if( m_SpectralFlags & nFLags_LogScale )
        range.Set( log( range.GetBegin() ), log( range.GetEnd() ) );

    // weight = Spectrum over lambdarange flag
    for( Int32 i=0; i< m_Samples.size(); i++ )
    {
        mask[i] = (Mask)0;
        // If this sample is somewhere in a valid lambdaRande, tag weight with 1

        if( m_Samples[i] >= range.GetBegin() && m_Samples[i] <= range.GetEnd() )
        {
            mask[i]= (Mask)1;
        }

    }
}

Bool CSpectrumSpectralAxis::ClampLambdaRange( const TFloat64Range& range, TFloat64Range& clampedRange ) const
{
    TFloat64Range otherRange = GetLambdaRange();

    clampedRange = TFloat64Range();

    if( range.GetIsEmpty() || otherRange.GetIsEmpty() )
        return false;

    // Clamp lambda start
    Float64 start = range.GetBegin();
    if ( start < otherRange.GetBegin() )
        start = otherRange.GetBegin();

    // CLamp lambda end
    Float64 end = range.GetEnd();
    if ( end > otherRange.GetEnd() )
        end = otherRange.GetEnd();

    clampedRange = TFloat64Range( start, end );
    return true;
}

Bool CSpectrumSpectralAxis::PlotResolution( const char* filePath ) const
{
    FILE* f = fopen( filePath, "w+" );
    if( f == NULL )
        return false;

    if( m_Samples.size() >= 2 )
    {
        for( int i=0;i<m_Samples.size()-1;i++)
        {
            fprintf( f, "%f %f\n", m_Samples[i], m_Samples[i+1] - m_Samples[i] );
        }
    }
    fclose( f );

    return true;
}

TInt32Range CSpectrumSpectralAxis::GetIndexesAtWaveLengthRange( const TFloat64Range& waveLengthRange ) const
{
    TInt32Range r;

    r.SetBegin( GetIndexAtWaveLength( waveLengthRange.GetBegin() ) );
    r.SetEnd( GetIndexAtWaveLength( waveLengthRange.GetEnd() ));

    return r;
}

Int32 CSpectrumSpectralAxis::GetIndexAtWaveLength( Float64 waveLength ) const
{
    Int32 m;
    Int32 lo = 0;
    Int32 hi = GetSamplesCount()-1;

    if( waveLength <= m_Samples[lo] )
        return lo;

    if( waveLength >= m_Samples[hi] )
        return hi;

    for (;;)
    {

        m = (lo + hi) / 2;

        if( waveLength < m_Samples[m] )
            hi = m - 1;
        else if( waveLength > m_Samples[m] )
            lo = m + 1;
        else
            return m;

        if (lo >= hi)
            return(lo);
    }

    return -1;
}

Bool CSpectrumSpectralAxis::ConvertToLinearScale()
{
    if( ( m_SpectralFlags & nFLags_LogScale ) == false )
        return true;

    for( Int32 i=0; i<GetSamplesCount(); i++ )
    {
        m_Samples[i] = exp( m_Samples[i] );
    }

    m_SpectralFlags &= ~ nFLags_LogScale;

    return true;
}

Bool CSpectrumSpectralAxis::ConvertToLogScale()
{
    if( m_SpectralFlags & nFLags_LogScale )
        return true;

    for( Int32 i=0; i<GetSamplesCount(); i++ )
    {
        m_Samples[i] = log( m_Samples[i] );
    }

    m_SpectralFlags |= nFLags_LogScale;

    return true;
}
