#include <epic/redshift/spectrum/spectralaxis.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/common/mask.h>

#include <math.h>

using namespace NSEpic;
using namespace std;

CSpectrumSpectralAxis::CSpectrumSpectralAxis() :
    m_SpectralFlags( 0 )
{

}

CSpectrumSpectralAxis::CSpectrumSpectralAxis( UInt32 n, Bool isLogScale ) :
    CSpectrumAxis( n ),
    m_SpectralFlags( 0 )
{
    if( isLogScale )
        m_SpectralFlags |= nFLags_LogScale;
}

CSpectrumSpectralAxis::CSpectrumSpectralAxis( const Float64* samples, UInt32 n, Bool isLogScale ) :
    CSpectrumAxis( samples, n ),
    m_SpectralFlags( 0 )
{
    if( isLogScale )
        m_SpectralFlags |= nFLags_LogScale;
}

CSpectrumSpectralAxis::~CSpectrumSpectralAxis()
{

}


CSpectrumSpectralAxis::CSpectrumSpectralAxis( const CSpectrumSpectralAxis& origin, Float64 wavelengthOffset, EShiftDirection direction  )
{
    ShiftByWaveLength( origin, wavelengthOffset, direction );
}

Void CSpectrumSpectralAxis::ShiftByWaveLength( const CSpectrumSpectralAxis& origin, Float64 wavelengthOffset, EShiftDirection direction )
{
    m_SpectralFlags = 0;

    DebugAssert( origin.GetSamplesCount() == GetSamplesCount() );

    SetSize( origin.GetSamplesCount() );
    const Float64* originSamples = origin.GetSamples();

    DebugAssert( direction == nShiftForward || direction == nShiftBackward );

    if( origin.IsInLogScale() )
    {
        wavelengthOffset = log( wavelengthOffset );
        m_SpectralFlags |= nFLags_LogScale;

        if( direction == nShiftForward )
        {
            for( Int32 i=0; i< origin.GetSamplesCount(); i++ )
            {
                m_Samples[i] = originSamples[i] + wavelengthOffset;
            }
        }
        else if( direction == nShiftBackward )
        {
            for( Int32 i=0; i< origin.GetSamplesCount(); i++ )
            {
                m_Samples[i] = originSamples[i] - wavelengthOffset;
            }
        }
    }
    else
    {
        if( direction == nShiftForward )
        {
            for( Int32 i=0; i< origin.GetSamplesCount(); i++ )
            {
                m_Samples[i] = originSamples[i] * wavelengthOffset;
            }
        }
        else if( direction == nShiftBackward )
        {
            for( Int32 i=0; i< origin.GetSamplesCount(); i++ )
            {
                m_Samples[i] = originSamples[i] / wavelengthOffset;
            }
        }
    }

}

Void CSpectrumSpectralAxis::CopyFrom( const CSpectrumSpectralAxis& other )
{
    m_SpectralFlags = other.m_SpectralFlags;

    SetSize( other.GetSamplesCount() );

    // Copy spectral data
    const Float64* otherData = other.GetSamples();
    for( UInt32 i=0; i<other.GetSamplesCount(); i++ )
    {
        m_Samples[i] = otherData[i];
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
    TFloat64Range effectiveRange = GetLambdaRange();

    clampedRange = TFloat64Range();

    if( range.GetIsEmpty() || effectiveRange.GetIsEmpty() )
        return false;

    // Clamp lambda start
    Float64 start = range.GetBegin();
    if ( start < effectiveRange.GetBegin() )
        start = effectiveRange.GetBegin();

    // CLamp lambda end
    Float64 end = range.GetEnd();
    if ( end > effectiveRange.GetEnd() )
        end = effectiveRange.GetEnd();

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