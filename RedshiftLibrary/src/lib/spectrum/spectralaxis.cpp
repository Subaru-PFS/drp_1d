#include <RedshiftLibrary/spectrum/spectralaxis.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/common/mask.h>

#include <math.h>

using namespace NSEpic;
using namespace std;

/**
 * Constructor, zeroes flags.
 */
CSpectrumSpectralAxis::CSpectrumSpectralAxis() :
    m_SpectralFlags( 0 )
{

}

/**
 * Constructor, flags log scale when set.
 */
CSpectrumSpectralAxis::CSpectrumSpectralAxis( UInt32 n, Bool isLogScale ) :
    CSpectrumAxis( n ),
    m_SpectralFlags( 0 )
{
    if( isLogScale )
        m_SpectralFlags |= nFLags_LogScale;
}

/**
 * Constructor, flags log scale when set.
 */
CSpectrumSpectralAxis::CSpectrumSpectralAxis( const Float64* samples, UInt32 n, Bool isLogScale ) :
    CSpectrumAxis( samples, n ),
    m_SpectralFlags( 0 )
{
    if( isLogScale )
        m_SpectralFlags |= nFLags_LogScale;
}

/**
 * Constructor.
 */
CSpectrumSpectralAxis::CSpectrumSpectralAxis( const Float64* samples, UInt32 n) :
    CSpectrumAxis( samples, n ),
    m_SpectralFlags( 0 )
{
}

CSpectrumSpectralAxis::CSpectrumSpectralAxis( const TFloat64List samples, UInt32 n) :
    CSpectrumAxis( samples, n ),
    m_SpectralFlags( 0 )
{
}

/**
 * Constructor, shifts origin along direction an offset distance.
 */
CSpectrumSpectralAxis::CSpectrumSpectralAxis( const CSpectrumSpectralAxis& origin, Float64 wavelengthOffset, EShiftDirection direction  ) :
  CSpectrumAxis( origin.GetSamplesCount() )
{
    ShiftByWaveLength( origin, wavelengthOffset, direction );
}

/**
 * Shift current axis the input offset in the input direction.
 */
void CSpectrumSpectralAxis::ShiftByWaveLength( Float64 wavelengthOffset, EShiftDirection direction )
{
	ShiftByWaveLength( *this, wavelengthOffset, direction );
}

/**
 * Copy the input axis samples and shift the axis the specified offset in the specidifed direction.
 */
void CSpectrumSpectralAxis::ShiftByWaveLength( const CSpectrumSpectralAxis& origin, Float64 wavelengthOffset, EShiftDirection direction )
{
    Int32 nSamples = origin.GetSamplesCount();
    m_SpectralFlags = 0;
    m_Samples.resize(nSamples);

    const Float64* originSamples = origin.GetSamples();

    DebugAssert( direction == nShiftForward || direction == nShiftBackward );

    if( wavelengthOffset == 0.0 )
    {
        for( Int32 i=0; i<nSamples; i++ )
        {
            m_Samples[i] = originSamples[i];
        }

        return;
    }

    if( origin.IsInLogScale() )
    {
        wavelengthOffset = log( wavelengthOffset );
        m_SpectralFlags |= nFLags_LogScale;

        if( direction == nShiftForward )
        {
            for( Int32 i=0; i<nSamples; i++ )
            {
                m_Samples[i] = originSamples[i] + wavelengthOffset;
            }
        }
        else if( direction == nShiftBackward )
        {
            for( Int32 i=0; i<nSamples; i++ )
            {
                m_Samples[i] = originSamples[i] - wavelengthOffset;
            }
        }
    }
    else
    {

        if( direction == nShiftForward )
        {
            for( Int32 i=0; i<nSamples; i++ )
            {
                m_Samples[i] = originSamples[i] * wavelengthOffset;
            }
        }
        else if( direction == nShiftBackward )
        {
            for( Int32 i=0; i<nSamples; i++ )
            {
                m_Samples[i] = originSamples[i] / wavelengthOffset;
            }
        }
    }

}

void CSpectrumSpectralAxis::ApplyOffset( Float64 wavelengthOffset )
{
    Int32 nSamples = m_Samples.size();
    for( Int32 i=0; i<nSamples ; i++ )
    {
        m_Samples[i] += wavelengthOffset;
    }

}

/**
 * Return the wavelength interval between two consecutive samples.
 */
Float64 CSpectrumSpectralAxis::GetResolution( Float64 atWavelength ) const
{
    if( GetSamplesCount() < 2 )
        return 0.0;

    Int32 i = 0;

    if( atWavelength >= 0.0 )
    {
        i = GetIndexAtWaveLength( atWavelength );

        if( i > m_Samples.size()-1 )
            i = m_Samples.size()-1;
        if( i < 1 )
            i = 1;

        return m_Samples[i] - m_Samples[i-1];
    }
    else
    {
        return m_Samples[1] - m_Samples[0];
    }

    return 0.0;
}

/**
 *
 */
Float64 CSpectrumSpectralAxis::GetMeanResolution() const
{
    if( GetSamplesCount() < 2 )
        return 0.0;

    Float64 resolution = 0.0;
    Int32 nsum = 0;
    Int32 nSamples = m_Samples.size()-1;
    for( Int32 i=0; i<nSamples; i++ )
    {
        resolution += (m_Samples[i+1]-m_Samples[i]);
        nsum++;
    }

    if(nsum>0){
        resolution /=nsum;
    }

    return resolution;
}

/**
 *
 */
Bool CSpectrumSpectralAxis::IsInLogScale() const
{
    return m_SpectralFlags & nFLags_LogScale;
}

/**
 *
 */
Bool CSpectrumSpectralAxis::IsInLinearScale() const
{
    return !(m_SpectralFlags & nFLags_LogScale);
}

/**
 *
 */
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

/**
 *
 */
void CSpectrumSpectralAxis::GetMask( const TFloat64Range& lambdaRange,  CMask& mask ) const
{
    TFloat64Range range = lambdaRange;

    if( m_SpectralFlags & nFLags_LogScale )
        range.Set( log( range.GetBegin() ), log( range.GetEnd() ) );

    mask.SetSize( m_Samples.size() );

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

/**
 *
 */
Float64 CSpectrumSpectralAxis::IntersectMaskAndComputeOverlapRate( const TFloat64Range& lambdaRange,  CMask& omask ) const
{
    TFloat64Range range = lambdaRange;

    if( m_SpectralFlags & nFLags_LogScale )
        range.Set( log( range.GetBegin() ), log( range.GetEnd() ) );

    Int32 selfRate=0;
    Int32 otherRate=0;
    const Mask* otherWeight = omask.GetMasks();
    Int32 nSamples = m_Samples.size();
    // weight = Spectrum over lambdarange flag
    for( Int32 i=0; i<nSamples; i++ )
    {
        //otherWeight[i] = 0;
        // If this sample is somewhere in a valid lambdaRande, tag weight with 1
        if( m_Samples[i] >= range.GetBegin() && m_Samples[i] <= range.GetEnd() )
        {
            if(otherWeight[i]){
                otherRate++;
            }
            selfRate++;
        }
    }

    if( selfRate == 0.0 )
        return 0;

    return (Float64)otherRate/(Float64)selfRate;
}

/**
 *
 */
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


/**
 *
 */
TInt32Range CSpectrumSpectralAxis::GetIndexesAtWaveLengthRange( const TFloat64Range& waveLengthRange ) const
{
    TInt32Range r;

    r.SetBegin( GetIndexAtWaveLength( waveLengthRange.GetBegin() ) );
    r.SetEnd( GetIndexAtWaveLength( waveLengthRange.GetEnd() ));

    return r;
}

/**
 *
 */
Int32 CSpectrumSpectralAxis::GetIndexAtWaveLength( Float64 waveLength ) const
{
    Int32 m;
    Int32 lo = 0;
    Int32 hi = GetSamplesCount()-1;

    if( waveLength <= m_Samples[lo] )
        return lo;

    if( waveLength >= m_Samples[hi] )
        return hi;

    auto it = std::lower_bound(m_Samples.begin(),m_Samples.end(),waveLength);

    return (it - m_Samples.begin());

}

/**
 *
 */
Bool CSpectrumSpectralAxis::ConvertToLinearScale()
{
    if( ( m_SpectralFlags & nFLags_LogScale ) == false )
        return true;
    Int32 nSamples = GetSamplesCount();

    for( Int32 i=0; i<nSamples; i++ )
    {
        m_Samples[i] = exp( m_Samples[i] );
    }

    m_SpectralFlags &= ~ nFLags_LogScale;

    return true;
}

/**
 *
 */
Bool CSpectrumSpectralAxis::ConvertToLogScale()
{
    if( m_SpectralFlags & nFLags_LogScale )
        return true;
    Int32 nSamples = GetSamplesCount();

    for( Int32 i=0; i<nSamples; i++ )
    {
        m_Samples[i] = log( m_Samples[i] );
    }

    m_SpectralFlags |= nFLags_LogScale;

    return true;
}
