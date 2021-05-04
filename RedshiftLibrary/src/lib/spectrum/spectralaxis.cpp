#include <RedshiftLibrary/spectrum/spectralaxis.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/common/mask.h>
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include <cmath>
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
CSpectrumSpectralAxis::CSpectrumSpectralAxis( const TFloat64List samples, Bool isLogScale ) :
    CSpectrumAxis(std::move(samples)),
    m_SpectralFlags( 0 )
{
    if( isLogScale )
        m_SpectralFlags |= nFLags_LogScale;
}

/**
 * Constructor.
 */
CSpectrumSpectralAxis::CSpectrumSpectralAxis( const TFloat64List samples) :
    CSpectrumAxis(std::move(samples)),
    m_SpectralFlags( 0 )
{
}
//only used by client
CSpectrumSpectralAxis::CSpectrumSpectralAxis( const Float64* samples, UInt32 n) :
    CSpectrumAxis( samples, n ),
    m_SpectralFlags( 0 )
{
}

CSpectrumSpectralAxis::CSpectrumSpectralAxis(const CSpectrumSpectralAxis & other):
    CSpectrumAxis(other),
    m_SpectralFlags(other.m_SpectralFlags),
    m_regularLogSamplingStep(other.m_regularLogSamplingStep),
    m_regularLogSamplingChecked(other.m_regularLogSamplingChecked)
{}
CSpectrumSpectralAxis::CSpectrumSpectralAxis(CSpectrumSpectralAxis && other):
    CSpectrumAxis(other),
    m_SpectralFlags(other.m_SpectralFlags),
    m_regularLogSamplingStep(other.m_regularLogSamplingStep),
    m_regularLogSamplingChecked(other.m_regularLogSamplingChecked)
{}

CSpectrumSpectralAxis& CSpectrumSpectralAxis::operator=(const CSpectrumSpectralAxis& other)
{
    m_Samples = other.m_Samples;
    m_SpectralFlags  = other.m_SpectralFlags;
    m_regularLogSamplingStep = other.m_regularLogSamplingStep;
    m_regularLogSamplingChecked = other.m_regularLogSamplingChecked;
    return *this;
}
CSpectrumSpectralAxis& CSpectrumSpectralAxis::operator=( CSpectrumSpectralAxis&& other)
{
    m_Samples = other.m_Samples;
    m_SpectralFlags  = other.m_SpectralFlags;
    m_regularLogSamplingStep = other.m_regularLogSamplingStep;
    m_regularLogSamplingChecked = other.m_regularLogSamplingChecked;
    return *this;
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

void CSpectrumSpectralAxis::SetLogScale() 
{
    m_SpectralFlags |= nFLags_LogScale; //not sure what is the difference between these two flags      
}

/**
 * Check if spectralAxis is well rebinned in log
*/
Bool CSpectrumSpectralAxis::CheckLoglambdaSampling() const
{
    Float64 logGridStep;
    if(IsInLogScale())
        logGridStep = m_Samples[1] - m_Samples[0];
    else
        logGridStep = log(m_Samples[1]/m_Samples[0]);

    Float64 relativelogGridStepTol = 1e-1;
    Float64 maxAbsRelativeError = 0.0; 
    Float64 lbda1 = m_Samples[0];
    for (Int32 t = 1; t < m_Samples.size(); t++)
    {
        Float64 lbda2 =  m_Samples[t];
        Float64 _logGridStep;
        if(IsInLogScale())
            _logGridStep = lbda2 - lbda1;
        else
            _logGridStep = log(lbda2/lbda1);

        Float64 relativeErrAbs = std::abs((_logGridStep - logGridStep) / logGridStep);
        maxAbsRelativeError = max(relativeErrAbs, maxAbsRelativeError);

        if (relativeErrAbs > relativelogGridStepTol)
        {//return without setting anything
            m_SpectralFlags &= ~ nFLags_LogSampled; // unsetting
            Log.LogDebug("   CSpectrumSpectralAxis::CheckLoglambdaSampling: Log-regular lambda check FAILED");

            return false;
        }
        lbda1 = lbda2;
    }
    //save step in a member variable
    m_regularLogSamplingStep = logGridStep; 
    m_SpectralFlags |= nFLags_LogSampled; //setting through a logical OR
    m_regularLogSamplingChecked = true;
    Log.LogDetail("   CSpectrumSpectralAxis::CheckLoglambdaSampling: max Abs Relative Error (log lbda step)= %f", maxAbsRelativeError);

    return true;
}

/**
 * Brief:
 * In the actual version we consider that a spectral axis can be log sampled while having its values in Angstrom (i.e., non-log)
 * 
*/
Bool CSpectrumSpectralAxis::IsLogSampled(Float64 logGridstep) const
{       
    if (!IsLogSampled() )
        return false;

    if (std::abs(m_regularLogSamplingStep - logGridstep)>1E-8)
    {
        Log.LogDetail("   CSpectrumSpectralAxis::IsLogSampled: Log-regular sampling with bad step");
        return false;
    }
    
    return true;
}

Bool CSpectrumSpectralAxis::IsLogSampled() const
{
    if(!m_regularLogSamplingChecked) 
        CheckLoglambdaSampling();

    return m_SpectralFlags & nFLags_LogSampled;
}


Float64 CSpectrumSpectralAxis::GetlogGridStep() const
{   
    if (!IsLogSampled())
    {
        Log.LogError("CSpectrumSpectralAxis::GetlogGridStep: axis is not logsampled");
        throw runtime_error("CSpectrumSpectralAxis::GetlogGridStep: axis is not logsampled");
    }

    return m_regularLogSamplingStep;
 }

//still TODO: check end-to-end redshift coverage
TFloat64List CSpectrumSpectralAxis::GetSubSamplingMask(UInt32 ssratio) const
{
    return GetSubSamplingMask(ssratio, TInt32Range(0, GetSamplesCount()-1));
}

/*@ssratio stands for sub-samplingRatio*/
TFloat64List CSpectrumSpectralAxis::GetSubSamplingMask(UInt32 ssratio, TInt32Range ilbda) const
{
    //if(std::isnan(m_regularLogSamplingStep))
    if(!IsLogSampled())
    {
        throw runtime_error("Cannot subsample spectrum!");
    }
    UInt32 s = GetSamplesCount();
    if(ssratio==1) return TFloat64List(s, 1.);
    TFloat64List mask(s, 0.); 
    //for(Int32 i = ilbda.GetBegin(); i>=ilbda.GetEnd(); i+=multi){ //here no need to reverse
    for(Int32 i=ilbda.GetEnd(); i>=ilbda.GetBegin();i-=ssratio){//ensure that z[0] remains the same
        mask[i]=1;
    }
    return mask;
}

/**
 * Brief: 
 * check that division of two floating values gives an int, and return modulo value
 *  Int32 logstep_int = Int32(logstep*1E12);
    Int32 rebinlogstep_int = Int32(rebinlogstep*1E12);
    auto lambda_redshift_modulo = logstep_int %rebinlogstep_int;
    auto modulo_2 = logstep_int - trunc(logstep_int/rebinlogstep_int)*rebinlogstep_int;
*/
UInt32 CSpectrumSpectralAxis::GetLogSamplingIntegerRatio(Float64 logstep, Float64& modulo) const
{
    if(!IsLogSampled())
    {
        throw runtime_error("  CSpectrumSpectralAxis::GetIntegerRatio: axis is not logsampled, thus cannot get integer ratio ");
    }

    UInt32 ratio = std::round(logstep/m_regularLogSamplingStep);
    modulo = logstep - ratio*m_regularLogSamplingStep;
    return ratio;
}