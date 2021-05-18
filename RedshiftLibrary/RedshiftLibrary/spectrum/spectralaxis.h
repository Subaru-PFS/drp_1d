#ifndef _REDSHIFT_SPECTRUM_SPECTRALAXIS_
#define _REDSHIFT_SPECTRUM_SPECTRALAXIS_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/common/range.h>

#include <vector>

namespace NSEpic
{

class CMask;

/**
 * \ingroup Redshift
 */
class CSpectrumSpectralAxis : public CSpectrumAxis
{

public:

    enum EFlags
    {
        nFLags_LogScale = 1 << 0,
        nFLags_LogSampled = 1 << 1 //using second bit(shifting once): 0000 0010
    };

    enum EShiftDirection
    {
        nShiftForward =0,
        nShiftBackward
    };

    CSpectrumSpectralAxis();
    CSpectrumSpectralAxis(const CSpectrumSpectralAxis & other) = default;
    CSpectrumSpectralAxis(CSpectrumSpectralAxis && other) = default;
    CSpectrumSpectralAxis( UInt32 n, Bool isLogScale =false);
    CSpectrumSpectralAxis( const TFloat64List & samples, Bool isLogScale=false  );
    CSpectrumSpectralAxis( TFloat64List && samples, Bool isLogScale=false  );
    CSpectrumSpectralAxis( const Float64* samples, UInt32 n);
    CSpectrumSpectralAxis( const CSpectrumSpectralAxis& origin, Float64 redshift, EShiftDirection direction );
    ~CSpectrumSpectralAxis() = default;

    CSpectrumSpectralAxis& operator=(const CSpectrumSpectralAxis& other) = default;
    CSpectrumSpectralAxis& operator=( CSpectrumSpectralAxis&& other) = default;

    Float64             GetResolution( Float64 atWavelength = -1.0 ) const;
    Float64             GetMeanResolution() const;

    void                ShiftByWaveLength(  const CSpectrumSpectralAxis& origin, Float64 wavelengthOffset, EShiftDirection direction );
    void                ShiftByWaveLength( Float64 wavelengthOffset, EShiftDirection direction );

    void                ApplyOffset(Float64 wavelengthOffset);

    Int32               GetIndexAtWaveLength( Float64 waveLength ) const;
    TInt32Range         GetIndexesAtWaveLengthRange( const TFloat64Range& waveLengthRange ) const;


    Bool                ConvertToLinearScale();
    Bool                ConvertToLogScale();
    Bool                IsInLogScale() const;
    Bool                IsInLinearScale() const;

    TLambdaRange        GetLambdaRange() const;
    Bool                ClampLambdaRange( const TFloat64Range& range, TFloat64Range& clampedRange ) const;
    void                GetMask( const TFloat64Range& range,  CMask& mask ) const;
    Float64             IntersectMaskAndComputeOverlapRate( const TFloat64Range& lambdaRange,  CMask& omask ) const;
    void                SetLogScale();
    Bool                CheckLoglambdaSampling()const;
    Bool                IsLogSampled(Float64 logGridstep)const;
    Bool                IsLogSampled()const;
    Float64             GetlogGridStep() const;
    void                RecomputePreciseLoglambda();
    TFloat64List        GetSubSamplingMask(UInt32 ssratio) const;
    TFloat64List        GetSubSamplingMask(UInt32 ssratio, TFloat64Range lambdarange) const;
    TFloat64List        GetSubSamplingMask(UInt32 ssratio, const TInt32Range & ilbda) const;
    UInt32              GetLogSamplingIntegerRatio(Float64 logstep, Float64& modulo) const;
private:

    mutable UInt32      m_SpectralFlags = 0;
    mutable Bool        m_regularLogSamplingChecked=false;
    mutable Float64     m_regularLogSamplingStep = NAN; //sampling log step with which sampling was validated in CheckLoglambdaSampling 

};

}

#endif
