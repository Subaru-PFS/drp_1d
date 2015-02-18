#ifndef _REDSHIFT_SPECTRUM_SPECTRALAXIS_
#define _REDSHIFT_SPECTRUM_SPECTRALAXIS_

#include <epic/redshift/common/datatypes.h>
#include <epic/redshift/spectrum/axis.h>
#include <epic/core/common/range.h>

#include <vector>

namespace NSEpic
{

class CMask;

class CSpectrumSpectralAxis : public CSpectrumAxis
{

public:

    enum EFlags
    {
        nFLags_LogScale = 1 << 0
    };

    enum EShiftDirection
    {
        nShiftForward =0,
        nShiftBackward
    };

    CSpectrumSpectralAxis();
    CSpectrumSpectralAxis( UInt32 n, Bool isLogScale );
    CSpectrumSpectralAxis( const Float64* samples, UInt32 n, Bool isLogScale  );
    CSpectrumSpectralAxis( const CSpectrumSpectralAxis& origin, Float64 redshift, EShiftDirection direction );
    ~CSpectrumSpectralAxis();

    Float64             GetResolution( Float64 atWavelength = -1.0 ) const;

    Void                ShiftByWaveLength(  const CSpectrumSpectralAxis& origin, Float64 wavelengthOffset, EShiftDirection direction );
    Void                ShiftByWaveLength( Float64 wavelengthOffset, EShiftDirection direction );

    Int32               GetIndexAtWaveLength( Float64 waveLength ) const;
    TInt32Range         GetIndexesAtWaveLengthRange( const TFloat64Range& waveLengthRange ) const;


    Bool                ConvertToLinearScale();
    Bool                ConvertToLogScale();
    Bool                IsInLogScale() const;
    Bool                IsInLinearScale() const;

    TLambdaRange        GetLambdaRange() const;
    Bool                ClampLambdaRange( const TFloat64Range& range, TFloat64Range& clampedRange ) const;
    Void                GetMask( const TFloat64Range& range,  CMask& mask ) const;

    Bool                PlotResolution( const char* filePath ) const;

    Void                CopyFrom( const CSpectrumSpectralAxis& other );

private:

    UInt32              m_SpectralFlags;
};

}

#endif
