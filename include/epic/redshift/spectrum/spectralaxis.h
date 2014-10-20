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

    CSpectrumSpectralAxis();
    CSpectrumSpectralAxis( UInt32 n );
    CSpectrumSpectralAxis( const Float64* samples, UInt32 n );
    ~CSpectrumSpectralAxis();

    Float64             GetResolution( Float64 atWavelength = -1.0 ) const;

    Void                ShiftByWaveLength( Float64 redShift );

    Int32               GetIndexAtWaveLength( Float64 waveLength ) const;
    TInt32Range         GetIndexesAtWaveLengthRange( const TFloat64Range& waveLengthRange ) const;

    Bool                HasUniformResolution() const;

    Bool                ConvertToLinearScale();
    Bool                ConvertToLogScale();
    Bool                IsInLogScale() const;
    Bool                IsInLinearScale() const;

    TLambdaRange        GetLambdaRange() const;
    Bool                ClampLambdaRange( const TFloat64Range& range, TFloat64Range& clampedRange ) const;
    Void                GetMask( const TFloat64Range& range,  CMask& mask ) const;

    Bool                PlotResolution( const char* filePath ) const;

private:

    UInt32              m_SpectralFlags;
};

}

#endif
