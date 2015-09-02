#ifndef _REDSHIFT_SPECTRUM_FLUXAXIS_
#define _REDSHIFT_SPECTRUM_FLUXAXIS_

#include <epic/redshift/common/datatypes.h>
#include <epic/redshift/spectrum/axis.h>
#include <epic/core/common/range.h>

namespace NSEpic
{

class CMask;
class CSpectrum;
class CSpectrumSpectralAxis;

class CSpectrumFluxAxis : public CSpectrumAxis
{

public:

    CSpectrumFluxAxis();
    CSpectrumFluxAxis( UInt32 n );
    CSpectrumFluxAxis( const Float64* samples, UInt32 n );
    ~CSpectrumFluxAxis();

    CSpectrumFluxAxis& operator=(const CSpectrumFluxAxis& other);

    const Float64*      GetError() const;
    Float64*            GetError();

    Void                SetSize( UInt32 s );

    Bool                ApplyMeanSmooth( UInt32 kernelHalfWidth );
    Bool                ApplyMedianSmooth( UInt32 kernelHalfWidth );


    Bool                ComputeMeanAndSDev( const CMask& mask, Float64& mean,  Float64& sdev, const Float64* error ) const;
    Float64             ComputeRMSDiff( const CSpectrumFluxAxis& other );
    Bool                Subtract(const CSpectrumFluxAxis& other);
    Bool                Invert();

    static Bool         Rebin( const TFloat64Range& range, const CSpectrumFluxAxis& sourceFluxAxis, const CSpectrumSpectralAxis& sourceSpectralAxis, const CSpectrumSpectralAxis& targetSpectralAxis,
                               CSpectrumFluxAxis& rebinedFluxAxis, CSpectrumSpectralAxis& rebinedSpectralAxis, CMask& rebinedMask );
    static Bool         Rebin2(const TFloat64Range& range, const CSpectrumFluxAxis& sourceFluxAxis, const Float64 *pfgTplBuffer, Float64 sourcez, const CSpectrumSpectralAxis& sourceSpectralAxis, const CSpectrumSpectralAxis& targetSpectralAxis,
                               CSpectrumFluxAxis& rebinedFluxAxis, CSpectrumSpectralAxis& rebinedSpectralAxis, CMask& rebinedMask );

private:

    Bool                ComputeMeanAndSDevWithoutError( const CMask& mask, Float64& mean,  Float64& sdev) const;
    Bool                ComputeMeanAndSDevWithError( const CMask& mask, Float64& mean, Float64& sdev, const Float64* error ) const;

    TFloat64List        m_StatError;

};

inline
Float64* CSpectrumFluxAxis::GetError()
{
    return m_StatError.data();
}

inline
const Float64* CSpectrumFluxAxis::GetError() const
{
    return m_StatError.data();
}

}

#endif
