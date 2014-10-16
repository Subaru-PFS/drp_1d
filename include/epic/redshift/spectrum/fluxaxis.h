#ifndef _REDSHIFT_SPECTRUM_FLUXAXIS_
#define _REDSHIFT_SPECTRUM_FLUXAXIS_

#include <epic/redshift/common/datatypes.h>
#include <epic/redshift/spectrum/axis.h>

namespace __NS__
{

class CMask;

class CSpectrumFluxAxis : public CSpectrumAxis
{

public:

    CSpectrumFluxAxis();
    CSpectrumFluxAxis( UInt32 n );
    CSpectrumFluxAxis( const Float64* samples, UInt32 n );
    ~CSpectrumFluxAxis();

    const Float64*      GetError() const;
    Float64*            GetError();

    Void                SetSize( UInt32 s );

    Bool                ApplyMeanSmooth( UInt32 kernelHalfWidth );
    Bool                ApplyMedianSmooth( UInt32 kernelHalfWidth );


    Bool                ComputeMeanAndSDev( const CMask& mask, Float64& mean,  Float64& sdev, const Float64* error ) const;

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
