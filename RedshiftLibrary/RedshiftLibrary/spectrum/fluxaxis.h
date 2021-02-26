#ifndef _REDSHIFT_SPECTRUM_FLUXAXIS_
#define _REDSHIFT_SPECTRUM_FLUXAXIS_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/spectrum/noiseaxis.h>
#include <RedshiftLibrary/common/range.h>

namespace NSEpic
{

class CMask;

/**
 * \ingroup Redshift
 */
class CSpectrumFluxAxis : public CSpectrumAxis
{

public:

    CSpectrumFluxAxis();
    explicit CSpectrumFluxAxis( UInt32 n );
    CSpectrumFluxAxis( const Float64* samples, UInt32 n );
    CSpectrumFluxAxis( const Float64* samples, UInt32 n, const Float64* error);
    ~CSpectrumFluxAxis();

    CSpectrumFluxAxis& operator=(const CSpectrumFluxAxis& other);

    const CSpectrumNoiseAxis&      GetError() const;
    CSpectrumNoiseAxis&            GetError();

    void                SetSize( UInt32 s );

    Bool                ApplyMeanSmooth( UInt32 kernelHalfWidth );
    Bool                ApplyMedianSmooth( UInt32 kernelHalfWidth );


    Bool                ComputeMeanAndSDev( const CMask& mask, Float64& mean,  Float64& sdev, const CSpectrumNoiseAxis error ) const;
    Float64             ComputeRMSDiff( const CSpectrumFluxAxis& other );
    Bool                Subtract(const CSpectrumFluxAxis& other);
    Bool                Invert();

private:

    Bool                ComputeMeanAndSDevWithoutError( const CMask& mask, Float64& mean,  Float64& sdev) const;
    Bool                ComputeMeanAndSDevWithError( const CMask& mask, Float64& mean, Float64& sdev, const CSpectrumNoiseAxis error ) const;

    CSpectrumNoiseAxis        m_StdError;//STD

};

inline
CSpectrumNoiseAxis& CSpectrumFluxAxis::GetError()
{
  return m_StdError;
}

inline
const CSpectrumNoiseAxis& CSpectrumFluxAxis::GetError() const
{
  return m_StdError;
}

}

#endif
