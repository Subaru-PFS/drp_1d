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

    CSpectrumFluxAxis() = default;
    //value =0. is a default value for flux and not for error.
    explicit CSpectrumFluxAxis( UInt32 n, Float64 value = 0.0);
    CSpectrumFluxAxis( CSpectrumAxis otherFlux, CSpectrumNoiseAxis otherError );
    CSpectrumFluxAxis( const Float64* samples, UInt32 n );
    CSpectrumFluxAxis( const TFloat64List & samples);
    CSpectrumFluxAxis( TFloat64List && samples);
    CSpectrumFluxAxis( const Float64* samples, UInt32 n, const Float64* error, const UInt32 m);

    const CSpectrumNoiseAxis&      GetError() const;
    CSpectrumNoiseAxis&            GetError();

    void                SetSize( UInt32 s );
    void                clear();
    Bool                ApplyMeanSmooth( UInt32 kernelHalfWidth );
    Bool                ApplyMedianSmooth( UInt32 kernelHalfWidth );


    Bool                ComputeMeanAndSDev( const CMask& mask, Float64& mean,  Float64& sdev) const;
    Float64             ComputeRMSDiff( const CSpectrumFluxAxis& other );
    Bool                Subtract(const CSpectrumFluxAxis& other);
    Bool                Invert();
    Int32               extractFrom(const CSpectrumFluxAxis& other, Int32 startIdx, Int32 endIdx);//this is mainly applied on m_StdError

private:

    Bool                ComputeMeanAndSDevWithoutError( const CMask& mask, Float64& mean,  Float64& sdev) const;
    Bool                ComputeMeanAndSDevWithError( const CMask& mask, Float64& mean, Float64& sdev) const;

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

inline
Int32 CSpectrumFluxAxis::extractFrom(const CSpectrumFluxAxis& other, Int32 startIdx, Int32 endIdx)
{
    (*this).extractFrom(other, startIdx, endIdx); 
    m_StdError.SetSize(endIdx-startIdx +1);
    for(Int32 i = startIdx; i < endIdx + 1; i++){
        m_StdError[i - startIdx] = other[i];
    }
    return 0;
}
}

#endif
