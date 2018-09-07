#ifndef _REDSHIFT_SPECTRUM_FLUXAXIS_
#define _REDSHIFT_SPECTRUM_FLUXAXIS_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/common/range.h>

namespace NSEpic
{

class CMask;
class CSpectrum;
class CSpectrumSpectralAxis;

/**
 * \ingroup Redshift
 */
class CSpectrumFluxAxis : public CSpectrumAxis
{

public:

    CSpectrumFluxAxis();
    CSpectrumFluxAxis( UInt32 n );
    CSpectrumFluxAxis( const Float64* samples, UInt32 n );
    CSpectrumFluxAxis( const Float64* samples, UInt32 n, const Float64* error, UInt32 m );
    ~CSpectrumFluxAxis();

    CSpectrumFluxAxis& operator=(const CSpectrumFluxAxis& other);

    const TFloat64List&      GetError() const;
    TFloat64List&            GetError();

    void                SetSize( UInt32 s );

    Bool                ApplyMeanSmooth( UInt32 kernelHalfWidth );
    Bool                ApplyMedianSmooth( UInt32 kernelHalfWidth );


    Bool                ComputeMeanAndSDev( const CMask& mask, Float64& mean,  Float64& sdev, const TFloat64List error ) const;
    Float64             ComputeRMSDiff( const CSpectrumFluxAxis& other );
    Bool                Subtract(const CSpectrumFluxAxis& other);
    Bool                Invert();

    static Bool         Rebin( const TFloat64Range& range, const CSpectrumFluxAxis& sourceFluxAxis, const CSpectrumSpectralAxis& sourceSpectralAxis, const CSpectrumSpectralAxis& targetSpectralAxis,
                               CSpectrumFluxAxis& rebinedFluxAxis, CSpectrumSpectralAxis& rebinedSpectralAxis, CMask& rebinedMask );
    static Bool         Rebin2(const TFloat64Range& range, const CSpectrumFluxAxis& sourceFluxAxis, const Float64 *pfgTplBuffer, Float64 sourcez, const CSpectrumSpectralAxis& sourceSpectralAxis, const CSpectrumSpectralAxis& targetSpectralAxis,
                               CSpectrumFluxAxis& rebinedFluxAxis, CSpectrumSpectralAxis& rebinedSpectralAxis, CMask& rebinedMask , const std::string opt_interp);
    static Bool         RebinVarianceWeighted( const CSpectrumFluxAxis& sourceFluxAxis, const CSpectrumSpectralAxis& sourceSpectralAxis, const CSpectrumFluxAxis& sourceError,
                                                   const CSpectrumSpectralAxis& targetSpectralAxis,
                                                   CSpectrumFluxAxis& rebinedFluxAxis, CSpectrumSpectralAxis& rebinedSpectralAxis, CSpectrumFluxAxis& rebinedError,
                                                   const std::string opt_interp );

private:

    Bool                ComputeMeanAndSDevWithoutError( const CMask& mask, Float64& mean,  Float64& sdev) const;
    Bool                ComputeMeanAndSDevWithError( const CMask& mask, Float64& mean, Float64& sdev, const TFloat64List error ) const;

    TFloat64List        m_StatError;

};

inline
TFloat64List& CSpectrumFluxAxis::GetError()
{
  return m_StatError;
}

inline
const TFloat64List& CSpectrumFluxAxis::GetError() const
{
  return m_StatError;
}

}

#endif
