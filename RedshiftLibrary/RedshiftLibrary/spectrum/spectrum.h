#ifndef _REDSHIFT_SPECTRUM_SPECTRUM_
#define _REDSHIFT_SPECTRUM_SPECTRUM_

#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/spectrum/fluxaxis.h>
#include <RedshiftLibrary/spectrum/spectralaxis.h>
#include <RedshiftLibrary/continuum/continuum.h>

#include <string>

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class CSpectrum
{

public:

    enum EFLags
    {

    };

    CSpectrum();
    CSpectrum(const CSpectrum& other, TFloat64List mask);
    ~CSpectrum();

    CSpectrum& operator=(const CSpectrum& other);

    Void  SetName( const char* name );

    Bool InvertFlux();

    const CSpectrumSpectralAxis&    GetSpectralAxis() const;
    const CSpectrumFluxAxis&        GetFluxAxis() const;

    const std::string               GetName() const;

    CSpectrumFluxAxis&              GetFluxAxis();
    CSpectrumSpectralAxis&          GetSpectralAxis();

    UInt32                          GetSampleCount() const;
    Float64                         GetResolution() const;
    Float64                         GetMeanResolution() const;
    TLambdaRange                    GetLambdaRange() const;

    bool                            GetMeanAndStdFluxInRange( TFloat64Range wlRange, Float64& mean, Float64& std ) const;
    bool                            GetLinearRegInRange( TFloat64Range wlRange,  Float64& a, Float64& b) const;

    Bool                            ConvertToLogScale();
    Bool                            ConvertToLinearScale();

    Bool                            RemoveContinuum( CContinuum& remover );
    const Bool                      IsFluxValid(Float64 LambdaMin, Float64 LambdaMax) const;
    const Bool                      IsNoiseValid(Float64 LambdaMin, Float64 LambdaMax) const;

    const std::string&       	    GetFullPath() const;
    const Int32                     GetDecompScales() const;
    void 			    SetFullPath(const char* nameP);
    void 			    SetDecompScales(Int32 decompScales);

private:

    std::string                     m_Name;
    CSpectrumFluxAxis               m_FluxAxis;
    CSpectrumSpectralAxis           m_SpectralAxis;
    std::string                     m_FullPath;
    Int32                           m_nbScales;
};

inline
UInt32 CSpectrum::GetSampleCount() const
{
    return m_SpectralAxis.GetSamplesCount();
}

inline
const CSpectrumSpectralAxis& CSpectrum::GetSpectralAxis() const
{
    return m_SpectralAxis;
}

inline
const CSpectrumFluxAxis& CSpectrum::GetFluxAxis() const
{
    return m_FluxAxis;
}

inline
CSpectrumSpectralAxis& CSpectrum::GetSpectralAxis()
{
    return m_SpectralAxis;
}

inline
CSpectrumFluxAxis& CSpectrum::GetFluxAxis()
{
    return m_FluxAxis;
}


}

#endif
