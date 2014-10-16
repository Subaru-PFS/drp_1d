#ifndef _REDSHIFT_SPECTRUM_SPECTRUM_
#define _REDSHIFT_SPECTRUM_SPECTRUM_

#include <epic/core/common/managedobject.h>
#include <epic/core/common/range.h>
#include <epic/redshift/spectrum/fluxaxis.h>
#include <epic/redshift/spectrum/spectralaxis.h>
#include <epic/redshift/continuum/continuum.h>

namespace __NS__
{

class CSpectrum : public CManagedObject
{
    DEFINE_MANAGED_OBJECT( CSpectrum )

public:

    enum EFLags
    {

    };

    CSpectrum();
    ~CSpectrum();

    const CSpectrumSpectralAxis&    GetSpectralAxis() const;
    const CSpectrumFluxAxis&        GetFluxAxis() const;

    CSpectrumFluxAxis&              GetFluxAxis();
    CSpectrumSpectralAxis&          GetSpectralAxis();

    UInt32                          GetSampleCount() const;
    Float64                         GetResolution() const;
    TLambdaRange                    GetLambdaRange() const;

    Bool                            ConvertToLogScale();
    Bool                            ConvertToLinearScale();

    template< typename ContinuumRemover >
    Bool                            RemoveContinuum();

private:

    Bool                    Serialize( CSerializer& ar );

    CSpectrumFluxAxis               m_FluxAxis;
    CSpectrumSpectralAxis           m_SpectralAxis;
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

template<typename ContinuumRemover>
Bool CSpectrum::RemoveContinuum()
{
    ContinuumRemover cr;
    CSpectrumFluxAxis fluxAxisWithoutContinuum;

    cr.RemoveContinuum( *this, fluxAxisWithoutContinuum );

    m_FluxAxis = fluxAxisWithoutContinuum;

    return true;
}

}

#endif
