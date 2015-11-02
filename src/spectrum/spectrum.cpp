#include <epic/redshift/spectrum/spectrum.h>

#include <epic/core/debug/assert.h>

#include <math.h>
#include <stdio.h>
#include <algorithm>

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CSpectrum )

CSpectrum::CSpectrum()
{

}


CSpectrum::CSpectrum(const CSpectrum& other, TFloat64List mask)
{
    TFloat64List tmpFlux;
    TFloat64List tmpError;
    TFloat64List tmpWave;

    CSpectrumSpectralAxis otherSpectral = other.GetSpectralAxis();
    CSpectrumFluxAxis otherFlux = other.GetFluxAxis();


    const Float64* error = otherFlux.GetError();

    for(Int32 i=0; i<mask.size(); i++){
        if(mask[i]!=0){
            tmpWave.push_back(otherSpectral[i]);
            tmpFlux.push_back(otherFlux[i]);
            if( error!= NULL ){
                tmpError.push_back(error[i]);
            }
        }
    }

    CSpectrumSpectralAxis *_SpectralAxis = new CSpectrumSpectralAxis(tmpWave.size(), otherSpectral.IsInLogScale());
    CSpectrumFluxAxis *_FluxAxis = new CSpectrumFluxAxis(tmpFlux.size());

    for(Int32 i=0; i<tmpFlux.size(); i++){
        (*_SpectralAxis)[i] = tmpWave[i];
        (*_FluxAxis)[i] = tmpFlux[i];
        if( error!= NULL ){
            (*_FluxAxis).GetError()[i] = tmpError[i];
        }
    }

    m_SpectralAxis = *_SpectralAxis;
    m_FluxAxis = *_FluxAxis;
}

CSpectrum::~CSpectrum()
{

}

CSpectrum& CSpectrum::operator=(const CSpectrum& other)
{
    m_SpectralAxis = other.GetSpectralAxis();
    m_FluxAxis = other.GetFluxAxis();
}

/**
 * Invert the flux axis
 */
Bool CSpectrum::InvertFlux()
{
    m_FluxAxis.Invert();
}

/**
 * Convert the spectral axis to a neperian logarithm scale
 */
Bool CSpectrum::ConvertToLogScale()
{
    return m_SpectralAxis.ConvertToLogScale();
}

/**
 * Convert the spectral axis to a linear scale
 */
Bool CSpectrum::ConvertToLinearScale()
{
    return m_SpectralAxis.ConvertToLinearScale();
}

Float64 CSpectrum::GetResolution() const
{
    return m_SpectralAxis.GetResolution();
}

Float64 CSpectrum::GetMeanResolution() const
{
    return m_SpectralAxis.GetMeanResolution();
}

/**
 * Return the lambda range of the entire spectrum.
 * Range is always expressed in linear scale NOT in log scale even if the underlying spcetrum is in log scale
 */
TLambdaRange CSpectrum::GetLambdaRange() const
{
    return m_SpectralAxis.GetLambdaRange();
}


const std::string CSpectrum::GetName() const
{
    return m_Name;
}

Void CSpectrum::SetName( const char* name )
{
    m_Name = name;
}
