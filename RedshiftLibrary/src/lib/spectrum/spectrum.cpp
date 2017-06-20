#include <RedshiftLibrary/spectrum/spectrum.h>


#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/debug/assert.h>

#include <math.h>
#include <stdio.h>
#include <algorithm>

#include <gsl/gsl_fit.h>

using namespace NSEpic;

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
    return *this;
}


Bool CSpectrum::RemoveContinuum(  CContinuum& remover )
{
    CSpectrumFluxAxis fluxAxisWithoutContinuum;

    remover.RemoveContinuum( *this, fluxAxisWithoutContinuum );

    m_FluxAxis = fluxAxisWithoutContinuum;

    return true;
}

/**
 * Invert the flux axis
 */
Bool CSpectrum::InvertFlux()
{
    m_FluxAxis.Invert();
    return true;
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

bool CSpectrum::GetMeanAndStdFluxInRange(TFloat64Range wlRange,  Float64& mean, Float64 &std) const
{
    //wlrange should be totally included in the spectrum lambdarange
    if(wlRange.GetBegin()<m_SpectralAxis.GetLambdaRange().GetBegin())
    {
        return false;
    }
    if(wlRange.GetEnd()>m_SpectralAxis.GetLambdaRange().GetEnd())
    {
        return false;
    }

    CMask mask;
    m_SpectralAxis.GetMask( wlRange, mask );
    const Float64* error = m_FluxAxis.GetError();
    Float64 _Mean = 0.0;
    Float64 _SDev = 0.0;
    m_FluxAxis.ComputeMeanAndSDev( mask,_Mean ,_SDev, error);

    mean = _Mean;
    std = _SDev;
    return true;
}

bool CSpectrum::GetLinearRegInRange(TFloat64Range wlRange,  Float64& a, Float64 &b) const
{
    //wlrange should be totally included in the spectrum lambdarange
    if(wlRange.GetBegin()<m_SpectralAxis.GetLambdaRange().GetBegin())
    {
        return false;
    }
    if(wlRange.GetEnd()>m_SpectralAxis.GetLambdaRange().GetEnd())
    {
        return false;
    }

    const Float64* error = m_FluxAxis.GetError();

    TInt32Range iRange = m_SpectralAxis.GetIndexesAtWaveLengthRange(wlRange);
    Int32 n = iRange.GetEnd()-iRange.GetBegin()+1;
    Float64* x = (Float64*)malloc(n*sizeof(Float64));
    Float64* y = (Float64*)malloc(n*sizeof(Float64));
    Float64* w = (Float64*)malloc(n*sizeof(Float64));

    for(Int32 k=0; k<n; k++)
    {
        Int32 ik=k+iRange.GetBegin();
        w[k] = 1.0/(error[ik]*error[ik]);
        x[k] = m_SpectralAxis[ik];
        y[k] = m_FluxAxis[ik];
    }

    double c0, c1, cov00, cov01, cov11, chisq;
    gsl_fit_wlinear (x, 1, w, 1, y, 1, n, &c0, &c1, &cov00, &cov01, &cov11, &chisq);

    a = c1;
    b = c0;
    //todo: use chisq (weighted sum of squares of the residuals) for fitting quality,
    free(x);
    free(y);
    free(w);
    return true;
}

const std::string CSpectrum::GetName() const
{
    return m_Name;
}

Void CSpectrum::SetName( const char* name )
{
    m_Name = name;
}

const Bool CSpectrum::IsFluxValid( Float64 LambdaMin,  Float64 LambdaMax ) const
{
    Bool allzero=true;

    const Float64* flux = m_FluxAxis.GetSamples();
    Int32 iMin = m_SpectralAxis.GetIndexAtWaveLength(LambdaMin);
    Int32 iMax = m_SpectralAxis.GetIndexAtWaveLength(LambdaMax);
    for(Int32 i=iMin; i<iMax; i++){
        if( flux[i] != 0.0 ){
            allzero = false;
        }
    }
    Bool valid = !allzero;
    return valid;
}

const Bool CSpectrum::IsNoiseValid( Float64 LambdaMin,  Float64 LambdaMax ) const
{
    Bool valid=true;

    const Float64* error = m_FluxAxis.GetError();
    Int32 iMin = m_SpectralAxis.GetIndexAtWaveLength(LambdaMin);
    Int32 iMax = m_SpectralAxis.GetIndexAtWaveLength(LambdaMax);
    for(Int32 i=iMin; i<iMax; i++){
        if( error[i] <= 0 ){
            valid = false;
        }
    }
    return valid;
}


const std::string& CSpectrum::GetFullPath() const
{
      return m_FullPath;
}

const Int32 CSpectrum::GetDecompScales() const
{
      return m_nbScales;
}

void CSpectrum::SetFullPath(const char* nameP)
{
	m_FullPath = nameP;
}

void CSpectrum::SetDecompScales( Int32 decompScales )
{
	m_nbScales = decompScales;
}
