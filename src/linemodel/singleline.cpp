
#include <epic/redshift/linemodel/element.h>
#include <epic/redshift/linemodel/singleline.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/core/log/log.h>

#include <algorithm>

using namespace NSEpic;

CSingleLine::CSingleLine(const CRay& r , Float64 nominalWidth, std::vector<Int32> catalogIndexes)
{
    m_Ray = r;
    m_NominalWidth = nominalWidth;
    m_FittedAmplitude = -1;

    m_NSigmaSupport = 8.0;
    m_Start = -1;
    m_End = -1;

    for(int i=0; i<catalogIndexes.size(); i++){
        m_LineCatalogIndexes.push_back(catalogIndexes[i]);
    }
}

CSingleLine::~CSingleLine()
{
}

void CSingleLine::fitAmplitude(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64  redshift)
{
    prepareSupport(spectralAxis, redshift);

    if(m_OutsideLambdaRange){
        return;
    }


    const Float64* flux = fluxAxis.GetSamples();
    const Float64* spectral = spectralAxis.GetSamples();
    const Float64* error = fluxAxis.GetError();
    Float64 mu = m_Ray.GetPosition()*(1+redshift);
    Float64 c = m_NominalWidth*(1.0+redshift);

    Float64 y = 0.0;
    Float64 x = 0.0;
    Float64 yg = 0.0;

    Float64 sumCross = 0.0;
    Float64 sumGauss = 0.0;
    Float64 err2 = 0.0;
    Int32 num = 0;

    //A estimation
    for ( Int32 i = m_Start; i < m_End; i++)
    {
        y = flux[i];
        x = spectral[i];
        yg = exp (-1.*(x-mu)*(x-mu)/(2*c*c));

        num++;
        err2 = 1.0 / (error[i] * error[i]);
        sumCross += yg*y*err2;
        sumGauss += yg*yg*err2;
    }

    if ( num==0 || sumCross==0 || sumGauss==0 )
    {
        return;
    }

    m_FittedAmplitude = std::max(0.0, sumCross / sumGauss);

    /*
    //SNR estimation
    Float64 sumErr  = 0.0;
    Float64 sumGaussA = 0.0;
    for ( Int32 i = start; i < end; i++)
    {
        x = spectral[i];
        sumGaussA += A*exp (-1.*(x-mu)*(x-mu)/(2*c*c));
        sumErr += (error[i] * error[i]);
    }
    Float64 SNRThres = 0.0001;
    if(sumGaussA/sumErr < SNRThres){
        A = 0.0;
    }
    */

    return;

}

void CSingleLine::addToSpectrumModel( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift )
{
    if(m_OutsideLambdaRange){
        return;
    }

    Float64 A = m_FittedAmplitude;
    Float64 mu = m_Ray.GetPosition()*(1+redshift);
    Float64 c = m_NominalWidth*(1.0+redshift);
    Float64* flux = modelfluxAxis.GetSamples();
    const Float64* spectral = modelspectralAxis.GetSamples();

    for ( Int32 i = m_Start; i < m_End; i++)
    {
        Float64 x = spectral[i];
        Float64 Yi = A * exp (-1.*(x-mu)*(x-mu)/(2*c*c));
        flux[i] += Yi;
    }

  return;
}

void CSingleLine::prepareSupport(const CSpectrumSpectralAxis& spectralAxis, Float64 redshift)
{
    Float64 mu = m_Ray.GetPosition()*(1+redshift);
    Float64 c = m_NominalWidth*(1.0+redshift);
    Float64 winsize = m_NSigmaSupport*c;
    m_Start = spectralAxis.GetIndexAtWaveLength(mu-winsize/2.0);
    m_End = spectralAxis.GetIndexAtWaveLength(mu+winsize/2.0);

    if(m_Start <= 0 || m_End >= spectralAxis.GetSamplesCount()-1){
        m_OutsideLambdaRange=true;
    }else{
        m_OutsideLambdaRange=false;
    }
}

 Float64 CSingleLine::GetFittedAmplitude(Int32 subeIdx){
     return m_FittedAmplitude;
 }
