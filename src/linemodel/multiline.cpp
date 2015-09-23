
#include <epic/redshift/linemodel/element.h>
#include <epic/redshift/linemodel/multiline.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/core/log/log.h>

#include <algorithm>

using namespace NSEpic;

CMultiLine::CMultiLine(std::vector<CRay> rs, std::vector<Float64> nominalAmplitudes, Float64 nominalWidth, std::vector<Int32> catalogIndexes)
{
    m_Rays = rs;
    m_NominalWidth = nominalWidth;
    m_NominalAmplitudes = nominalAmplitudes;

    m_NSigmaSupport = 8.0;

    for(int i=0; i<catalogIndexes.size(); i++){
        m_LineCatalogIndexes.push_back(catalogIndexes[i]);
    }
}

CMultiLine::~CMultiLine()
{
}


void CMultiLine::prepareSupport(const CSpectrumSpectralAxis& spectralAxis, Float64 redshift)
{
    m_OutsideLambdaRange=true;
    m_Start.resize(m_Rays.size());
    m_End.resize(m_Rays.size());
    m_OutsideLambdaRangeList.resize(m_Rays.size());
    for(Int32 i=0; i<m_Rays.size(); i++){
        Float64 mu = m_Rays[i].GetPosition()*(1+redshift);
        Float64 c = m_NominalWidth*(1.0+redshift);
        Float64 winsize = m_NSigmaSupport*c;
        m_Start[i] = spectralAxis.GetIndexAtWaveLength(mu-winsize/2.0);
        m_End[i] = spectralAxis.GetIndexAtWaveLength(mu+winsize/2.0);

        if(m_Start[i] <= 0 || m_End[i] >= spectralAxis.GetSamplesCount()-1 || m_End[i] <=0 || m_Start[i] >=spectralAxis.GetSamplesCount()-1 ){
            m_OutsideLambdaRangeList[i]=true;
        }else{
            m_OutsideLambdaRangeList[i]=false;
        }

        // set the global outside lambda range
        m_OutsideLambdaRange = m_OutsideLambdaRange && m_OutsideLambdaRangeList[i];
    }

    for(Int32 i=0; i<m_Rays.size(); i++){
        if(m_OutsideLambdaRangeList[i]){
            continue;
        }
        for(Int32 j=0; j<m_Rays.size(); j++){
            if(m_OutsideLambdaRangeList[j]){
                continue;
            }
            if(i==j){
                continue;
            }
            if(m_Rays[i].GetPosition()<m_Rays[j].GetPosition() && m_Start[j]<=m_End[i]){
                m_Start[j]=m_End[i]+1;
            }
        }
    }
}


Float64 CMultiLine::GetFittedAmplitude(Int32 subeIdx){
    return m_FittedAmplitudes[subeIdx];
}


void CMultiLine::fitAmplitude(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64  redshift)
{
    prepareSupport(spectralAxis, redshift);

    m_FittedAmplitudes.resize(m_Rays.size());
    for(Int32 k=0; k<m_Rays.size(); k++){
        m_FittedAmplitudes[k] = -1.0;
    }

    if(m_OutsideLambdaRange){
        return;
    }

    const Float64* flux = fluxAxis.GetSamples();
    const Float64* spectral = spectralAxis.GetSamples();
    const Float64* error = fluxAxis.GetError();

    Float64 y = 0.0;
    Float64 x = 0.0;
    Float64 yg = 0.0;

    Float64 sumCross = 0.0;
    Float64 sumGauss = 0.0;
    Float64 err2 = 0.0;
    Int32 num = 0;
    for(Int32 k=0; k<m_Rays.size(); k++){ //loop for the intervals
        if(m_OutsideLambdaRangeList[k]){
            continue;
        }

        //A estimation
        for ( Int32 i = m_Start[k]; i <= m_End[k]; i++)
        {
            y = flux[i];
            x = spectral[i];

            yg = 0.0;
            for(Int32 k2=0; k2<m_Rays.size(); k2++){ //loop for the signal synthesis
                Float64 mu = m_Rays[k2].GetPosition()*(1+redshift);
                Float64 c = m_NominalWidth*(1.0+redshift);
                yg += m_NominalAmplitudes[k2] * exp (-1.*(x-mu)*(x-mu)/(2*c*c));
            }
            num++;
            err2 = 1.0 / (error[i] * error[i]);
            sumCross += yg*y*err2;
            sumGauss += yg*yg*err2;
        }

    }
    if ( num==0 || sumCross==0 || sumGauss==0 )
    {
        return;
    }

    Float64 A = std::max(0.0, sumCross / sumGauss);

    for(Int32 k=0; k<m_Rays.size(); k++){
        m_FittedAmplitudes[k] = A*m_NominalAmplitudes[k];
    }
    //

    return;
}


void CMultiLine::addToSpectrumModel( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift )
{
    if(m_OutsideLambdaRange){
        return;
    }

    Float64 c = m_NominalWidth*(1.0+redshift);
    Float64* flux = modelfluxAxis.GetSamples();
    const Float64* spectral = modelspectralAxis.GetSamples();
    for(Int32 k=0; k<m_Rays.size(); k++){ //loop on the interval
        if(m_OutsideLambdaRangeList[k]){
            continue;
        }

        for ( Int32 i = m_Start[k]; i <= m_End[k]; i++)
        {
            Float64 x = spectral[i];

            Float64 Yi=0.0;
            for(Int32 k2=0; k2<m_Rays.size(); k2++) //loop
            {
                Float64 A = m_FittedAmplitudes[k2];
                Float64 mu = m_Rays[k2].GetPosition()*(1+redshift);
                Yi += A * exp (-1.*(x-mu)*(x-mu)/(2*c*c));
            }
            flux[i] += Yi;
        }
    }
  return;
}

Int32 CMultiLine::FindElementIndex(std::string LineTagStr)
{
    Int32 idx = -1;
    for( UInt32 iElts=0; iElts<m_Rays.size(); iElts++ )
    {
        std::string name = m_Rays[iElts].GetName();
        std::size_t foundstra = name.find(LineTagStr.c_str());

        if (foundstra!=std::string::npos){
            idx = iElts;
            break;
        }
    }

    return idx;
}

void CMultiLine::LimitFittedAmplitude(Int32 subeIdx, Float64 limit){
    //..
    return;
}
