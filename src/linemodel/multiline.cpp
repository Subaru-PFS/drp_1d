
#include <epic/redshift/linemodel/element.h>
#include <epic/redshift/linemodel/multiline.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/core/log/log.h>

#include <algorithm>

using namespace NSEpic;

CMultiLine::CMultiLine(std::vector<CRay> rs, Int32 widthType, std::vector<Float64> nominalAmplitudes, Float64 nominalWidth, std::vector<Int32> catalogIndexes):CLineModelElement(widthType)
{

    m_ElementType = "CMultiLine";
    m_Rays = rs;

    m_SignFactors.resize(m_Rays.size());
    for(Int32 i=0; i<m_Rays.size(); i++){
        if( m_Rays[i].GetType()==CRay::nType_Emission ){
            m_SignFactors[i] = 1.0;
        }else{
            m_SignFactors[i] = -1.0;
        }
    }
    m_NominalWidth = nominalWidth;
    m_NominalAmplitudes = nominalAmplitudes;

    m_NSigmaSupport = 8.0;

    for(int i=0; i<catalogIndexes.size(); i++){
        m_LineCatalogIndexes.push_back(catalogIndexes[i]);
    }

    SetFittedAmplitude(-1, -1);
}

CMultiLine::~CMultiLine()
{
}

std::string CMultiLine::GetRayName(Int32 subeIdx)
{
    if(subeIdx>=m_Rays.size())
    {
        return "-1";
    }

    return m_Rays[subeIdx].GetName();
}

Float64 CMultiLine::GetLineWidth(Float64 lambda, Float64 z)
{
    Float64 instrumentSigma = -1;
    if( m_LineWidthType == nWidthType_PSFInstrumentDriven){
        instrumentSigma = lambda/m_Resolution/m_FWHM_factor;
    }else if( m_LineWidthType == nWidthType_ZDriven){
        instrumentSigma = m_NominalWidth*(1+z);
    }else if( m_LineWidthType == nWidthType_Fixed){
        instrumentSigma = m_NominalWidth;
    }

//    Float64 v = 200;
//    Float64 c = 300000.0;
//    Float64 velocitySigma = v/c*lambda;//, useless /(1+z)*(1+z);
//    Float64 sigma = sqrt(instrumentSigma*instrumentSigma + velocitySigma*velocitySigma);

    return instrumentSigma;
}

void CMultiLine::prepareSupport(const CSpectrumSpectralAxis& spectralAxis, Float64 redshift, const TFloat64Range &lambdaRange)
{    
    m_OutsideLambdaRange=true;
    m_Start.resize(m_Rays.size());
    m_End.resize(m_Rays.size());
    m_OutsideLambdaRangeList.resize(m_Rays.size());
    for(Int32 i=0; i<m_Rays.size(); i++){
        Float64 mu = m_Rays[i].GetPosition()*(1+redshift);
        Float64 c = GetLineWidth(mu, redshift);
        Float64 winsize = m_NSigmaSupport*c;

        Float64 lambda_start = mu-winsize/2.0;
        if(lambda_start < lambdaRange.GetBegin()){
            lambda_start = lambdaRange.GetBegin();
        }
        m_Start[i] = spectralAxis.GetIndexAtWaveLength(lambda_start);

        Float64 lambda_end = mu+winsize/2.0;
        if(lambda_end > lambdaRange.GetEnd()){
            lambda_end = lambdaRange.GetEnd();
        }
        m_End[i] = spectralAxis.GetIndexAtWaveLength(lambda_end);

        Int32 minLineOverlap = m_OutsideLambdaRangeOverlapThreshold*winsize;
        if( m_Start[i] >= (spectralAxis.GetSamplesCount()-1-minLineOverlap) || m_End[i] <=minLineOverlap){
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

TInt32RangeList CMultiLine::getSupport()
{
    TInt32RangeList support;
    if(m_OutsideLambdaRange==false){
        for(Int32 i=0; i<m_Rays.size(); i++){
            if(m_OutsideLambdaRangeList[i]){
                continue;
            }

            support.push_back(TInt32Range(m_Start[i], m_End[i]));
        }

    }
    return support;
}

Float64 CMultiLine::GetFittedAmplitude(Int32 subeIdx){
    return m_FittedAmplitudes[subeIdx];
}

Float64 CMultiLine::GetFittedAmplitudeErrorSigma(Int32 subeIdx){
    return m_FittedAmplitudeErrorSigmas[subeIdx];
}

Float64 CMultiLine::GetElementAmplitude(){
    if(m_OutsideLambdaRange){

        return -1;
    }else{

        return m_FittedAmplitudes[0]/m_NominalAmplitudes[0];
    }
}

Float64 CMultiLine::GetNominalAmplitude(Int32 subeIdx){
    return m_NominalAmplitudes[subeIdx];
}

void CMultiLine::SetFittedAmplitude(Float64 A, Float64 SNR)
{
    m_FittedAmplitudes.resize(m_Rays.size());
    m_FittedAmplitudeErrorSigmas.resize(m_Rays.size());
    if(m_OutsideLambdaRange){
        for(Int32 k=0; k<m_Rays.size(); k++){
            m_FittedAmplitudes[k] = -1;
            m_FittedAmplitudeErrorSigmas[k] = -1;
        }
        return;
    }else{
        A = std::max(0.0, A);
        for(Int32 k=0; k<m_Rays.size(); k++){
            if(m_OutsideLambdaRangeList[k]){
                m_FittedAmplitudes[k] = -1;
            }
            m_FittedAmplitudes[k] = A*m_NominalAmplitudes[k];
            m_FittedAmplitudeErrorSigmas[k] = SNR*m_NominalAmplitudes[k]; //todo: check correct formulation for Error
        }
    }

}


void CMultiLine::fitAmplitude(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64  redshift)
{
    m_FittedAmplitudes.resize(m_Rays.size());
    m_FittedAmplitudeErrorSigmas.resize(m_Rays.size());
    for(Int32 k=0; k<m_Rays.size(); k++){
        m_FittedAmplitudes[k] = -1.0;
        m_FittedAmplitudeErrorSigmas[k] = -1.0;
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

    std::vector<Float64> mu;
    mu.resize(m_Rays.size());
    std::vector<Float64> c;
    c.resize(m_Rays.size());
    for(Int32 k2=0; k2<m_Rays.size(); k2++){ //loop for the signal synthesis
        mu[k2] = m_Rays[k2].GetPosition()*(1+redshift);
        c[k2] = GetLineWidth(mu[k2], redshift);
    }

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
                if(m_OutsideLambdaRangeList[k2]){
                    continue;
                }
                yg += m_SignFactors[k2] * m_NominalAmplitudes[k2] * exp (-1.*(x-mu[k2])*(x-mu[k2])/(2*c[k2]*c[k2]));
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
        if(m_OutsideLambdaRangeList[k]){
            continue;
        }
        m_FittedAmplitudes[k] = A*m_NominalAmplitudes[k];
        if(A==0){
            m_FittedAmplitudeErrorSigmas[k] = 0.0;
        }else{
            m_FittedAmplitudeErrorSigmas[k] = 1.0/sqrt(sumGauss);

        }
    }
    //

    return;
}


void CMultiLine::addToSpectrumModel( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift )
{
    if(m_OutsideLambdaRange){
        return;
    }


    const Float64* spectral = modelspectralAxis.GetSamples();
    Float64* flux = modelfluxAxis.GetSamples();
    for(Int32 k=0; k<m_Rays.size(); k++){ //loop on the interval
        if(m_OutsideLambdaRangeList[k]){
            continue;
        }

        for ( Int32 i = m_Start[k]; i <= m_End[k]; i++)
        {
            Float64 lambda = spectral[i];
            Float64 Yi=getModelAtLambda(lambda, redshift);
            flux[i] += Yi;
        }
    }
  return;
}

Float64 CMultiLine::getModelAtLambda(Float64 lambda, Float64 redshift )
{
    if(m_OutsideLambdaRange){
        return 0.0;
    }
    Float64 Yi=0.0;

    Float64 x = lambda;

    for(Int32 k2=0; k2<m_Rays.size(); k2++) //loop on rays
    {
        if(m_OutsideLambdaRangeList[k2]){
            continue;
        }
        Float64 A = m_FittedAmplitudes[k2];
        Float64 mu = m_Rays[k2].GetPosition()*(1+redshift);
        Float64 c = GetLineWidth(mu, redshift);
        Yi += m_SignFactors[k2] * A * exp (-1.*(x-mu)*(x-mu)/(2*c*c));
    }
    return Yi;
}

void CMultiLine::initSpectrumModel( CSpectrumFluxAxis &modelfluxAxis, CSpectrumFluxAxis &continuumfluxAxis )
{
    if(m_OutsideLambdaRange){
        return;
    }

    Float64* flux = modelfluxAxis.GetSamples();
    for(Int32 k=0; k<m_Rays.size(); k++){ //loop on the interval
        if(m_OutsideLambdaRangeList[k]){
            continue;
        }

        for ( Int32 i = m_Start[k]; i <= m_End[k]; i++)
        {
            flux[i] = continuumfluxAxis[i];
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

    if(m_FittedAmplitudes[subeIdx] > limit){
        m_FittedAmplitudes[subeIdx] = std::max(0.0, limit);
    }
    return;
}


bool CMultiLine::IsOutsideLambdaRange(Int32 subeIdx){
    return m_OutsideLambdaRangeList[subeIdx];
}
