
#include <epic/redshift/linemodel/element.h>
#include <epic/redshift/linemodel/singleline.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/core/log/log.h>

#include <float.h>
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
        m_FittedAmplitude = -1;
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

    if ( num==0 || sumGauss==0 )
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

//    Float64 SNRThres = 0.1;
//    if(m_FittedAmplitude/sqrt(sumGauss) < SNRThres){
//        m_FittedAmplitude = 0.0;
//    }


    return;

}


//Float64 CSingleLine::FitAmplitudeIterative( const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64 lambda, Float64 width, Int32 start, Int32 end)
//{
//    Float64 A = boost::numeric::bounds<float>::lowest();
//    const Float64* flux = fluxAxis.GetSamples();
//    const Float64* spectral = spectralAxis.GetSamples();
//    const Float64* error = fluxAxis.GetError();

//    //A first guess
//    for ( Int32 i = start; i < end; i++)
//    {
//        Float64 y = flux[i];
//        if(y>A){
//            A = y;
//        }
//    }

//    if(A<=0){
//        return 0.0;
//    }
//    //A fitting iteration loop
//    A = A*1.5;
//    Float64 mu = lambda;
//    Float64 c = width;
//    Float64 thres = 1e-5;
//    Int32 maxIteration = 100;
//    Float64 AstepDown = A/((Float64)(maxIteration+1));
//    Float64 sum2 = boost::numeric::bounds<float>::highest();
//    Float64 sum2prev = boost::numeric::bounds<float>::highest();
//    Int32 icmpt = 0;
//    while( sum2prev>=sum2 && sum2>thres && icmpt<maxIteration){
//        sum2prev = sum2;
//        sum2 = 0.0;
//        for ( Int32 i = start; i < end; i++)
//        {
//            Float64 x = spectral[i];
//            Float64 Yi = A * exp (-1.*(x-mu)*(x-mu)/(2*c*c));
//            //sum2 += Yi-flux[i];
//            sum2 += pow( Yi - flux[i] , 2.0 ) / pow( error[i], 2.0 );
//        }
//        //sum2 /= (Float64)(end-start+1);
//        icmpt++;
//        A = A-AstepDown;
//    }

//    if(A<0){
//        A=0;
//    }
//    return A;
//}

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

    if(m_Start <= 0 || m_End >= spectralAxis.GetSamplesCount()-1 || m_End <=0 || m_Start >=spectralAxis.GetSamplesCount()-1 ){
        m_OutsideLambdaRange=true;
    }else{
        m_OutsideLambdaRange=false;
    }
}

Float64 CSingleLine::GetFittedAmplitude(Int32 subeIdx){
//    if(m_OutsideLambdaRange){
//        m_FittedAmplitude = -1;
//    }

    return m_FittedAmplitude;
}

void CSingleLine::LimitFittedAmplitude(Int32 subeIdx, Float64 limit){
    if(m_FittedAmplitude > limit){
        m_FittedAmplitude = std::max(0.0, limit);
    }
}

Int32 CSingleLine::FindElementIndex(std::string LineTagStr)
{
    Int32 idx = -1;

    std::string name = m_Ray.GetName();
    std::size_t foundstra = name.find(LineTagStr.c_str());

    if (foundstra!=std::string::npos){
        idx = 0;
    }

    return idx;
}
