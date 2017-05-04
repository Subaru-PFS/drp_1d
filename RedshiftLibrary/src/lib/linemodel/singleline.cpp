
#include <RedshiftLibrary/linemodel/element.h>
#include <RedshiftLibrary/linemodel/singleline.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/log/log.h>

#include <float.h>
#include <algorithm>

using namespace NSEpic;

/**
 * Attribution constructor.
 */
CSingleLine::CSingleLine(const CRay& r , const std::string& widthType, const Float64 resolution, const Float64 velocityEmission, const Float64 velocityAbsorption, Float64 nominalWidth, std::vector<Int32> catalogIndexes):CLineModelElement(widthType, resolution, velocityEmission, velocityAbsorption)
{
    m_ElementType = "CSingleLine";
    m_Ray = r;

    if( m_Ray.GetType()==CRay::nType_Emission )
      {
        m_SignFactor = 1.0;
      }
    else
      {
        m_SignFactor = -1.0;
      }

    m_NominalWidth = nominalWidth;
    m_FittedAmplitude = -1;
    m_FittedAmplitudeErrorSigma = -1;

    m_Start = -1;
    m_End = -1;

    for( int i=0; i<catalogIndexes.size(); i++ )
      {
        m_LineCatalogIndexes.push_back( catalogIndexes[i] );
      }

    //initialize fitted amplitude
    SetFittedAmplitude(-1.0, -1.0 );
}

/**
 * Empty destructor.
 */
CSingleLine::~CSingleLine()
{

}

/**
 * Wrapper around m_Ray.GetName().
 */
std::string CSingleLine::GetRayName( Int32 subeIdx )
{
    return m_Ray.GetName();
}

/**
 * Gets m_SignFactor.
 */
Float64 CSingleLine::GetSignFactor( Int32 subeIdx )
{
    return m_SignFactor;
}

/**
 * Calculates and returns the line width for the input redshift.
 */
Float64 CSingleLine::GetWidth( Int32 subeIdx, Float64 redshift )
{
    Float64 mu = m_Ray.GetPosition()*(1+redshift);
    Float64 c = GetLineWidth( mu, redshift, m_Ray.GetIsEmission(), m_Ray.GetProfile() );

    return c;
}

/**
 * \brief Returns a vector containing m_Ray.
 **/
std::vector<CRay> CSingleLine::GetRays()
{
    std::vector<CRay> rays;
    rays.push_back(m_Ray);
    return rays;
}

/**
 * /brief Calculates the profile of a fit of this line with the specified redshift.
 * Sets some m_ variables depending on the input parameters.
 * Find the wavelength this line would have at the specified redshift.
 * Calculate the line width this line would have at the specified redshift.
 * This defines a gaussian profile, which is recorded by GetNSigmaSupport.
 * Check if this profile is within the specified lambdaRange. If it isn't, invalidate the current support.
 */
void CSingleLine::prepareSupport(const CSpectrumSpectralAxis& spectralAxis, Float64 redshift, const TFloat64Range &lambdaRange)
{
    Float64 mu = m_Ray.GetPosition()*(1+redshift);
    Float64 c = GetLineWidth(mu, redshift, m_Ray.GetIsEmission(), m_Ray.GetProfile());
    std::string profile = m_Ray.GetProfile();
    Float64 winsize = GetNSigmaSupport(profile)*c;
    Float64 lambda_start = mu-winsize/2.0;
    if( lambda_start<lambdaRange.GetBegin() )
      {
        lambda_start = lambdaRange.GetBegin();
      }
    m_Start = spectralAxis.GetIndexAtWaveLength( lambda_start );

    Float64 lambda_end = mu+winsize/2.0;
    if( lambda_end>lambdaRange.GetEnd() )
      {
        lambda_end = lambdaRange.GetEnd();
      }
    m_End = spectralAxis.GetIndexAtWaveLength( lambda_end );

    Int32 minLineOverlap = m_OutsideLambdaRangeOverlapThreshold*winsize;
    if( m_Start>=(spectralAxis.GetSamplesCount()-1-minLineOverlap) || m_End<=minLineOverlap )
      {
        m_OutsideLambdaRange=true;
      }
    else 
      {
        m_OutsideLambdaRange=false;
      }

    if( m_OutsideLambdaRange )
      {
        m_FittedAmplitude = -1.0;
        m_FittedAmplitudeErrorSigma = -1.0;
        return;
      }
}

/**
 * Returns an empty support if m_OutsideLambdaRange is false, or a TInt32Range ( m_Start, m_End ) otherwise.
 */
TInt32RangeList CSingleLine::getSupport()
{
    TInt32RangeList support;
    if( m_OutsideLambdaRange==false )
      {
        support.push_back( TInt32Range( m_Start, m_End ) );
      }
    return support;
}

/**
 * Returns an empty range if m_OtsideLambdaRange is false, and a range between m_Start and m_End otherwise.
 */
TInt32Range CSingleLine::getSupportSubElt( Int32 subeIdx )
{
    TInt32Range support;
     if( m_OutsideLambdaRange==false )
       {
         support = TInt32Range( m_Start, m_End );
       }
     return support;
}

/**
 * Calls getSupportSubElt()
 */
TInt32Range CSingleLine::getTheoreticalSupportSubElt( Int32 subeIdx )
{
     return getSupportSubElt(subeIdx);
}




/**
 * /brief Fits a gaussian considering the flux in a width calculated by the prepareSupport method.
 * Sets m_FittedAmplitude and m_FittedAmplitudeErrorSigma based on the parameters.
 */
void CSingleLine::fitAmplitude(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64  redshift, Int32 lineIdx )
{
    if( m_OutsideLambdaRange )
      {
        m_FittedAmplitude = -1.0;
        m_FittedAmplitudeErrorSigma = -1.0;
        return;
      }

    const Float64* flux = fluxAxis.GetSamples();
    const Float64* spectral = spectralAxis.GetSamples();
    const Float64* error = fluxAxis.GetError();
    Float64 mu = m_Ray.GetPosition()*(1+redshift);
    Float64 c = GetLineWidth(mu, redshift, m_Ray.GetIsEmission(), m_Ray.GetProfile());
    std::string profile = m_Ray.GetProfile();

    Float64 y = 0.0;
    Float64 x = 0.0;
    Float64 yg = 0.0;

    Float64 sumCross = 0.0;
    Float64 sumGauss = 0.0;
    Float64 err2 = 0.0;
    Int32 num = 0;

    //A estimation
    for ( Int32 i = m_Start; i <= m_End; i++ )
    {
        y = flux[i];
        x = spectral[i];
        yg = m_SignFactor * GetLineProfile(profile, x, mu, c);

        num++;
        err2 = 1.0 / (error[i] * error[i]);
        sumCross += yg*y*err2;
        sumGauss += yg*yg*err2;
    }

    if ( num==0 || sumGauss==0 )
    {
        return;
    }

    m_FittedAmplitude = sumCross / sumGauss;
    if( m_FittedAmplitude<0 )
      {
        m_FittedAmplitude = 0.0;
        m_FittedAmplitudeErrorSigma = 1.0/sqrt(sumGauss) - sumCross / sumGauss; //warning: sigma estimated = rough approx.
      }
    else
      {
        m_FittedAmplitudeErrorSigma = 1.0/sqrt(sumGauss);
//        //SNR estimation
//        Float64 SNRThres = 1.0;
//        if(m_FittedAmplitudeErrorSigma/m_FittedAmplitudeErrorSigma < SNRThres){
//            m_FittedAmplitudeErrorSigma = 0;
//            m_FittedAmplitude=0;
//        }
    }


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

/**
 * If applicable, add to each x the value from getModelAtLambda for the input redshift.
 */
void CSingleLine::addToSpectrumModel(const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift, Int32 lineIdx )
{
    if( m_OutsideLambdaRange )
    {
        return;
    }

    Float64 A = m_FittedAmplitude;
    if( A==0 )
    {
        return;
    }

    Float64* flux = modelfluxAxis.GetSamples();
    const Float64* spectral = modelspectralAxis.GetSamples();

    for ( Int32 i = m_Start; i<=m_End; i++ )
    {
        Float64 x = spectral[i];
        Float64 Yi = getModelAtLambda( x, redshift );
        flux[i] += Yi;
    }

    return;
}

/**
 * If applicable, add to each x the value from GetModelDerivSigmaAtLambda for the input redshift.
 */
void CSingleLine::addToSpectrumModelDerivSigma(const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift)
{
    if( m_OutsideLambdaRange )
    {
        return;
    }

    Float64 A = m_FittedAmplitude;
    if( A==0 )
    {
        return;
    }

    Float64* flux = modelfluxAxis.GetSamples();
    const Float64* spectral = modelspectralAxis.GetSamples();

    for ( Int32 i = m_Start; i<=m_End; i++ )
    {
        Float64 x = spectral[i];
        Float64 Yi = GetModelDerivSigmaAtLambda( x, redshift );
        flux[i] += Yi;
    }

    return;
}


/**
 * \brief Returns the amplitude for the argument lambda, at the argument redshift.
 * If outside lambda range, return 0.
 * Return the amplitude according to sign, fitted amplitude and line profile.
 */
Float64 CSingleLine::getModelAtLambda(Float64 lambda, Float64 redshift , Int32 kRaySupport)
{
    if(m_OutsideLambdaRange){
        return 0.0;
    }
    Float64 Yi=0.0;

    Float64 x = lambda;

    Float64 A = m_FittedAmplitude;
    Float64 mu = m_Ray.GetPosition()*(1+redshift);
    Float64 c = GetLineWidth(mu, redshift, m_Ray.GetIsEmission(), m_Ray.GetProfile());
    std::string profile = m_Ray.GetProfile();

    Yi = m_SignFactor * A * GetLineProfile(profile, x, mu, c);

    return Yi;
}

/**
 *
 */
Float64 CSingleLine::GetModelDerivAmplitudeAtLambda(Float64 lambda, Float64 redshift )
{
    if(m_OutsideLambdaRange){
        return 0.0;
    }
    Float64 Yi=0.0;

    Float64 x = lambda;

    Float64 mu = m_Ray.GetPosition()*(1+redshift);
    Float64 c = GetLineWidth(mu, redshift, m_Ray.GetIsEmission(), m_Ray.GetProfile());
    std::string profile = m_Ray.GetProfile();

    Yi = m_SignFactor * GetLineProfile(profile, x, mu, c);

    return Yi;
}

/**
 *
 */
Float64 CSingleLine::GetModelDerivSigmaAtLambda(Float64 lambda, Float64 redshift )
{
    if(m_OutsideLambdaRange){
        return 0.0;
    }
    Float64 Yi=0.0;

    Float64 x = lambda;

    Float64 A = m_FittedAmplitude;
    Float64 mu = m_Ray.GetPosition()*(1+redshift);
    Float64 c = GetLineWidth(mu, redshift, m_Ray.GetIsEmission(), m_Ray.GetProfile());
    std::string profile = m_Ray.GetProfile();

    Yi = m_SignFactor * A * GetLineProfileDerivSigma(profile, x, mu, c);

    return Yi;
}

/**
 * If outside lambda range, return.
 * Copy the continuum flux into the model flux.
 */
void CSingleLine::initSpectrumModel(CSpectrumFluxAxis& modelfluxAxis, CSpectrumFluxAxis &continuumfluxAxis, Int32 lineIdx)
{
    if( m_OutsideLambdaRange )
      {
        return;
      }
    Float64* flux = modelfluxAxis.GetSamples();
    for ( Int32 i=m_Start; i<=m_End; i++ )
    {
        flux[i] = continuumfluxAxis[i];
    }

  return;
}

/**
 * \brief Returns m_FittedAmplitude.
 */
Float64 CSingleLine::GetFittedAmplitude(Int32 subeIdx)
{
    return m_FittedAmplitude;
}

/**
 * \brief Returns m_FittedAmplitudeErrorSigma.
 */
Float64 CSingleLine::GetFittedAmplitudeErrorSigma(Int32 subeIdx)
{
    return m_FittedAmplitudeErrorSigma;
}

/**
 * \brief Returns 1.0.
 */
Float64 CSingleLine::GetNominalAmplitude(Int32 subeIdx){
    return 1.0;
}

/**
 * \brief Returns 1.0.
 */
bool CSingleLine::SetNominalAmplitude(Int32 subeIdx, Float64 nominalamp){
    return true;
}

/**
 * \brief Returns m_FittedAmplitude.
 */
Float64 CSingleLine::GetElementAmplitude(){
    return m_FittedAmplitude;
}

/**
 * \brief Sets the amplitude and sigma members.
 * If outside the lambda range, set invalid amplitude and sigma.
 * Else, set amplitude to A (or 0 if A is negative) and sigma to SNR.
 */
void CSingleLine::SetFittedAmplitude(Float64 A, Float64 SNR)
{
  if( m_OutsideLambdaRange )
    {
      m_FittedAmplitude = -1;
      m_FittedAmplitudeErrorSigma = -1;
    }
  else
    {
      m_FittedAmplitude = std::max(0.0, A);
      m_FittedAmplitudeErrorSigma = SNR;
    }
}

/**
 * \brief If the current amplitude is above the argument limit, update it to the limit (or zero if the limit is negative).
 */
void CSingleLine::LimitFittedAmplitude( Int32 subeIdx, Float64 limit )
{
  if( m_FittedAmplitude>limit )
    {
      m_FittedAmplitude = std::max(0.0, limit);
    }
}

/**
 * \brief Returns -1 if the argument is not a line name, and 0 if it is.
 */
Int32 CSingleLine::FindElementIndex(std::string LineTagStr)
{
    Int32 idx = -1;

    std::string name = m_Ray.GetName();
    std::size_t foundstra = name.find(LineTagStr.c_str());

    if (foundstra!=std::string::npos)
      {
        idx = 0;
      }

    return idx;
}

/**
 * \brief Returns m_OutsideLambdaRange.
 */
bool CSingleLine::IsOutsideLambdaRange(Int32 subeIdx)
{
    return m_OutsideLambdaRange;
}
