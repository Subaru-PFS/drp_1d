
#include <epic/redshift/linemodel/element.h>
#include <epic/redshift/linemodel/multiline.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/core/log/log.h>

#include <algorithm>

using namespace NSEpic;

/**
 * \brief Constructs the object setting memebers according to arguments and defaults.
 **/
CMultiLine::CMultiLine( std::vector<CRay> rs,
			const std::string& widthType,
			const Float64 resolution,
			const Float64 velocityEmission,
			const Float64 velocityAbsorption,
			std::vector<Float64> nominalAmplitudes,
			Float64 nominalWidth,
			std::vector<Int32> catalogIndexes ) : CLineModelElement ( widthType, resolution, velocityEmission, velocityAbsorption )
{
    m_ElementType = "CMultiLine";
    m_Rays = rs;

    Int32 nRays = m_Rays.size();

    m_SignFactors.resize(nRays);
    for(Int32 i=0; i<nRays; i++)
      {
        if( m_Rays[i].GetType()==CRay::nType_Emission )
	  {
            m_SignFactors[i] = 1.0;
	  }
	else
	  {
            m_SignFactors[i] = -1.0;
	  }
      }
    m_NominalWidth = nominalWidth;
    m_NominalAmplitudes = nominalAmplitudes;

    for(int i=0; i<catalogIndexes.size(); i++)
      {
        m_LineCatalogIndexes.push_back(catalogIndexes[i]);
      }

    if(m_FittedAmplitudes.size()!=nRays)
    {
        m_FittedAmplitudes.resize(nRays);
        m_FittedAmplitudeErrorSigmas.resize(nRays);

        mBuffer_mu.resize(nRays);
        mBuffer_c.resize(nRays);
        m_profile.resize(nRays);
    }
    for(Int32 k2=0; k2<nRays; k2++)
    {
        m_profile[k2] = m_Rays[k2].GetProfile();
    }

    SetFittedAmplitude(-1, -1);
}

/**
 * \brief Empty destructor.
 **/
CMultiLine::~CMultiLine()
{
}

/**
 * \brief If the argument is greater than or equal to the size of m_Rays, returns the string "-1". Otherwise returns a call to the m_Rays GetName.
 **/
std::string CMultiLine::GetRayName(Int32 subeIdx)
{
    if(subeIdx>=m_Rays.size())
    {
        return "-1";
    }
    return m_Rays[subeIdx].GetName();
}

/**
 * \brief Returns the content of the m_SignFactors with index equal to the argument.
 **/
Float64 CMultiLine::GetSignFactor(Int32 subeIdx)
{
    return m_SignFactors[subeIdx];
}

/**
 * \brief Limits each m_Rays element within the argument lambdaRange, and sets the m_FittedAmplitudes to -1.
 * Sets the global outside lambda range.
 * Inits the fitted amplitude values.
 **/
void CMultiLine::prepareSupport(const CSpectrumSpectralAxis& spectralAxis, Float64 redshift, const TFloat64Range &lambdaRange)
{    
    m_OutsideLambdaRange=true;
    m_RayIsActiveOnSupport.resize( m_Rays.size() , std::vector<Int32>( m_Rays.size() , 0.0 ) );
    m_Start.resize(m_Rays.size());
    m_End.resize(m_Rays.size());
    m_OutsideLambdaRangeList.resize(m_Rays.size());
    for(Int32 i=0; i<m_Rays.size(); i++){
        Float64 mu = m_Rays[i].GetPosition()*(1+redshift);
        Float64 c = GetLineWidth(mu, redshift, m_Rays[i].GetIsEmission());
        Float64 winsize = GetNSigmaSupport(m_profile[i])*c;

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

        //set the rays active on their own support
        m_RayIsActiveOnSupport[i][i] = 1;
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
                //set the rays active on the overlapping support
                m_RayIsActiveOnSupport[i][j] = 1;
                m_RayIsActiveOnSupport[j][i] = 1;
                //append all the previously overlapping rays as active on the support
                for(Int32 i2=0; i2<i; i2++){
                    if(m_OutsideLambdaRangeList[i2]){
                        continue;
                    }
                    if(m_RayIsActiveOnSupport[i][i2] == 1)
                    {
                        m_RayIsActiveOnSupport[i2][j] = 1;
                        m_RayIsActiveOnSupport[j][i2] = 1;
                    }
                }
                //append all the previously overlapping rays as active on the support
                for(Int32 j2=0; j2<j; j2++){
                    if(m_OutsideLambdaRangeList[j2]){
                        continue;
                    }
                    if(m_RayIsActiveOnSupport[j][j2] == 1)
                    {
                        m_RayIsActiveOnSupport[j2][i] = 1;
                        m_RayIsActiveOnSupport[i][j2] = 1;
                    }
                }
            }
        }
    }

    //init the fitted amplitude values
    for(Int32 k=0; k<m_Rays.size(); k++){
        if(m_OutsideLambdaRangeList[k]){
            m_FittedAmplitudes[k] = -1.0;
        }
    }
}

/**
 * \brief Creates an empty list of ranges as the return value. If not m_OutsideLambdaRange, for each m_Rays element which is also not outside lambda range, add its support to the return value.
 **/
TInt32RangeList CMultiLine::getSupport()
{
    TInt32RangeList support;
    if(m_OutsideLambdaRange==false)
      {
        for(Int32 i=0; i<m_Rays.size(); i++)
	  {
            if(m_OutsideLambdaRangeList[i])
	      {
                continue;
	      }
            support.push_back(TInt32Range(m_Start[i], m_End[i]));
	  }
      }
    return support;
}

/**
 * \brief Creates an empty list of ranges as the return value. If not m_OutsideLambdaRange, for each m_Rays element belonging to the argument subeIdx which is also not outside lambda range, add its support to the return value.
 **/
TInt32Range CMultiLine::getSupportSubElt(Int32 subeIdx)
{
    TInt32Range support;
    if(m_OutsideLambdaRange==false){
        if(m_OutsideLambdaRangeList[subeIdx]){
            support = TInt32Range(-1, -1);
        }
        support = TInt32Range(m_Start[subeIdx], m_End[subeIdx]);
    }
    return support;
}

/**
 * \brief Calls GetLineWidth using the arguments and a calculated argument mu.
 **/
Float64 CMultiLine::GetWidth(Int32 subeIdx, Float64 redshift)
{
    Float64 mu = m_Rays[subeIdx].GetPosition()*(1+redshift);
    Float64 c = GetLineWidth(mu, redshift, m_Rays[subeIdx].GetIsEmission());
    return c;
}

/**
 * \brief Returns a copy of m_Rays.
 **/
std::vector<CRay> CMultiLine::GetRays()
{
    std::vector<CRay> rays;
    for(Int32 k=0; k<m_Rays.size(); k++){
        rays.push_back(m_Rays[k]);
    }
    return rays;
}

/**
 * \brief Returns the fitted amplitude for the argument index.
 **/
Float64 CMultiLine::GetFittedAmplitude(Int32 subeIdx)
{
    return m_FittedAmplitudes[subeIdx];
}

/**
 * \brief Returns the fitted amplitude error for the argument index.
 **/
Float64 CMultiLine::GetFittedAmplitudeErrorSigma(Int32 subeIdx)
{
    return m_FittedAmplitudeErrorSigmas[subeIdx];
}

/**
 * \brief Returns -1 if m_OutsideLambdaRange, and the fitted amplitude / nominal amplitude of the first element otherwise.
 **/
Float64 CMultiLine::GetElementAmplitude()
{
  if( m_OutsideLambdaRange )
    {
      return -1;
    }
  return m_FittedAmplitudes[0]/m_NominalAmplitudes[0];
}

/**
 * \brief Returns the nominal amplitude with index subeIdx.
 **/
Float64 CMultiLine::GetNominalAmplitude(Int32 subeIdx)
{
    return m_NominalAmplitudes[subeIdx];
}

/**
 * \brief If outside lambda range, sets fitted amplitudes and errors to -1. If inside, sets each ray's fitted amplitude and error to -1 if ray outside lambda range, or amplitude to A * nominal amplitude and error to SNR * nominal amplitude.
 **/
void CMultiLine::SetFittedAmplitude(Float64 A, Float64 SNR)
{
    if(m_OutsideLambdaRange)
      {
        for(Int32 k=0; k<m_Rays.size(); k++)
	  {
	    // This separation in distinct negative values is to facilitate unit testing. All negative amplitudes should be considered invalid.
	    if( A<0 )
	      {
		m_FittedAmplitudes[k] = A;
		m_FittedAmplitudeErrorSigmas[k] = A;
	      }
	    else
	      {
		m_FittedAmplitudes[k] = -1 -A;
		m_FittedAmplitudeErrorSigmas[k] = -1 -A;
	      }
		
	  }
        return;
      }
    A = std::max(0.0, A);
    for(Int32 k=0; k<m_Rays.size(); k++)
      {
	if(m_OutsideLambdaRangeList[k])
	  {
	    m_FittedAmplitudes[k] = -1;
	  }
	m_FittedAmplitudes[k] = A*m_NominalAmplitudes[k];
	m_FittedAmplitudeErrorSigmas[k] = SNR*m_NominalAmplitudes[k]; //todo: check correct formulation for Error
      }
}

/**
 * \brief Estimates an amplitude and fit the relevant rays using this estimate weighed by each ray's nominal amplitude.
 * Loop for the signal synthesis.
 * Loop for the intervals.
 * A estimation.
 * Loop for the signal synthesis.
 **/
void CMultiLine::fitAmplitude(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64  redshift, Int32 lineIdx )
{
    Float64 nRays = m_Rays.size();

    for(Int32 k=0; k<nRays; k++)
      {
        m_FittedAmplitudes[k] = -1.0;
        m_FittedAmplitudeErrorSigmas[k] = -1.0;
      }

    if(m_OutsideLambdaRange)
      {
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



    for(Int32 k2=0; k2<nRays; k2++)
      {
        mBuffer_mu[k2] = m_Rays[k2].GetPosition()*(1+redshift);
        mBuffer_c[k2] = GetLineWidth(mBuffer_mu[k2], redshift, m_Rays[k2].GetIsEmission());
      }

    for(Int32 k=0; k<nRays; k++)
      { //loop for the intervals
        if(m_OutsideLambdaRangeList[k])
	  {
            continue;
	  }
        if( lineIdx>-1 && !(m_RayIsActiveOnSupport[k][lineIdx]))
        {
            continue;
        }

        //A estimation
        for ( Int32 i = m_Start[k]; i <= m_End[k]; i++)
        {
            y = flux[i];
            x = spectral[i];

            yg = 0.0;
            for(Int32 k2=0; k2<nRays; k2++)
            { //loop for the signal synthesis
                if(m_OutsideLambdaRangeList[k2])
                {
                    continue;
                }
                if(m_RayIsActiveOnSupport[k2][k]==0)
                {
                    continue;
                }
                yg += m_SignFactors[k2] * m_NominalAmplitudes[k2] * GetLineProfile(m_profile[k2], x-mBuffer_mu[k2], mBuffer_c[k2]);
            }
            num++;
            err2 = 1.0 / (error[i] * error[i]);
            sumCross += yg*y*err2;
            sumGauss += yg*yg*err2;
        }

      }

    if ( num==0 || sumGauss==0 )
      {
        return;
      }

    Float64 A = std::max(0.0, sumCross / sumGauss);

    for(Int32 k=0; k<nRays; k++)
    {
        if(m_OutsideLambdaRangeList[k])
        {
            continue;
        }
        m_FittedAmplitudes[k] = A*m_NominalAmplitudes[k];
        if(A==0)
        {
            m_FittedAmplitudeErrorSigmas[k] = 0.0;
        }
        else
        {
            m_FittedAmplitudeErrorSigmas[k] = m_NominalAmplitudes[k]*1.0/sqrt(sumGauss); //Achtung! To be discussed with Didier V.
        }
    }
    return;
}

/**
 * \brief Adds to the model's flux, at each ray not outside lambda range, the value contained in the corresponding lambda for each catalog line.
 **/
void CMultiLine::addToSpectrumModel( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift, Int32 lineIdx )
{
    if(m_OutsideLambdaRange)
      {
        return;
      }

    const Float64* spectral = modelspectralAxis.GetSamples();
    Float64* flux = modelfluxAxis.GetSamples();
    for(Int32 k=0; k<m_Rays.size(); k++)
      { //loop on the interval
        if(m_OutsideLambdaRangeList[k])
	  {
            continue;
	  }

        if( lineIdx>-1 && !(m_RayIsActiveOnSupport[k][lineIdx]))
        {
            continue;
        }

        for ( Int32 i = m_Start[k]; i <= m_End[k]; i++)
      {
            Float64 lambda = spectral[i];
            Float64 Yi=getModelAtLambda(lambda, redshift, k);
            flux[i] += Yi;
        }
    }
  return;
}

void CMultiLine::addToSpectrumModelDerivSigma( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift )
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
            Float64 Yi=GetModelDerivSigmaAtLambda(lambda, redshift);
            flux[i] += Yi;
        }
    }
  return;
}

/**
 * \brief Returns the sum of the amplitude of each ray on redshifted lambda.
 **/
Float64 CMultiLine::getModelAtLambda( Float64 lambda, Float64 redshift, Int32 kRaySupport)
{
    if(m_OutsideLambdaRange)
      {
        return 0.0;
      }
    Float64 Yi=0.0;

    Float64 x = lambda;

    for(Int32 k2=0; k2<m_Rays.size(); k2++) //loop on rays
    {
        if(m_OutsideLambdaRangeList[k2])
	  {
            continue;
	  }
        if( kRaySupport>=0 && m_RayIsActiveOnSupport[k2][kRaySupport]==0 )
        {
            continue;
        }
        Float64 A = m_FittedAmplitudes[k2];
        Float64 mu = m_Rays[k2].GetPosition()*(1+redshift);
        Float64 c = GetLineWidth(mu, redshift, m_Rays[k2].GetIsEmission());

        Yi += m_SignFactors[k2] * A * GetLineProfile(m_profile[k2], x-mu, c);
    }
    return Yi;
}

Float64 CMultiLine::GetModelDerivAmplitudeAtLambda(Float64 lambda, Float64 redshift )
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

        Float64 mu = m_Rays[k2].GetPosition()*(1+redshift);
        Float64 c = GetLineWidth(mu, redshift, m_Rays[k2].GetIsEmission());

        Yi += m_SignFactors[k2] * GetLineProfile(m_profile[k2], x-mu, c);
    }
    return Yi;
}

Float64 CMultiLine::GetModelDerivSigmaAtLambda(Float64 lambda, Float64 redshift )
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
        Float64 c = GetLineWidth(mu, redshift, m_Rays[k2].GetIsEmission());

        Yi += m_SignFactors[k2] * A * GetLineProfileDerivSigma(m_profile[k2], x, mu, c);
    }
    return Yi;
}


/**
 * \brief For rays inside lambda range, sets the flux to the continuum flux.
 **/
void CMultiLine::initSpectrumModel( CSpectrumFluxAxis &modelfluxAxis, CSpectrumFluxAxis &continuumfluxAxis, Int32 lineIdx )
{
    if(m_OutsideLambdaRange)
      {
        return;
      }

    Float64* flux = modelfluxAxis.GetSamples();
    for(Int32 k=0; k<m_Rays.size(); k++)
      { //loop on the interval
        if(m_OutsideLambdaRangeList[k])
	  {
            continue;
	  }

        if( lineIdx>-1 && !(m_RayIsActiveOnSupport[k][lineIdx]))
        {
            continue;
        }

        for ( Int32 i = m_Start[k]; i <= m_End[k]; i++)
	  {
            flux[i] = continuumfluxAxis[i];
	  }
      }
    return;
}

/**
 * \brief Returns the index corresponding to the first ray whose GetName method returns LineTagStr.
 **/
Int32 CMultiLine::FindElementIndex(std::string LineTagStr)
{
    Int32 idx = -1;
    for( UInt32 iElts=0; iElts<m_Rays.size(); iElts++ )
    {
        std::string name = m_Rays[iElts].GetName();
        std::size_t foundstra = name.find(LineTagStr.c_str());

        if (foundstra!=std::string::npos)
	  {
            idx = iElts;
            break;
	  }
    }
    return idx;
}

/**
 * \brief If the fitted amplitude of ray with index subeIdx is above the limit, sets it to either that limit or zero, whichever is greater.
 **/
void CMultiLine::LimitFittedAmplitude(Int32 subeIdx, Float64 limit)
{

    if(m_FittedAmplitudes[subeIdx] > limit)
      {
        m_FittedAmplitudes[subeIdx] = std::max(0.0, limit);

        //now update the amplitude of the other lines
        Float64 amplitudeRef = m_FittedAmplitudes[subeIdx]/m_NominalAmplitudes[subeIdx];
        for(Int32 k=0; k<m_Rays.size(); k++)
        {
            m_FittedAmplitudes[k] = m_NominalAmplitudes[k]*amplitudeRef;
        }
      }
    return;
}

/**
 * \brief Returns whether the ray with index subeIdx is outside the lambda range.
 **/
bool CMultiLine::IsOutsideLambdaRange(Int32 subeIdx)
{
    return m_OutsideLambdaRangeList[subeIdx];
}
