
#include <RedshiftLibrary/linemodel/element.h>
#include <RedshiftLibrary/linemodel/multiline.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/log/log.h>

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
    m_c_kms = 300000.0; //to be defined in a better location ?

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

    m_absLinesLimit = -1.0; //-1: disable the ABS lines amplitude cut, any other value is used as a limit for the abs line coeff (typically: 1.0)

    m_sumCross = 0.0;
    m_sumGauss = 0.0;
    m_dtmFree = 0.0;

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

bool CMultiLine::SetAbsLinesLimit(Float64 limit)
{
    m_absLinesLimit = limit;
    return true;
}

/**
 * @brief GetContinuumAtCenterProfile
 * @param subeIdx
 * @param spectralAxis
 * @param redshift
 * @param lambdaRange
 * @param continuumfluxAxis
 * @return the continuum flux val at the sub element center wavelength. Error returns -999/-9999 if center profile not in range
 *
 */
Float64 CMultiLine::GetContinuumAtCenterProfile(Int32 subeIdx, const CSpectrumSpectralAxis& spectralAxis, Float64 redshift, CSpectrumFluxAxis &continuumfluxAxis)
{
    Float64 mu = GetObservedPosition(subeIdx, redshift);

    Int32 IdxCenterProfile = spectralAxis.GetIndexAtWaveLength(mu);
    if(IdxCenterProfile<0 || IdxCenterProfile>continuumfluxAxis.GetSamplesCount()-1)
    {
        return -9999.0;
    }

    Float64 cont = continuumfluxAxis[IdxCenterProfile];

    return cont;
}


/**
 * \brief Returns the theoretical support range for the line
 **/
TInt32Range CMultiLine::EstimateTheoreticalSupport(Int32 subeIdx, const CSpectrumSpectralAxis& spectralAxis, Float64 redshift,  const TFloat64Range &lambdaRange)
{
    Float64 mu = GetObservedPosition(subeIdx, redshift);
    Float64 c = GetLineWidth(mu, redshift, m_Rays[subeIdx].GetIsEmission(), m_profile[subeIdx]);
    Float64 winsize = GetNSigmaSupport(m_profile[subeIdx])*c;

    TInt32Range supportRange = EstimateIndexRange(subeIdx, spectralAxis, redshift, lambdaRange, winsize);

    return supportRange;
}

/**
 * \brief Returns the index range for a given window size (Angstrom)
 **/
TInt32Range CMultiLine::EstimateIndexRange(Int32 subeIdx, const CSpectrumSpectralAxis& spectralAxis, Float64 redshift,  const TFloat64Range &lambdaRange, Float64 winsizeAngstrom)
{
    TInt32Range supportRange;
    Float64 mu = GetObservedPosition(subeIdx, redshift);
    Float64 winsize = winsizeAngstrom;

    Float64 lambda_start = mu-winsize/2.0;
    if(lambda_start < lambdaRange.GetBegin()){
        lambda_start = lambdaRange.GetBegin();
    }
    supportRange.SetBegin(spectralAxis.GetIndexAtWaveLength(lambda_start));

    Float64 lambda_end = mu+winsize/2.0;
    if(lambda_end > lambdaRange.GetEnd()){
        lambda_end = lambdaRange.GetEnd();
    }
    supportRange.SetEnd(spectralAxis.GetIndexAtWaveLength(lambda_end));

    //correct the end value if not higher or equal to the begin value
    if(supportRange.GetEnd()<supportRange.GetBegin())
    {
        supportRange.SetEnd(supportRange.GetBegin()-1);
    }

    return supportRange;
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
    m_StartNoOverlap.resize(m_Rays.size());
    m_EndNoOverlap.resize(m_Rays.size());
    m_StartTheoretical.resize(m_Rays.size());
    m_EndTheoretical.resize(m_Rays.size());
    m_OutsideLambdaRangeList.resize(m_Rays.size());
    for(Int32 i=0; i<m_Rays.size(); i++){
        TInt32Range supportRange = EstimateTheoreticalSupport(i, spectralAxis, redshift, lambdaRange);
        m_StartTheoretical[i] = supportRange.GetBegin();
        m_EndTheoretical[i] = supportRange.GetEnd();
        m_StartNoOverlap[i] = supportRange.GetBegin();
        m_EndNoOverlap[i] = supportRange.GetEnd();

        if(supportRange.GetBegin()>supportRange.GetEnd()) //in this case the line is completely outside the lambdarange
        {
            m_OutsideLambdaRangeList[i]=true;
        }else{  //in this case the line is completely inside the lambdarange or with partial overlap

            Float64 mu = GetObservedPosition(i, redshift);
            Float64 c = GetLineWidth(mu, redshift, m_Rays[i].GetIsEmission(), m_profile[i]);
            Float64 winsize = GetNSigmaSupport(m_profile[i])*c;
            Int32 minLineOverlap = m_OutsideLambdaRangeOverlapThreshold*winsize;
            Float64 startLbda = spectralAxis[m_StartNoOverlap[i]];
            Float64 endLbda = spectralAxis[m_EndNoOverlap[i]];

            if( startLbda >= (lambdaRange.GetEnd()-minLineOverlap) || endLbda<=(lambdaRange.GetBegin()+minLineOverlap) ){
                m_OutsideLambdaRangeList[i]=true;
            }else{
                m_OutsideLambdaRangeList[i]=false;
            }
        }

        // set the global outside lambda range
        m_OutsideLambdaRange = m_OutsideLambdaRange && m_OutsideLambdaRangeList[i];

        //set the rays active on their own support
        m_RayIsActiveOnSupport[i][i] = 1;
    }

    bool supportNoOverlap_has_duplicates=true;
    Int32 x1=0;
    Int32 y1=0;
    Int32 x2=0;
    Int32 y2=0;
    Int32 icmpt=0;
    Int32 ncmpt=20;
    while(supportNoOverlap_has_duplicates && icmpt<ncmpt)
    {
        icmpt++;
        for(Int32 i=0; i<m_Rays.size(); i++){
            if(m_OutsideLambdaRangeList[i]){
                continue;
            }
            if(m_StartNoOverlap[i]>m_EndNoOverlap[i])
            {
                continue;
            }
            for(Int32 j=0; j<m_Rays.size(); j++){
                if(m_OutsideLambdaRangeList[j]){
                    continue;
                }
                if(m_StartNoOverlap[j]>m_EndNoOverlap[j])
                {
                    continue;
                }
                if(i==j){
                    continue;
                }
                bool rayActiveSupportToBeCorrected=false;
                //
                x1 = m_StartNoOverlap[i];
                x2 = m_EndNoOverlap[i];
                y1 = m_StartNoOverlap[j];
                y2 = m_EndNoOverlap[j];
                Int32 max = std::max(x1,y1);
                Int32 min = std::min(x2,y2);
                if( max-min < 0 )
                {
                    m_StartNoOverlap[i]=std::min(x1,y1);
                    m_EndNoOverlap[i]=std::max(x2,y2);
                    m_StartNoOverlap[j]=m_EndNoOverlap[i]; //deactivate j when end is start -1
                    m_EndNoOverlap[j]=m_EndNoOverlap[i]-1; //deactivate j

                    rayActiveSupportToBeCorrected = true;
                }

                if(rayActiveSupportToBeCorrected)
                {
                    //set the rays active on the overlapping support
                    m_RayIsActiveOnSupport[i][j] = 1;
                    m_RayIsActiveOnSupport[j][i] = 1;
                    //append all the previously overlapping rays as active on the support
                    for(Int32 i2=0; i2<m_Rays.size(); i2++){
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
                    for(Int32 j2=0; j2<m_Rays.size(); j2++){
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

        supportNoOverlap_has_duplicates = false;

        //check that there are no overlapping sub-supports in the list
        for(Int32 i=0; i<m_Rays.size(); i++)
        {
            if(supportNoOverlap_has_duplicates)
            {
                break;
            }
            if(m_OutsideLambdaRangeList[i])
            {
                continue;
            }
            for(Int32 j=0; j<m_Rays.size(); j++)
            {
                if(m_OutsideLambdaRangeList[j])
                {
                    continue;
                }
                if(i==j){
                    continue;
                }

                x1 = m_StartNoOverlap[i];
                x2 = m_EndNoOverlap[i];
                y1 = m_StartNoOverlap[j];
                y2 = m_EndNoOverlap[j];
                Int32 max = std::max(x1,y1);
                Int32 min = std::min(x2,y2);
                if( max-min < 0 )
                {
                    supportNoOverlap_has_duplicates=true;
                    break;
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
            support.push_back(TInt32Range(m_StartNoOverlap[i], m_EndNoOverlap[i]));
	  }
      }
    return support;
}

/**
 * \brief Creates an empty list of ranges as the return value. If not m_OutsideLambdaRange, for each m_Rays element belonging to the argument subeIdx which is also not outside lambda range, add its support to the return value.
 **/
TInt32Range CMultiLine::getSupportSubElt(Int32 subeIdx)
{
    /*
    TInt32Range support;
    if(m_OutsideLambdaRange==false){
        if( !m_OutsideLambdaRangeList[subeIdx] ){
            support = TInt32Range(m_StartNoOverlap[subeIdx], m_EndNoOverlap[subeIdx]);
        }
    }
    */
    TInt32Range support = TInt32Range(m_StartNoOverlap[subeIdx], m_EndNoOverlap[subeIdx]);
    return support;
}

/**
 * \brief Returns the theoretical support of the line (sub-element).
 **/
TInt32Range CMultiLine::getTheoreticalSupportSubElt(Int32 subeIdx)
{
    /*
    TInt32Range support;
    if(m_OutsideLambdaRange==false){
        if( !m_OutsideLambdaRangeList[subeIdx] ){
            support = TInt32Range(m_StartTheoretical[subeIdx], m_EndTheoretical[subeIdx]);
        }
    }
    */
    TInt32Range support = TInt32Range(m_StartTheoretical[subeIdx], m_EndTheoretical[subeIdx]);
    return support;
}

/**
 * \brief Calls GetLineWidth using the arguments and a calculated argument mu.
 **/
Float64 CMultiLine::GetWidth(Int32 subeIdx, Float64 redshift)
{
    Float64 mu = GetObservedPosition(subeIdx, redshift);
    Float64 c = GetLineWidth(mu, redshift, m_Rays[subeIdx].GetIsEmission(), m_profile[subeIdx]);
    return c;
}

/**
 * \brief Get the observed position of the sub-element subeIdx for a given redshift
 **/
Float64 CMultiLine::GetObservedPosition(Int32 subeIdx, Float64 redshift)
{
    Float64 dzOffset = m_Rays[subeIdx].GetOffset()/m_c_kms;
    Float64 mu = m_Rays[subeIdx].GetPosition()*(1+redshift)*(1+dzOffset);
    return mu;
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
  for(Int32 k=0; k<m_Rays.size(); k++)
  {
      if(!m_OutsideLambdaRangeList[k] && m_NominalAmplitudes[k]!= 0.0)
      {
          return m_FittedAmplitudes[k]/m_NominalAmplitudes[k];
      }
  }
  return -1;
}

/**
 * \brief Returns the nominal amplitude with index subeIdx.
 **/
Float64 CMultiLine::GetNominalAmplitude(Int32 subeIdx)
{
    return m_NominalAmplitudes[subeIdx];
}

/**
 * \brief Set the nominal amplitude with index subeIdx.
 **/
bool CMultiLine::SetNominalAmplitude(Int32 subeIdx, Float64 nominalamp)
{
    m_NominalAmplitudes[subeIdx]=nominalamp;
    return true;
}

//WARNING: setfittedamplitude should not be applied to sub element individually ?
void CMultiLine::SetFittedAmplitude(Int32 subeIdx, Float64 A, Float64 SNR)
{
    if(m_OutsideLambdaRangeList[subeIdx])
    {
        m_FittedAmplitudes[subeIdx] = -1;
    }
    m_FittedAmplitudes[subeIdx] = A*m_NominalAmplitudes[subeIdx];
    // limit the absorption to 0.0-1.0, so that it's never <0
    //*
    if(m_SignFactors[subeIdx]==-1 && m_absLinesLimit>0.0 && m_FittedAmplitudes[subeIdx]>m_absLinesLimit){
        m_FittedAmplitudes[subeIdx]=m_absLinesLimit;
    }
    //*/
    m_FittedAmplitudeErrorSigmas[subeIdx] = SNR*m_NominalAmplitudes[subeIdx]; //todo: check correct formulation for Error

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
    m_fitAmplitude = A;
    for(Int32 k=0; k<m_Rays.size(); k++)
    {
        if(m_OutsideLambdaRangeList[k])
        {
            m_FittedAmplitudes[k] = -1;
        }
        m_FittedAmplitudes[k] = A*m_NominalAmplitudes[k];
        // limit the absorption to 0.0-1.0, so that it's never <0
        //*
        if(m_SignFactors[k]==-1 && m_absLinesLimit>0.0 && m_FittedAmplitudes[k]>m_absLinesLimit){
            m_FittedAmplitudes[k]=m_absLinesLimit;
        }
        //*/
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
void CMultiLine::fitAmplitude(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& noContinuumfluxAxis, const CSpectrumFluxAxis &continuumfluxAxis, Float64  redshift, Int32 lineIdx )
{
    Float64 nRays = m_Rays.size();

    m_sumCross = 0.0;
    m_sumGauss = 0.0;
    m_dtmFree = 0.0;

    for(Int32 k=0; k<nRays; k++)
      {
        m_FittedAmplitudes[k] = -1.0;
        m_FittedAmplitudeErrorSigmas[k] = -1.0;
      }

    if(m_OutsideLambdaRange)
      {
        return;
      }
    const Float64* fluxNoContinuum = noContinuumfluxAxis.GetSamples();
    const Float64* spectral = spectralAxis.GetSamples();
    const Float64* error = noContinuumfluxAxis.GetError();
    const Float64* fluxContinuum = continuumfluxAxis.GetSamples();

    Float64 y = 0.0;
    Float64 x = 0.0;
    Float64 yg = 0.0;
    Float64 c = 1.0;


    Float64 err2 = 0.0;
    Int32 num = 0;


    for(Int32 k2=0; k2<nRays; k2++)
      {
        mBuffer_mu[k2] = GetObservedPosition(k2, redshift);
        mBuffer_c[k2] = GetLineWidth(mBuffer_mu[k2], redshift, m_Rays[k2].GetIsEmission(), m_profile[k2]);
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
        for ( Int32 i = m_StartNoOverlap[k]; i <= m_EndNoOverlap[k]; i++)
        {
            c = fluxContinuum[i];
            y = fluxNoContinuum[i];
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

                if(m_SignFactors[k2]==-1){
                    yg += m_SignFactors[k2] * c * m_NominalAmplitudes[k2] * GetLineProfile(m_profile[k2], x, mBuffer_mu[k2], mBuffer_c[k2]);
                }else{
                    yg += m_SignFactors[k2] * m_NominalAmplitudes[k2] * GetLineProfile(m_profile[k2], x, mBuffer_mu[k2], mBuffer_c[k2]);
                }
            }
            num++;
            err2 = 1.0 / (error[i] * error[i]);
            m_dtmFree += yg*y*err2;
            m_sumGauss += yg*yg*err2;
        }

      }

    if ( num==0 || m_sumGauss==0 )
      {
        return;
      }

    m_sumCross = std::max(0.0, m_dtmFree);
    Float64 A = m_sumCross / m_sumGauss;
    m_fitAmplitude = A; //todo: warning m_fitAmplitude should be updated when modifying sub-elements amplitudes: ex. rules.
    //Float64 A = std::max(0.0, m_sumCross / m_sumGauss);

    for(Int32 k=0; k<nRays; k++)
    {
        if(m_OutsideLambdaRangeList[k])
        {
            continue;
        }
        m_FittedAmplitudes[k] = A*m_NominalAmplitudes[k];

        // limit the absorption to 0.0-1.0, so that it's never <0
        //*
        if(m_SignFactors[k]==-1 && m_absLinesLimit>0.0 && m_FittedAmplitudes[k]>m_absLinesLimit){
            m_FittedAmplitudes[k]=m_absLinesLimit;
        }
        //*/

//        if(A==0)
//        {
//            m_FittedAmplitudeErrorSigmas[k] = 0.0; //why would this be useful ?
//        }
//        else
        {
            m_FittedAmplitudeErrorSigmas[k] = m_NominalAmplitudes[k]*1.0/sqrt(m_sumGauss);
        }
    }
    return;
}

/**
 * \brief Adds to the model's flux, at each ray not outside lambda range, the value contained in the corresponding lambda for each catalog line.
 **/
void CMultiLine::addToSpectrumModel( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, CSpectrumFluxAxis &continuumfluxAxis, Float64 redshift, Int32 lineIdx )
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

        for ( Int32 i = m_StartNoOverlap[k]; i <= m_EndNoOverlap[k]; i++)
      {
            Float64 lambda = spectral[i];
            Float64 Yi=getModelAtLambda(lambda, redshift, continuumfluxAxis[i], k);
            flux[i] += Yi;
        }
    }
  return;
}

void CMultiLine::addToSpectrumModelDerivVel( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis,  CSpectrumFluxAxis& continuumfluxAxis, Float64 redshift , bool emissionRay)
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
        if((emissionRay  ^ m_Rays[k].GetIsEmission())){
            continue;
        }

        for ( Int32 i = m_StartNoOverlap[k]; i <= m_EndNoOverlap[k]; i++)
        {

            Float64 x = spectral[i];
            Float64 A = m_FittedAmplitudes[k];
            Float64 mu = GetObservedPosition(k, redshift);
            Float64 sigma = GetLineWidth(mu, redshift, m_Rays[k].GetIsEmission(), m_profile[k]);

            if(m_SignFactors[k]==-1){
                flux[i] += m_SignFactors[k] * A * continuumfluxAxis[i] * GetLineProfileDerivVel(m_profile[k], x, mu, sigma,  m_Rays[k].GetIsEmission());
            }else{
                flux[i] += m_SignFactors[k] * A * GetLineProfileDerivVel(m_profile[k], x, mu, sigma,  m_Rays[k].GetIsEmission());
            }

        }
    }
  return;
}

/**
 * \brief Returns the sum of the amplitude of each ray on redshifted lambda.
 **/
Float64 CMultiLine::getModelAtLambda(Float64 lambda, Float64 redshift, Float64 continuumFlux, Int32 kRaySupport)
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
        Float64 mu = GetObservedPosition(k2, redshift);
        Float64 sigma = GetLineWidth(mu, redshift, m_Rays[k2].GetIsEmission(), m_profile[k2]);

        //*
        if(m_SignFactors[k2]==-1){
            Yi += m_SignFactors[k2] * continuumFlux * A * GetLineProfile(m_profile[k2], x, mu, sigma);
        }else{
            Yi += m_SignFactors[k2] * A * GetLineProfile(m_profile[k2], x, mu, sigma);
        }
        //*/
    }
    return Yi;
}

Float64 CMultiLine::GetModelDerivAmplitudeAtLambda(Float64 lambda, Float64 redshift, Float64 continuumFlux )
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

        Float64 mu = GetObservedPosition(k2, redshift);
	Float64 sigma = GetLineWidth(mu, redshift, m_Rays[k2].GetIsEmission(), m_profile[k2]);

				if(m_SignFactors[k2]==-1){
            Yi += m_SignFactors[k2] * m_NominalAmplitudes[k2] * continuumFlux * GetLineProfile(m_profile[k2], x, mu, sigma);
        }else{
            Yi += m_SignFactors[k2] * m_NominalAmplitudes[k2] * GetLineProfile(m_profile[k2], x, mu, sigma);
        }
    }
    return Yi;
}

Float64  CMultiLine::GetModelDerivContinuumAmpAtLambda(Float64 lambda, Float64 redshift, Float64 continuumFluxUnscale  ){
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

			if(m_SignFactors[k2]==1){
				continue;
			}

			Float64 A = m_FittedAmplitudes[k2];
			Float64 dzOffset = m_Rays[k2].GetOffset()/m_c_kms;
			Float64 mu = m_Rays[k2].GetPosition()*(1+redshift)*(1+dzOffset);
			Float64 sigma = GetLineWidth(mu, redshift, m_Rays[k2].GetIsEmission(), m_profile[k2]);

			Yi += m_SignFactors[k2] *continuumFluxUnscale* A * GetLineProfile(m_profile[k2], x, mu, sigma);

	}
	return Yi;
}
/* Given the value of the partial deriv of the flux of this multiline at the given lamda when
 * The continuum is not a variable of z
 */
Float64 CMultiLine::GetModelDerivZAtLambdaNoContinuum(Float64 lambda, Float64 redshift, Float64 continuumFlux){
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
        Float64 mu = GetObservedPosition(k2, redshift);
        Float64 dzOffset = m_Rays[k2].GetOffset()/m_c_kms;
        Float64 lamdba0 = m_Rays[k2].GetPosition() * (1+dzOffset);
        Float64 sigma = GetLineWidth(mu, redshift, m_Rays[k2].GetIsEmission(), m_profile[k2]);

        if(m_SignFactors[k2]==1){
            Yi += m_SignFactors[k2] * A * GetLineProfileDerivZ(m_profile[k2], x, lamdba0, redshift, sigma);
        }else{
            Yi += m_SignFactors[k2] * A * continuumFlux * GetLineProfileDerivZ(m_profile[k2], x, lamdba0, redshift, sigma);
        }
    }
    return Yi;
}
/* Given the value of the partial deriv of the flux of this multiline at the given lamda when
 * The continuum is a variable of z
 */
Float64 CMultiLine::GetModelDerivZAtLambda(Float64 lambda, Float64 redshift, Float64 continuumFlux, Float64 continuumFluxDerivZ){
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
        Float64 mu = GetObservedPosition(k2, redshift);
        Float64 dzOffset = m_Rays[k2].GetOffset()/m_c_kms;
        Float64 lamdba0 = m_Rays[k2].GetPosition() * (1+dzOffset);
        Float64 sigma = GetLineWidth(mu, redshift, m_Rays[k2].GetIsEmission(), m_profile[k2]);

        if(m_SignFactors[k2]==1){
            Yi += m_SignFactors[k2] * A * GetLineProfileDerivZ(m_profile[k2], x, lamdba0, redshift, sigma);
        }else{
            Yi += m_SignFactors[k2] * A * continuumFlux * GetLineProfileDerivZ(m_profile[k2], x, lamdba0, redshift, sigma)
                    + m_SignFactors[k2] * A * continuumFluxDerivZ * GetLineProfile(m_profile[k2], x, mu, sigma);
        }
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

        for ( Int32 i = m_StartNoOverlap[k]; i <= m_EndNoOverlap[k]; i++)
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
