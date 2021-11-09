// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================

#include "RedshiftLibrary/linemodel/element.h"
#include "RedshiftLibrary/linemodel/multiline.h"

#include "RedshiftLibrary/debug/assert.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/log/log.h"

#include <float.h>
#include <algorithm>

using namespace NSEpic;

/**
 * \brief Constructs the object setting memebers according to arguments and defaults.
 **/
CMultiLine::CMultiLine( std::vector<CRay> rs,
                        const std::string& widthType,
                        const Float64 velocityEmission,
                        const Float64 velocityAbsorption,
                        std::vector<Float64> nominalAmplitudes,
                        Float64 nominalWidth,//corresponds to the lsf of type constant width
                        std::vector<UInt32> catalogIndexes ) : 
CLineModelElement ( widthType, velocityEmission, velocityAbsorption ),
m_NominalAmplitudes(nominalAmplitudes)
{
    //TODO:below variables should be initialized throw CLineModelElt rather than here
    m_Rays = rs;
    m_ElementType = "CMultiLine";
    m_absLinesLimit=-1.0;//-1: disable the ABS lines amplitude cut, any other value is used as a limit for the abs line coeff (typically: 1.0)
    m_sumCross = 0.;
    m_sumGauss = 0.;
    m_dtmFree = 0.;
    m_NominalWidth = nominalWidth;

    Int32 nRays = m_Rays.size();
    m_SignFactors.resize(nRays);
    for(Int32 i=0; i<nRays; i++)
    {
        if( m_Rays[i].GetType()==CRay::nType_Emission )
            m_SignFactors[i] = 1.0;
        else
            m_SignFactors[i] = -1.0;
    }

    for(int i=0; i<catalogIndexes.size(); i++)
    {
        m_LineCatalogIndexes.push_back(catalogIndexes[i]);
    }

    if(m_FittedAmplitudes.size()!=nRays)
    {
        m_FittedAmplitudes.resize(nRays);
        m_FittedAmplitudeErrorSigmas.resize(nRays);
        m_profile.resize(nRays);
    }
    for(Int32 k2=0; k2<nRays; k2++)
    {
        m_profile[k2] = m_Rays[k2].GetProfile();
        if(m_profile[k2]->isAsymFit())
            m_asymLineIndices.push_back(k2);
    }
    SetFittedAmplitude(-1.0, -1.0);
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
TInt32Range CMultiLine::EstimateTheoreticalSupport(Int32 subeIdx, const CSpectrumSpectralAxis& spectralAxis, Float64 redshift, const TFloat64Range &lambdaRange)
{
    Float64 mu = GetObservedPosition(subeIdx, redshift);
    Float64 c = GetLineWidth(mu, redshift, m_Rays[subeIdx].GetIsEmission());
    Float64 winsize = m_profile[subeIdx]->GetNSigmaSupport()*c;

    TInt32Range supportRange = EstimateIndexRange(spectralAxis, mu, lambdaRange, winsize);

    m_StartTheoretical[subeIdx] = supportRange.GetBegin();
    m_EndTheoretical[subeIdx] = supportRange.GetEnd();
    m_StartNoOverlap[subeIdx] = supportRange.GetBegin();
    m_EndNoOverlap[subeIdx] = supportRange.GetEnd();

    if(supportRange.GetBegin()>supportRange.GetEnd()) //in this case the line is completely outside the lambdarange
    {
        m_OutsideLambdaRangeList[subeIdx] = true;
    }else{  //in this case the line is completely inside the lambdarange or with partial overlap

        Int32 minLineOverlap = m_OutsideLambdaRangeOverlapThreshold*winsize;
        Float64 startLbda = spectralAxis[m_StartNoOverlap[subeIdx]];
        Float64 endLbda = spectralAxis[m_EndNoOverlap[subeIdx]];

        if( startLbda >= (lambdaRange.GetEnd()-minLineOverlap) || endLbda<=(lambdaRange.GetBegin()+minLineOverlap) ){
            m_OutsideLambdaRangeList[subeIdx] = true;
        }else{
            m_OutsideLambdaRangeList[subeIdx] = false;
        }
    }

    return supportRange;
}

/**
 * \brief Returns the index range for a given window size (Angstrom)
 **/
TInt32Range CMultiLine::EstimateIndexRange(const CSpectrumSpectralAxis& spectralAxis, Float64 mu,  const TFloat64Range &lambdaRange, Float64 winsizeAngstrom)
{
    TInt32Range supportRange;
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

    //correct the end value if higher then lambdaRange end
    //Log.LogDebug( "    multiline: spectralAxis[supportRange.GetEnd()] = %f", spectralAxis[supportRange.GetEnd()]);
    if(spectralAxis[supportRange.GetEnd()]>lambdaRange.GetEnd())
    {
        supportRange.SetEnd(supportRange.GetEnd()-1);
    }
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
    Int32 nRays = m_Rays.size();
    m_OutsideLambdaRange = true;
    m_RayIsActiveOnSupport.resize( nRays , std::vector<Int32>( nRays , 0.0 ) );
    m_StartNoOverlap.resize(nRays);
    m_EndNoOverlap.resize(nRays);
    m_StartTheoretical.resize(nRays);
    m_EndTheoretical.resize(nRays);
    m_OutsideLambdaRangeList.resize(nRays);
    for(Int32 i=0; i<nRays; i++){
        TInt32Range supportRange = EstimateTheoreticalSupport(i, spectralAxis, redshift, lambdaRange);

        // set the global outside lambda range
        m_OutsideLambdaRange = m_OutsideLambdaRange && m_OutsideLambdaRangeList[i];

        //set the rays active on their own support
        m_RayIsActiveOnSupport[i][i] = 1;
    }

    bool supportNoOverlap_has_duplicates = true;
    Int32 x1=0;
    Int32 y1=0;
    Int32 x2=0;
    Int32 y2=0;
    Int32 icmpt=0;
    Int32 ncmpt=20;
    while(supportNoOverlap_has_duplicates && icmpt<ncmpt)
    {
        icmpt++;
        for(Int32 i=0; i<nRays; i++){
            if(m_OutsideLambdaRangeList[i]){
                continue;
            }
            if(m_StartNoOverlap[i]>m_EndNoOverlap[i])
            {
                continue;
            }
            for(Int32 j=0; j<nRays; j++){
                if(m_OutsideLambdaRangeList[j]){
                    continue;
                }
                if(m_StartNoOverlap[j]>m_EndNoOverlap[j]){
                    continue;
                }
                if(i==j){
                    continue;
                }
                bool rayActiveSupportToBeCorrected = false;
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
                    for(Int32 i2=0; i2<nRays; i2++){
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
                    for(Int32 j2=0; j2<nRays; j2++){
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
        for(Int32 i=0; i<nRays; i++)
        {
            if(supportNoOverlap_has_duplicates)
            {
                break;
            }
            if(m_OutsideLambdaRangeList[i])
            {
                continue;
            }
            for(Int32 j=0; j<nRays; j++)
            {
                if(m_OutsideLambdaRangeList[j]){
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
                    supportNoOverlap_has_duplicates = true;
                    break;
                }
            }
        }
    }

    //init the fitted amplitude values
    for(Int32 k=0; k<nRays; k++){
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
    Int32 nRays = m_Rays.size();
    if(m_OutsideLambdaRange == false)
    {
        for(Int32 i=0; i<nRays; i++)
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

TInt32RangeList CMultiLine::getTheoreticalSupport()
{
    TInt32RangeList support;
    Int32 nRays = m_Rays.size();

    if(m_OutsideLambdaRange == false)
    {
        for(Int32 i=0; i<nRays; i++)
        {
            if(m_OutsideLambdaRangeList[i])
            {
                continue;
            }
            support.push_back(TInt32Range(m_StartTheoretical[i], m_EndTheoretical[i]));
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
    Float64 mu = GetObservedPosition(subeIdx, redshift, false);
    Float64 c = GetLineWidth(mu, redshift, m_Rays[subeIdx].GetIsEmission());
    return c;
}

/**
 * \brief Get the observed position of the sub-element subeIdx for a given redshift
 **/
Float64 CMultiLine::GetObservedPosition(Int32 subeIdx, Float64 redshift, Bool doAsymfitdelta)
{
    Float64 dzOffset = m_Rays[subeIdx].GetOffset()/m_speedOfLightInVacuum;

    Float64 mu = m_Rays[subeIdx].GetPosition()*(1+redshift)*(1+dzOffset);

    // deals with delta of asym profile
    if (doAsymfitdelta)
    {
        std::shared_ptr<CLineProfile> profile = m_Rays[subeIdx].GetProfile();
        mu -= profile->GetAsymDelta();
    }
    return mu;
}

/**
 * \brief Returns the line profile of the sub-element subIdx at wavelength x, for a given redshift.
 **/
Float64 CMultiLine::GetLineProfileAtRedshift(Int32 subeIdx, Float64 redshift, Float64 x)
{
    Float64 mu = GetObservedPosition(subeIdx, redshift, false); // do not apply Lya asym offset
    Float64 sigma = GetLineWidth(mu, redshift, m_Rays[subeIdx].GetIsEmission());
    return m_profile[subeIdx]->GetLineProfile(x, mu, sigma);
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
 * \brief Returns -1.0 if m_OutsideLambdaRange, and the fitted amplitude / nominal amplitude of the first element otherwise.
 **/
Float64 CMultiLine::GetElementAmplitude()
{
  if( m_OutsideLambdaRange )
  {
      return -1.0;
  }
  for(Int32 k=0; k<m_Rays.size(); k++)
  {
      if(!m_OutsideLambdaRangeList[k] && m_NominalAmplitudes[k]!= 0.0)
      {
          return m_FittedAmplitudes[k]/m_NominalAmplitudes[k];
      }
  }
  return -1.0;
}


/**
 * \brief Returns -1.0 if m_OutsideLambdaRange, and the fitted error / nominal amplitude of the first element otherwise.
 **/
Float64 CMultiLine::GetElementError()
{
  if( m_OutsideLambdaRange )
  {
      return -1.0;
  }
  for(Int32 k=0; k<m_Rays.size(); k++)
  {
      if(!m_OutsideLambdaRangeList[k] && m_NominalAmplitudes[k]!= 0.0)
      {
          return m_FittedAmplitudeErrorSigmas[k]/m_NominalAmplitudes[k];
      }
  }
  return -1.0;
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
            if(A<0.0)
            {
                m_FittedAmplitudes[k] = A;
                m_FittedAmplitudeErrorSigmas[k] = A;
            }
            else
            {
                m_FittedAmplitudes[k] = -1.0 -A;
                m_FittedAmplitudeErrorSigmas[k] = -1.0 -A;
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
            m_FittedAmplitudes[k] = -1.0;
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
 * @brief CMultiLine::fitAmplitudeAndLambdaOffset
 * Fit the amplitudes and a unique Offset for all lines by using the fitAmplitude() in a loop to tabulate the offset
 * @param spectralAxis
 * @param noContinuumfluxAxis
 * @param continuumfluxAxis
 * @param redshift
 * @param lineIdx
 */
void CMultiLine::fitAmplitudeAndLambdaOffset(const CSpectrumSpectralAxis& spectralAxis,
                                             const CSpectrumFluxAxis& noContinuumfluxAxis,
                                             const CSpectrumFluxAxis &continuumfluxAxis,
                                             Float64  redshift,
                                             Int32 lineIdx,
                                             bool enableOffsetFitting,
                                             Float64 step,
                                             Float64 min,
                                             Float64 max)
{
    Int32 nRays = m_Rays.size();
    Int32 nSteps = int((max-min)/step+0.5);

    bool atLeastOneOffsetToFit = false;
    if(enableOffsetFitting)
    {
        for(Int32 iR=0; iR<nRays; iR++)
        {
            //check if the line is to be fitted
            if(m_Rays[iR].GetOffsetFitEnabled())
            {
                atLeastOneOffsetToFit = true;
                break;
            }
        }
    }

    if(!atLeastOneOffsetToFit){
        if(m_verbose)
        {
            Log.LogDebug( "    multiline: no offsets to fit");
        }
        nSteps = 1;
    }else{
        if(m_verbose)
        {
            Log.LogDebug( "    multiline: offsets to fit n=%d", nSteps);
        }
    }

    Float64 bestMerit = DBL_MAX;
    Int32 idxBestMerit = -1;
    for(Int32 iO=0; iO<nSteps; iO++)
    {
        //set offset value
        if(atLeastOneOffsetToFit)
        {
            Float64 offset = min+step*iO;
            for(Int32 iR=0; iR<nRays; iR++)
            {
                if(m_Rays[iR].GetOffsetFitEnabled())
                {
                    m_Rays[iR].SetOffset(offset);
                }
            }
        }

        //fit for this offset
        fitAmplitude(spectralAxis, noContinuumfluxAxis, continuumfluxAxis, redshift, lineIdx );

        //check fitting
        if(atLeastOneOffsetToFit)
        {
            Float64 dtm = GetSumCross();
            Float64 mtm = GetSumGauss();
            Float64 a = GetFitAmplitude();
            Float64 term1 = a*a*mtm;
            Float64 term2 = - 2.*a*dtm;
            Float64 fit = term1 + term2;
            if(fit<bestMerit)
            {
                bestMerit = fit;
                idxBestMerit = iO;
            }
        }
    }

    if(idxBestMerit>=0 && atLeastOneOffsetToFit)
    {
        //set offset value
        if(atLeastOneOffsetToFit)
        {
            Float64 offset = min+step*idxBestMerit;
            if(m_verbose)
            {
                Log.LogDebug( "    multiline: offset best found=%f", offset);
            }
            for(Int32 iR=0; iR<nRays; iR++)
            {
                if(m_Rays[iR].GetOffsetFitEnabled())
                {
                    m_Rays[iR].SetOffset(offset);
                }
            }
        }
        //fit again for this offset
        fitAmplitude(spectralAxis, noContinuumfluxAxis, continuumfluxAxis, redshift, lineIdx );
    }
}

/**
 * \brief Estimates an amplitude and fit the relevant rays using this estimate weighed by each ray's nominal amplitude.
 * Loop for the signal synthesis.
 * Loop for the intervals.
 * A estimation.
 * Loop for the signal synthesis.
 **/
void CMultiLine::fitAmplitude(const CSpectrumSpectralAxis& spectralAxis,
                              const CSpectrumFluxAxis& noContinuumfluxAxis,
                              const CSpectrumFluxAxis &continuumfluxAxis,
                              Float64 redshift,
                              Int32 lineIdx )
{
    Int32 nRays = m_Rays.size();

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
    const CSpectrumNoiseAxis& error = noContinuumfluxAxis.GetError();

    if(m_verbose)
    {
        Log.LogDebug("    error[0]:%e", error[0]);
    }
    const Float64* fluxContinuum = continuumfluxAxis.GetSamples();

    Float64 y = 0.0;
    Float64 x = 0.0;
    Float64 yg = 0.0;
    Float64 c = 1.0;

    Float64 err2 = 0.0;
    Int32 num = 0;
    if(m_verbose)
    {
        Log.LogDebug("    multiline: nLines=%d", nRays);
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
        //#pragma omp parallel for
        for(Int32 i = m_StartNoOverlap[k]; i <= m_EndNoOverlap[k]; i++)
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
                    yg += m_SignFactors[k2] * c * m_NominalAmplitudes[k2] * GetLineProfileAtRedshift(k2, redshift, x);
                }else{
                    yg += m_SignFactors[k2] * m_NominalAmplitudes[k2] * GetLineProfileAtRedshift(k2, redshift, x);
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
        if(m_verbose)
        {
            Log.LogDebug("    multiline failed:     num=%d, mtm=%f", num, m_sumGauss);
            for(Int32 k2=0; k2<nRays; k2++)
            {
                Log.LogDebug("    multiline failed:     subE=%d, nominal_amp=%f", k2, m_NominalAmplitudes[k2]);
            }
        }
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
    Int32 nRays = m_Rays.size();
    for(Int32 k=0; k<nRays; k++)
    { //loop on the interval
        if(m_OutsideLambdaRangeList[k])
        {
            continue;
        }

        if( lineIdx>-1 && !(m_RayIsActiveOnSupport[k][lineIdx]))
        {
            continue;
        }

        for(Int32 i = m_StartNoOverlap[k]; i <= m_EndNoOverlap[k]; i++)
        {
            Float64 lambda = spectral[i];
            Float64 Yi = getModelAtLambda(lambda, redshift, continuumfluxAxis[i], k);
            flux[i] += Yi;
            if(isnan(flux[i])){
                throw GlobalException(INTERNAL_ERROR,Formatter()<<"addToSpectrumModel has a NaN flux Ray"<< k<<": ContinuumFlux "<< continuumfluxAxis[i]<<", ModelAtLambda Yi = "<< Yi<<" for range ["<<m_StartNoOverlap[k]<<", "<< m_EndNoOverlap[k]<<"]");
            }
        }
    }
    return;
}

void CMultiLine::addToSpectrumModelDerivVel(const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, CSpectrumFluxAxis& continuumfluxAxis, Float64 redshift, bool emissionRay)
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
        if((emissionRay ^ m_Rays[k].GetIsEmission())){
            continue;
        }

        for(Int32 i = m_StartNoOverlap[k]; i <= m_EndNoOverlap[k]; i++)
        {

            Float64 x = spectral[i];
            Float64 A = m_FittedAmplitudes[k];
            Float64 mu = GetObservedPosition(k, redshift, false);
            Float64 sigma = GetLineWidth(mu, redshift, m_Rays[k].GetIsEmission());

            if(m_SignFactors[k]==-1){
                flux[i] += m_SignFactors[k] * A * continuumfluxAxis[i] * GetLineProfileDerivVel(m_profile[k], x, mu, sigma, m_Rays[k].GetIsEmission());
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
    Float64 Yi = 0.0;

    Float64 x = lambda;
    Int32 nRays = m_Rays.size();

    for(Int32 k2=0; k2<nRays; k2++) //loop on rays
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
        if(A>0.0)
        {
            if(m_SignFactors[k2]==-1){
                Yi += m_SignFactors[k2] * continuumFlux * A * GetLineProfileAtRedshift(k2, redshift, x);
            }else{
                Yi += m_SignFactors[k2] * A * GetLineProfileAtRedshift(k2, redshift, x);
            }
            if(isnan(Yi)){
                Log.LogError("Ray nb: %d adn GetLineProfileAtRedshift: %f", k2, GetLineProfileAtRedshift(k2, redshift, x));
            }
        }
    }
    return Yi;
}

Float64 CMultiLine::GetModelDerivAmplitudeAtLambda(Float64 lambda, Float64 redshift, Float64 continuumFlux)
{
    if(m_OutsideLambdaRange){
        return 0.0;
    }
    Float64 Yi = 0.0;

    Float64 x = lambda;

    for(Int32 k2=0; k2<m_Rays.size(); k2++) //loop on rays
    {
        if(m_OutsideLambdaRangeList[k2]){
            continue;
        }

        if(m_SignFactors[k2]==-1){
            Yi += m_SignFactors[k2] * m_NominalAmplitudes[k2] * continuumFlux * GetLineProfileAtRedshift(k2, redshift, x);
        }else{
            Yi += m_SignFactors[k2] * m_NominalAmplitudes[k2] * GetLineProfileAtRedshift(k2, redshift, x);
        }
    }
    return Yi;
}

Float64 CMultiLine::GetModelDerivContinuumAmpAtLambda(Float64 lambda, Float64 redshift, Float64 continuumFluxUnscale)
{
    if(m_OutsideLambdaRange){
        return 0.0;
    }
    Float64 Yi = 0.0;

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
        Yi += m_SignFactors[k2] *continuumFluxUnscale* A * GetLineProfileAtRedshift(k2, redshift, x);

    }
    return Yi;
}

/* Given the value of the partial deriv of the flux of this multiline at the given lamda when
 * The continuum is not a variable of z
 */
Float64 CMultiLine::GetModelDerivZAtLambdaNoContinuum(Float64 lambda, Float64 redshift, Float64 continuumFlux)
{
    if(m_OutsideLambdaRange){
        return 0.0;
    }
    Float64 Yi = 0.0;

    Float64 x = lambda;

    for(Int32 k2=0; k2<m_Rays.size(); k2++) //loop on rays
    {
        if(m_OutsideLambdaRangeList[k2]){
            continue;
        }
        Float64 A = m_FittedAmplitudes[k2];
        Float64 mu = GetObservedPosition(k2, redshift, false);
        Float64 dzOffset = m_Rays[k2].GetOffset()/m_speedOfLightInVacuum;
        Float64 lamdba0 = m_Rays[k2].GetPosition() * (1+dzOffset);
        Float64 sigma = GetLineWidth(mu, redshift, m_Rays[k2].GetIsEmission());

        if(m_SignFactors[k2]==1){
            Yi += m_SignFactors[k2] * A * m_profile[k2]->GetLineProfileDerivZ( x, lamdba0, redshift, sigma);
        }else{
            Yi += m_SignFactors[k2] * A * continuumFlux * m_profile[k2]->GetLineProfileDerivZ( x, lamdba0, redshift, sigma);
        }
    }
    return Yi;
}

/* Given the value of the partial deriv of the flux of this multiline at the given lamda when
 * The continuum is a variable of z
 */
Float64 CMultiLine::GetModelDerivZAtLambda(Float64 lambda, Float64 redshift, Float64 continuumFlux, Float64 continuumFluxDerivZ)
{
    if(m_OutsideLambdaRange){
        return 0.0;
    }
    Float64 Yi = 0.0;

    Float64 x = lambda;

    for(Int32 k2=0; k2<m_Rays.size(); k2++) //loop on rays
    {
        if(m_OutsideLambdaRangeList[k2]){
                continue;
        }
        Float64 A = m_FittedAmplitudes[k2];
        Float64 mu = GetObservedPosition(k2, redshift, false);
        Float64 dzOffset = m_Rays[k2].GetOffset()/m_speedOfLightInVacuum;
        Float64 lamdba0 = m_Rays[k2].GetPosition() * (1+dzOffset);
        Float64 sigma = GetLineWidth(mu, redshift, m_Rays[k2].GetIsEmission());

        if(m_SignFactors[k2]==1){
            Yi += m_SignFactors[k2] * A * m_profile[k2]->GetLineProfileDerivZ( x, lamdba0, redshift, sigma);
        }else{
            Yi += m_SignFactors[k2] * A * continuumFlux * m_profile[k2]->GetLineProfileDerivZ( x, lamdba0, redshift, sigma)
                  + m_SignFactors[k2] * A * continuumFluxDerivZ * m_profile[k2]->GetLineProfile( x, mu, sigma);
        }
    }
    return Yi;
}

/**
 * \brief For rays inside lambda range, sets the flux to the continuum flux.
 **/
void CMultiLine::initSpectrumModel(CSpectrumFluxAxis &modelfluxAxis, const CSpectrumFluxAxis &continuumfluxAxis, Int32 lineIdx)
{
    if(m_OutsideLambdaRange)
    {
        return;
    }

    Float64* flux = modelfluxAxis.GetSamples();
    for(Int32 k=0; k<m_Rays.size(); k++)
    { //loop on the interval
        if(m_OutsideLambdaRangeList[k])
            continue;

        if( lineIdx>-1 && !(m_RayIsActiveOnSupport[k][lineIdx]))
            continue;

        for(Int32 i = m_StartNoOverlap[k]; i <= m_EndNoOverlap[k]; i++)
            flux[i] = continuumfluxAxis[i];
    }
    return;
}

/**
 * \brief Returns the index corresponding to the first ray whose GetName method returns LineTagStr.
 **/
Int32 CMultiLine::FindElementIndex(std::string LineTagStr)
{
    Int32 idx = -1;
    Int32 rays = m_Rays.size();
    for( UInt32 iElts=0; iElts<rays; iElts++ )
    {
        std::string name = m_Rays[iElts].GetName();
        std::size_t foundstra = name.find(LineTagStr.c_str());

        if(foundstra!=std::string::npos)
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
