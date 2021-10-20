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
#include "RedshiftLibrary/operator/templatefittingBase.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
using namespace NSEpic;
using namespace std;

COperatorTemplateFittingBase::COperatorTemplateFittingBase(const CSpectrum& spectrum, const TFloat64Range& lambdaRange, const TFloat64List & redshifts):
    m_spectrum(spectrum),
    m_redshifts(redshifts)
{    
    m_spectrum.GetSpectralAxis().ClampLambdaRange(lambdaRange, m_lambdaRange );
}

/**
 * \brief this function estimates the likelihood_cstLog term withing the
 * wavelength range
 **/
Float64 COperatorTemplateFittingBase::EstimateLikelihoodCstLog(const CSpectrum &spectrum, 
                                                              const TFloat64Range &lambdaRange) 
{
    const CSpectrumSpectralAxis &spcSpectralAxis = spectrum.GetSpectralAxis();
    const TFloat64List &error = spectrum.GetFluxAxis().GetError().GetSamplesVector();;

    Int32 numDevs = 0;
    Float64 cstLog = 0.0;
    Float64 sumLogNoise = 0.0;

    Float64 imin = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetBegin());
    Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetEnd());
    for (UInt32 j = imin; j < imax; j++)
    {
        numDevs++;
        sumLogNoise += log(error[j]);
    }
    cstLog = -numDevs * 0.5 * log(2 * M_PI) - sumLogNoise;
    return cstLog;
}

std::shared_ptr<CModelSpectrumResult>   COperatorTemplateFittingBase::ComputeSpectrumModel( 
                                                            const std::shared_ptr<const CTemplate> & tpl,
                                                            Float64 redshift,
                                                            Float64 EbmvCoeff,
                                                            Int32 meiksinIdx,
                                                            Float64 amplitude,
                                                            std::string opt_interp,
                                                            const Float64 overlapThreshold)
{
    Log.LogDetail("  Operator-COperatorTemplateFitting: building spectrum model templateFitting for candidate Zcand=%f", redshift);
    
    Float64 overlapRate = 0.0;
    TFloat64Range currentRange;
    RebinTemplate(  tpl,
                    redshift, 
                    opt_interp,
                    currentRange,
                    overlapRate,
                    overlapThreshold);

    const TAxisSampleList & Xspc = m_spcSpectralAxis_restframe.GetSamplesVector();
    
    if ((EbmvCoeff>0.) || (meiksinIdx>-1)){
        m_templateRebined_bf.InitIsmIgmConfig(currentRange, redshift, tpl->m_ismCorrectionCalzetti, tpl->m_igmCorrectionMeiksin);
    }
    
    if (EbmvCoeff>0.)
    {
        if (m_templateRebined_bf.CalzettiInitFailed())
        {
            throw GlobalException(INTERNAL_ERROR,"no calzetti initialization for template");
        }
        Int32 idxEbmv = -1;
        idxEbmv = m_templateRebined_bf.m_ismCorrectionCalzetti->GetEbmvIndex(EbmvCoeff); 

        if (idxEbmv!=-1)
            m_templateRebined_bf.ApplyDustCoeff(idxEbmv);
    }

    if(meiksinIdx>-1)
    {
        if (m_templateRebined_bf.MeiksinInitFailed())
        {
            throw GlobalException(INTERNAL_ERROR,"no meiksin initialization for template");
        }
        m_templateRebined_bf.ApplyMeiksinCoeff(meiksinIdx);
    } 
    m_templateRebined_bf.ScaleFluxAxis(amplitude);
    //shift the spectralaxis to sync with the spectrum lambdaAxis
    const CSpectrumFluxAxis & modelflux =  m_templateRebined_bf.GetFluxAxis();
    CSpectrumSpectralAxis modelwav = m_templateRebined_bf.GetSpectralAxis(); // needs a copy to be shifted
    modelwav.ShiftByWaveLength((1.0+redshift), CSpectrumSpectralAxis::nShiftForward );

    return std::make_shared<CModelSpectrumResult>(CSpectrum(std::move(modelwav), modelflux));
}

void  COperatorTemplateFittingBase::RebinTemplate(  const std::shared_ptr<const CTemplate>& tpl, 
                                                    Float64 redshift,
                                                    const std::string & opt_interp,
                                                    TFloat64Range& currentRange,
                                                    Float64& overlapRate,
                                                    const Float64 overlapThreshold)
{
    Float64 onePlusRedshift = 1.0 + redshift;

    //shift lambdaRange backward to be in restframe
    TFloat64Range spcLambdaRange_restframe;
    TFloat64Range lambdaRange_restframe( m_lambdaRange.GetBegin() / onePlusRedshift,
                                         m_lambdaRange.GetEnd() / onePlusRedshift );

    //redshift in restframe the tgtSpectralAxis, i.e., division by (1+Z)
    m_spcSpectralAxis_restframe.ShiftByWaveLength(m_spectrum.GetSpectralAxis(), onePlusRedshift, CSpectrumSpectralAxis::nShiftBackward);
    m_spcSpectralAxis_restframe.ClampLambdaRange( lambdaRange_restframe, spcLambdaRange_restframe );
                                         
    // Compute clamped lambda range over template in restframe
    TFloat64Range tplLambdaRange;
    const CSpectrumSpectralAxis& tplSpectralAxis = tpl->GetSpectralAxis();
    tplSpectralAxis.ClampLambdaRange( lambdaRange_restframe, tplLambdaRange );
    // Compute the intersected range
    TFloat64Range intersectedLambdaRange( 0.0, 0.0 );
    TFloat64Range::Intersect( tplLambdaRange, spcLambdaRange_restframe, intersectedLambdaRange );
    Bool b = tpl->Rebin( intersectedLambdaRange, m_spcSpectralAxis_restframe, m_templateRebined_bf, m_mskRebined_bf, opt_interp);   
    if(!b) throw GlobalException(INTERNAL_ERROR,"COperatorTemplateFittingBase::RebinTemplate: error in rebinning tpl");
    //overlapRate
    overlapRate = m_spcSpectralAxis_restframe.IntersectMaskAndComputeOverlapRate( lambdaRange_restframe, m_mskRebined_bf );

    // Check for overlap rate
    if( overlapRate < overlapThreshold || overlapRate<=0.0 )
    {
        //status = nStatus_NoOverlap; 
        throw GlobalException(OVERLAPRATE_NOTACCEPTABLE,Formatter()<<"overlaprate of "<<overlapRate);
    }

    TFloat64Range logIntersectedLambdaRange( log( intersectedLambdaRange.GetBegin() ), log( intersectedLambdaRange.GetEnd() ) );
    //the spectral axis should be in the same scale
    currentRange = logIntersectedLambdaRange;
    if( m_spcSpectralAxis_restframe.IsInLinearScale() != tplSpectralAxis.IsInLinearScale() )
    {
        //status = nStatus_DataError;
        throw GlobalException(INTERNAL_ERROR,"COperatorTemplateFittingBase::RebinTemplate: data and model not in the same scale (lin/log)");
    }
    if(m_spcSpectralAxis_restframe.IsInLinearScale()){
        currentRange = intersectedLambdaRange;
    }
    return;
}

//get z at which igm starts given that LyA starts at lbda_rest=1216
Float64 COperatorTemplateFittingBase::GetIGMStartingRedshiftValue(Float64 spcLbda0)
{
    Float64 lbdarest_lya = 1216.;
    return spcLbda0/lbdarest_lya-1; //check the rounding thing     
}