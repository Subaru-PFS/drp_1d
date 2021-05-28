#include <RedshiftLibrary/operator/templatefittingBase.h>
#include <RedshiftLibrary/operator/modelspectrumresult.h>

using namespace NSEpic;
using namespace std;

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

Int32   COperatorTemplateFittingBase::ComputeSpectrumModel( const CSpectrum& spectrum,
                                                            const CTemplate& tpl,
                                                            Float64 redshift,
                                                            Float64 EbmvCoeff,
                                                            Int32 meiksinIdx,
                                                            Float64 amplitude,
                                                            std::string opt_interp,
                                                            const TFloat64Range& lambdaRange,
                                                            Float64 overlapThreshold,
                                                            std::shared_ptr<CModelSpectrumResult> & spcPtr)
{
    Log.LogDetail("  Operator-COperatorTemplateFitting: building spectrum model templateFitting for candidate Zcand=%f", redshift);
    EStatus status;
    TFloat64Range currentRange;
    Float64 overlapRate = 0.0;

    Int32 ret = RebinTemplate(  spectrum, 
                                tpl,
                                redshift, 
                                lambdaRange,
                                opt_interp,
                                currentRange,
                                overlapRate,
                                overlapThreshold);
    if( ret == -1 ){
        status = nStatus_NoOverlap; 
        return -1;
    }
    if( ret == -2 ){
        status = nStatus_DataError;
        return -1;
    }
    const TAxisSampleList & Xspc = m_spcSpectralAxis_restframe.GetSamplesVector();
    
    if ((EbmvCoeff>0.) || (meiksinIdx>-1)){
        m_templateRebined_bf.InitIsmIgmConfig(tpl.m_ismCorrectionCalzetti, tpl.m_igmCorrectionMeiksin);
        m_templateRebined_bf.SetIsmIgmLambdaRange(currentRange);
    }
    
    if (EbmvCoeff>0.)
    {
        if (m_templateRebined_bf.CalzettiInitFailed())
        {
            Log.LogError("  Operator-TemplateFitting: asked model with Dust extinction with no calzetti calib. file loaded in template" );
            return -1;
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
            Log.LogError("  Operator-TemplateFitting: asked model with IGM extinction with no Meikin calib. file loaded in template" );
            return -1;
        }
        Bool igmCorrectionAppliedOnce = m_templateRebined_bf.ApplyMeiksinCoeff(meiksinIdx, redshift);
    } 
    m_templateRebined_bf.ScaleFluxAxis(amplitude);
    //shift the spectralaxis to sync with the spectrum lambdaAxis
    const CSpectrumFluxAxis & modelflux =  m_templateRebined_bf.GetFluxAxis();
    CSpectrumSpectralAxis modelwav = m_templateRebined_bf.GetSpectralAxis(); // needs a copy to be shifted
    modelwav.ShiftByWaveLength((1.0+redshift), CSpectrumSpectralAxis::nShiftForward ) ;
    spcPtr = std::make_shared<CModelSpectrumResult>(CSpectrum(std::move(modelwav), modelflux));
    return 0;
}

Int32  COperatorTemplateFittingBase::RebinTemplate( const CSpectrum& spectrum,
                                                    const CTemplate& tpl, 
                                                    Float64 redshift,
                                                    const TFloat64Range& lambdaRange,
                                                    std::string opt_interp,
                                                    //return variables
                                                    TFloat64Range& currentRange,
                                                    Float64& overlapRate,
                                                    Float64 overlapThreshold)// const
{
    Float64 onePlusRedshift = 1.0 + redshift;

    //shift lambdaRange backward to be in restframe
    TFloat64Range spcLambdaRange_restframe;
    TFloat64Range lambdaRange_restframe( lambdaRange.GetBegin() / onePlusRedshift,
                                         lambdaRange.GetEnd() / onePlusRedshift );

    //redshift in restframe the tgtSpectralAxis, i.e., division by (1+Z)
    m_spcSpectralAxis_restframe.ShiftByWaveLength(spectrum.GetSpectralAxis(), onePlusRedshift, CSpectrumSpectralAxis::nShiftBackward);
    m_spcSpectralAxis_restframe.ClampLambdaRange( lambdaRange_restframe, spcLambdaRange_restframe );
                                         
    // Compute clamped lambda range over template in restframe
    TFloat64Range tplLambdaRange;
    const CSpectrumSpectralAxis& tplSpectralAxis = tpl.GetSpectralAxis();
    tplSpectralAxis.ClampLambdaRange( lambdaRange_restframe, tplLambdaRange );
    // Compute the intersected range
    TFloat64Range intersectedLambdaRange( 0.0, 0.0 );
    TFloat64Range::Intersect( tplLambdaRange, spcLambdaRange_restframe, intersectedLambdaRange );
    //m_templateRebined_bf.ResetNoIsmIgmFlux();//reset 
    Bool b = tpl.Rebin( intersectedLambdaRange, m_spcSpectralAxis_restframe, m_templateRebined_bf, m_mskRebined_bf, opt_interp);   
    if(!b) throw runtime_error("problems rebinning tpl");
    //overlapRate
    overlapRate = m_spcSpectralAxis_restframe.IntersectMaskAndComputeOverlapRate( lambdaRange_restframe, m_mskRebined_bf );

    // Check for overlap rate
    if( overlapRate < overlapThreshold || overlapRate<=0.0 )
    {
        //status = nStatus_NoOverlap; 
        return -1 ;
    }

    TFloat64Range logIntersectedLambdaRange( log( intersectedLambdaRange.GetBegin() ), log( intersectedLambdaRange.GetEnd() ) );
    //the spectral axis should be in the same scale
    currentRange = logIntersectedLambdaRange;
    if( m_spcSpectralAxis_restframe.IsInLinearScale() != tplSpectralAxis.IsInLinearScale() )
    {
        Log.LogError( "    chisquare operator: data and model not in the same scale (lin/log) ! Aborting.");
        //status = nStatus_DataError;
        return -2;
    }
    if(m_spcSpectralAxis_restframe.IsInLinearScale()){
        currentRange = intersectedLambdaRange;
    }
    return 0;
}
