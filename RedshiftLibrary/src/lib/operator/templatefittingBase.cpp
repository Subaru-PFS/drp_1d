#include <RedshiftLibrary/operator/templatefittingBase.h>
#include <RedshiftLibrary/operator/modelspectrumresult.h>

using namespace NSEpic;
using namespace std;

COperatorTemplateFittingBase::COperatorTemplateFittingBase()
{

}

COperatorTemplateFittingBase::~COperatorTemplateFittingBase()
{

}

Int32   COperatorTemplateFittingBase::ComputeSpectrumModel(const CSpectrum& spectrum,
                                           const CTemplate& tpl,
                                           Float64 redshift,
                                           Float64 DustCoeff,
                                           Int32 meiksinIdx,
                                           Float64 amplitude,
                                           std::string opt_interp,
					                       std::string opt_extinction,
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
    m_templateRebined_bf.SetIsmIgmLambdaRange(currentRange);

    if ((DustCoeff>0.) || (meiksinIdx>0)){
        m_templateRebined_bf.InitIsmIgmConfig(tpl.m_ismCorrectionCalzetti, tpl.m_igmCorrectionMeiksin);
    }

    if (DustCoeff>0.)
    {
        if (m_templateRebined_bf.CalzettiInitFailed())
        {
            Log.LogError("  Operator-TemplateFitting: asked model with Dust extinction with no calzetti calib. file loaded in template" );
            return -1;
        }
        Int32 idxDust = -1;
        for(Int32 kDust=0; m_templateRebined_bf.m_ismCorrectionCalzetti->GetNPrecomputedDustCoeffs(); kDust++)
        {
            Float64 kDustCoeff= m_templateRebined_bf.m_ismCorrectionCalzetti->GetEbmvValue(kDust);
            if(DustCoeff==kDustCoeff)
            {
                idxDust = kDust;
                break;
            }
        }

        if (idxDust!=-1)
            m_templateRebined_bf.ApplyDustCoeff(idxDust);
    }

    if(opt_extinction == "yes")
    {
        if (m_templateRebined_bf.MeiksinInitFailed())
        {
            Log.LogError("  Operator-TemplateFitting: asked model with IGM extinction with no Meikin calib. file loaded in template" );
            return -1;
        }
        Bool igmCorrectionAppliedOnce = m_templateRebined_bf.ApplyMeiksinCoeff(meiksinIdx, redshift);
    } 
    m_templateRebined_bf.ScaleFluxAxis(amplitude);
    spcPtr = std::make_shared<CModelSpectrumResult>(m_templateRebined_bf);
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

    tpl.Rebin( intersectedLambdaRange, m_spcSpectralAxis_restframe, m_templateRebined_bf, m_mskRebined_bf, opt_interp);   

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
