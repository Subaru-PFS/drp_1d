#include <RedshiftLibrary/processflow/inputcontext.h>
#include <RedshiftLibrary/processflow/parameterstore.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/logrebinning.h>
using namespace NSEpic;

CInputContext::CInputContext(std::shared_ptr<CSpectrum> spc,
                             std::shared_ptr<CTemplateCatalog> tmplCatalog,
                             std::shared_ptr<CRayCatalog> rayCatalog,
                             std::shared_ptr<CParameterStore> paramStore):
  m_Spectrum(std::move(spc)),
  m_TemplateCatalog(std::move(tmplCatalog)),
  m_RayCatalog(std::move(rayCatalog)),
  m_ParameterStore(std::move(paramStore))
{
    m_Spectrum->InitSpectrum(*m_ParameterStore);
    
    /*//todo use ::has prior to reading
    bool fft_processing_qso = fm_ParameterStore->Get<std::string>( "qso.templatefittingsolve.fftprocessing") == "yes";
    bool fft_processing = m_ParameterStore->Get<std::string>("galaxy.templatefittingsolve.fftprocessing") == "yes";
    bool fft_processinglmContinuum = m_ParameterStore->Get<std::string>("galaxy.linemodelsolve.linemodel.continuumfit.fftprocessing")=="yes";
    bool fft_processinglm = m_ParameterStore->Get<std::string>("galaxy.linemodelsolve.linemodel.fftprocessing")=="yes";//new param to decide if we should use fft for linemodel 
    if( fft_processing || fft_processing_qso || fft_processinglmContinuum ||fft_processinglm)
        RebinInputs();
    */
    std::string fft_processing, fft_processing_qso, fft_processinglmContinuum, fft_processinglm;
    m_ParameterStore->Get("qsosolve.templatefittingsolve.fftprocessing", fft_processing_qso, "yes" );  //check if it has changed
    m_ParameterStore->Get("galaxy.templatefittingsolve.fftprocessing", fft_processing, "no");
    m_ParameterStore->Get("galaxy.linemodelsolve.linemodel.continuumfit.fftprocessing", fft_processinglmContinuum, "yes");//this hasnt changed yet
    m_ParameterStore->Get("galaxy.linemodelsolve.linemodel.fftprocessing", fft_processinglm, "no");//new param to decide if we should use fft for linemodel 
    if( fft_processing=="yes" || fft_processing_qso=="yes" || fft_processinglmContinuum=="yes" ||fft_processinglm=="yes")
        RebinInputs();

    // Calzetti ISM & Meiksin IGM initialization, for both rebinned and original templates
    std::string calibrationPath =  m_ParameterStore->Get<std::string>( "calibrationDir");  
    m_TemplateCatalog->InitIsmIgm(calibrationPath);

}
/*
Two cases exist:
1. input spectrum is linear-sampled
2. input spectrum is log-sampled

C1 _Case1: 
1. To log-sample it: 
** determine log-step based on zgrid step (from parameter.json)
** lambda reference = lambdaRange.GetBegin() , thus independent from spectrum spectral axis cause anyway we are gona rebin the spectrum
** use this log-step and lambda reference for template rebin as well

C2 _Case2:
1. Check that the spectrum is well rebinned, using @logReb::CheckLoglambdaRebinSpectrum@
** determine logLambda step and log(Z+1) step and the reference lambda value
** use these params to rebin templates. 
** in the case of running bunches, check that all spectra are rebinned the same way, otherwise we have to rebin 
** the template catalog differently for each spectrum. 

Template rebinning:
* for the same reference lambda value and logLambda step, rebin ONCE the template catalogs
* If any param changes, redo the template rebinning:
** Add getters for templates as follows: @GetTemplate(..., "log", logstep, ref lambda)@
** the getter: 1)returns rebinnedTemplates OR 2) redo the template rebinning

Rebinning parameters for _Case2 should be extracted from m_Spectrum object, thus the client has the responsibility to add these info to each spectrum
*/

void CInputContext::RebinInputs() 
{
    TFloat64Range lambdaRange, spclambdaRange, redshiftRange;
    Float64       redshiftStep;
    m_ParameterStore->Get( "lambdarange", lambdaRange );
    m_ParameterStore->Get( "redshiftrange", redshiftRange );
    // we have a probem here, cause only considering redshift range of galaxies
    //if we want to rebin stars we should call again computlogstep with stars redshiftrange!! same for qso
    m_ParameterStore->Get( "redshiftstep", redshiftStep );
 
    CSpectrumLogRebinning logReb;
    std::string errorRebinMethod = "rebinVariance";//rebin error axis as well
    m_Spectrum->GetSpectralAxis().ClampLambdaRange( lambdaRange, spclambdaRange );

    m_Spectrum->GetSpectralAxis().IsLogSampled();
    logReb.Computelogstep(*m_Spectrum, lambdaRange, redshiftStep, redshiftRange);
    m_rebinnedSpectrum = logReb.LoglambdaRebinSpectrum(m_Spectrum, spclambdaRange, errorRebinMethod);        

    m_zrange = logReb.m_zrange;
    m_logGridStep = logReb.m_logGridStep;
    //rebin templates using previously identified parameters,
    //TODO: rebin only if parameters to use are different from previously used params
    Float64 margin = 1E-8;
    for(std::string s : m_TemplateCatalog->GetCategoryList())
    { 
        TTemplateRefList  TplList = m_TemplateCatalog->GetTemplate(TStringList{s});
        for (auto tpl : TplList)
        {   
            if( tpl->GetSpectralAxis().IsLogSampled(logReb.m_logGridStep) )//still missing another condition?!
            {
                break;
            }
            
            Bool overlapFull = true;
            if (spclambdaRange.GetBegin() < tpl->GetSpectralAxis()[0] * (1. + m_zrange.GetEnd()))
                overlapFull = false;
            if (spclambdaRange.GetEnd() > tpl->GetSpectralAxis()[tpl->GetSampleCount() - 1] * (1. + m_zrange.GetBegin()))
                overlapFull = false;
            if (!overlapFull)
            {
                Log.LogError("Operator-TemplateFittingLog: overlap found to be lower than 1.0 for this redshift range.");
                throw std::runtime_error("Operator-TemplateFittingLog: overlap found to be lower than 1.0 for this redshift range.");
            }
            std::shared_ptr<CTemplate> rebinnedTpl = logReb.LoglambdaRebinTemplate(tpl);
            m_TemplateCatalog->Add(rebinnedTpl, "log");
            logReb.verifyLogRebinningResults(m_rebinnedSpectrum->GetSpectralAxis(), rebinnedTpl->GetSpectralAxis());
        } 
    }  
    m_TemplateCatalog->SetScale("log"); 
    m_LogRebinningCompleted = 1;
}


