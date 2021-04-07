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
    
    //non clamped lambdaRange: to be clamped depending on used spectra
    m_lambdaRange = m_ParameterStore->Get<TFloat64Range>("lambdarange");  
    
    //TODO: It could be relevant to add a new function, ::hasFFTProcessing containing the below code, to the paramStore
    //The drawback of the below code is that we are violating object and method scopes by trying to read "deep" info
    //Another option could be to move fftprocessing to the object scope (vs object.method scope)
    bool fft_processing = false, fft_processing_qso = false, fft_processing_star = false, fft_processinglmContinuum = false;
    if(m_ParameterStore->Has<std::string>( "qso.templatefittingsolve.fftprocessing"))
        fft_processing_qso = m_ParameterStore->Get<std::string>( "qso.templatefittingsolve.fftprocessing") == "yes";
    if(m_ParameterStore->Has<std::string>( "star.templatefittingsolve.fftprocessing"))
        fft_processing_star = m_ParameterStore->Get<std::string>( "star.templatefittingsolve.fftprocessing") == "yes";
    if(m_ParameterStore->Has<std::string>("galaxy.templatefittingsolve.fftprocessing"))
        fft_processing = m_ParameterStore->Get<std::string>("galaxy.templatefittingsolve.fftprocessing") == "yes";
    if(m_ParameterStore->Has<std::string>("galaxy.linemodelsolve.linemodel.continuumfit.fftprocessing"))
        fft_processinglmContinuum = m_ParameterStore->Get<std::string>("galaxy.linemodelsolve.linemodel.continuumfit.fftprocessing")=="yes";
    /*if(m_ParameterStore->Has<std::string>("galaxy.linemodelsolve.linemodel.fftprocessing"))
        fft_processinglm = m_ParameterStore->Get<std::string>("galaxy.linemodelsolve.linemodel.fftprocessing")=="yes";//new param to decide if we should use fft for linemodel 
    */
    if(fft_processing_qso || fft_processing_star)
    {
        Log.LogError("FFT processing is not yet supported for stars or qso");
        throw std::runtime_error("FFT processing is not yet supported for star or qso");
    }
    m_use_LogLambaSpectrum = fft_processing || fft_processing_qso || fft_processing_star || fft_processinglmContinuum;
    
    /*
    std::string fft_processing, fft_processing_qso, fft_processinglmContinuum, fft_processinglm;
    m_ParameterStore->Get("qsosolve.templatefittingsolve.fftprocessing", fft_processing_qso, "yes" );  //check if it has changed
    m_ParameterStore->Get("galaxy.templatefittingsolve.fftprocessing", fft_processing, "no");
    m_ParameterStore->Get("galaxy.linemodelsolve.linemodel.continuumfit.fftprocessing", fft_processinglmContinuum, "yes");//this hasnt changed yet
    m_ParameterStore->Get("galaxy.linemodelsolve.linemodel.fftprocessing", fft_processinglm, "no");//new param to decide if we should use fft for linemodel 
    m_use_LogLambaSpectrum = fft_processing=="yes" || fft_processing_qso=="yes" || fft_processinglmContinuum=="yes" ||fft_processinglm=="yes" ;
    */
    if(m_use_LogLambaSpectrum)
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

    TFloat64Range   m_redshiftRange = m_ParameterStore->Get<TFloat64Range>("redshiftrange");
    Float64         m_redshiftStep = m_ParameterStore->Get<Float64>( "redshiftstep" );
   
    // we have a probem here, cause only considering redshift range of galaxies
    //if we want to rebin stars we should call again computlogstep with stars redshiftrange!! same for qso
    CSpectrumLogRebinning logReb;
    std::string errorRebinMethod = "rebinVariance";//rebin error axis as well

    logReb.SetupRebinning(*m_Spectrum, m_lambdaRange, m_redshiftStep, m_redshiftRange);
    if(m_Spectrum->GetSpectralAxis().IsLogSampled()){
        m_rebinnedSpectrum = m_Spectrum;
    }else
        m_rebinnedSpectrum = logReb.LoglambdaRebinSpectrum(m_Spectrum, errorRebinMethod);        
    
    m_redshiftRangeFFT = logReb.m_zrange;
    m_redshiftStepFFT = logReb.m_logGridStep;
    //rebin templates using previously identified parameters,
    //TODO: rebin only if parameters to use are different from previously used params
    Float64 margin = 1E-8;
    for(std::string s : m_TemplateCatalog->GetCategoryList()) //should retstrct to galaxy templates for now... (else depends on the fftprocessing by object type)
    { 
        TTemplateRefList  TplList = m_TemplateCatalog->GetTemplate(TStringList{s});
        for (auto tpl : TplList)
        {   
            if( tpl->GetSpectralAxis().IsLogSampled(logReb.m_logGridStep) )
            {
                break;
            } 
            Bool overlapFull = true;
            if (logReb.m_lambdaRange_tpl.GetBegin() < tpl->GetSpectralAxis()[0])
                overlapFull = false;
            if (logReb.m_lambdaRange_tpl.GetEnd() > tpl->GetSpectralAxis()[tpl->GetSampleCount() - 1])
                overlapFull = false;
            if (!overlapFull)
            {
                Log.LogError("Operator-TemplateFittingLog: overlap found to be lower than 1.0 for this redshift range.");
                throw std::runtime_error("Operator-TemplateFittingLog: overlap found to be lower than 1.0 for this redshift range.");
            }
            std::shared_ptr<CTemplate> rebinnedTpl = logReb.LoglambdaRebinTemplate(tpl);
            m_TemplateCatalog->Add(rebinnedTpl, "log");
        } 
    }  
    m_TemplateCatalog->SetScale("log"); // this is wrong , should probably be done later.
    m_LogRebinningCompleted = 1;
}


