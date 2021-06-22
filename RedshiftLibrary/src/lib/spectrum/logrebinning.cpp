#include <RedshiftLibrary/spectrum/logrebinning.h>
#include <RedshiftLibrary/log/log.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <gsl/gsl_fit.h>
#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;

using namespace NSEpic;
using namespace std;


void CSpectrumLogRebinning::RebinInputs(CInputContext& inputContext)
{
    // we have a probem here, cause only considering redshift range of galaxies
    //if we want to rebin stars we should call again computlogstep with stars redshiftrange!! same for qso
    std::string errorRebinMethod = "rebinVariance";//rebin error axis as well

    TFloat64Range   redshiftRange = inputContext.GetParameterStore()->Get<TFloat64Range>("galaxy.redshiftrange");
    Float64         redshiftStep = inputContext.GetParameterStore()->Get<Float64>( "galaxy.redshiftstep" );
    UInt32          SSratio = 1;
    if(inputContext.GetParameterStore()->Has<UInt32>( "galaxy.linemodelsolve.linemodel.firstpass.largegridstepratio"))
        SSratio = inputContext.GetParameterStore()->Get<UInt32>( "galaxy.linemodelsolve.linemodel.firstpass.largegridstepratio");

    std::shared_ptr<CSpectrum> spc;

    if(inputContext.GetSpectrum()->GetSpectralAxis().IsLogSampled()){
        spc = make_shared<CSpectrum>(inputContext.GetSpectrum()->GetName());
        CSpectrumSpectralAxis  spcWav = inputContext.GetSpectrum()->GetSpectralAxis();
        spcWav.RecomputePreciseLoglambda(); // in case input spectral values have been rounded
        spc->SetSpectralAndFluxAxes(std::move(spcWav), inputContext.GetSpectrum()->GetFluxAxis());
    }else{
        spc = inputContext.GetSpectrum();
    }

    SetupRebinning(*spc, 
                   inputContext.m_lambdaRange, 
                   redshiftStep, 
                   redshiftRange,
                   SSratio);
                   
    if(inputContext.GetSpectrum()->GetSpectralAxis().IsLogSampled()){
        inputContext.SetRebinnedSpectrum(spc);
    }else
        inputContext.SetRebinnedSpectrum(LoglambdaRebinSpectrum(inputContext.GetSpectrum(), errorRebinMethod));        
    
    inputContext.m_redshiftRangeFFT = m_zrange;
    inputContext.m_redshiftStepFFT = m_logGridStep;
    //rebin templates using previously identified parameters,
    // rebin only if rebinning parameters are different from previously used ones
    for(std::string s : inputContext.GetTemplateCatalog()->GetCategoryList()) //should retstrict to galaxy templates for now... (else depends on the fftprocessing by object type)
    { 
        //rebin only galaxy templates
        if(s!="galaxy")
            continue;
        // check existence of already  & correctly logsampled templates
        inputContext.GetTemplateCatalog()->m_logsampling = true;
        TTemplateRefList  TplList = inputContext.GetTemplateCatalog()->GetTemplate(TStringList{s});
        inputContext.GetTemplateCatalog()->m_logsampling = false;
        if (!TplList.empty())
        {            
            for (auto tpl : TplList)
            {   
                Bool needrebinning = false;
                if( !tpl->GetSpectralAxis().IsLogSampled(m_logGridStep) )
                    needrebinning = true;
                if (m_lambdaRange_tpl.GetBegin() < tpl->GetSpectralAxis()[0])
                    needrebinning = true;
                if (m_lambdaRange_tpl.GetEnd() > tpl->GetSpectralAxis()[tpl->GetSampleCount() - 1])
                    needrebinning = true;
                if (!CheckTemplateAlignment(tpl))
                    needrebinning = true;

                if (needrebinning)
                {
                    Log.LogDetail(" CInputContext::RebinInputs: need to rebin again the template: %s", tpl->GetName().c_str());
                    std::shared_ptr<const CTemplate> input_tpl = inputContext.GetTemplateCatalog()->GetTemplateByName(TStringList{s}, tpl->GetName());  
                    tpl = LoglambdaRebinTemplate(input_tpl); // assigin the tpl pointer to a new rebined template
                }
            } 
            continue; // next category
        }

        // no rebined templates in the category: rebin all templates
        TplList = inputContext.GetTemplateCatalog()->GetTemplate(TStringList{s});
        inputContext.GetTemplateCatalog()->m_logsampling = true;
        for (auto tpl : TplList)
        {   
            std::shared_ptr<CTemplate> rebinnedTpl = LoglambdaRebinTemplate(tpl);
            inputContext.GetTemplateCatalog()->Add(rebinnedTpl);
        }
        inputContext.GetTemplateCatalog()->m_logsampling = false;
    }  
}

/**
 * Brief: Get loglambdastep and update the zrange accordingly
 * Below code relies on the fact that both loglambda grid and the log(Redshift+1) grid follows the same arithmetic progession with a common step
 * if spectrum is already rebinned, then it imposes the rebinning and the creation of zGrid 
 * Otherwise, it's the input zrange that decides on the rebinning param. 
 * Note : If we want the log step to be log(1.0+redshiftstep) then spreadOverlog should 
  be modified to use 1+redshiftstep for the common ratio (or construct the grid using arithmetic log progression).
*/
void CSpectrumLogRebinning::SetupRebinning( CSpectrum &spectrum,
                                            const TFloat64Range &lambdaRange, 
                                            Float64 zInputStep,
                                            const TFloat64Range & zInputRange,
                                            UInt32 SubSamplingRatio)
{   
    TFloat64Range lambdaRange_ref;

    if(spectrum.GetSpectralAxis().IsLogSampled() )
    { 
        m_logGridStep = spectrum.GetSpectralAxis().GetlogGridStep();
 
        // compute reference lambda range
        // (the effective lambda range of log-sampled spectrum when initial spectrum overlaps lambdaRange 
        // assuming all spectra are aligned (this to avoid rebining the template several times)
        TLambdaRange lambda_range_spc = spectrum.GetLambdaRange();
        Float64 loglambda_start_spc = log(lambda_range_spc.GetBegin());
        Float64 loglambda_end_spc = log(lambda_range_spc.GetEnd());
        Float64 loglambda_start_ref = log(lambdaRange.GetBegin());
        Float64 loglambda_end_ref = log(lambdaRange.GetEnd());
        loglambda_start_ref = loglambda_start_spc + ceil((loglambda_start_ref - loglambda_start_spc)/m_logGridStep)*m_logGridStep; 
        loglambda_end_ref = loglambda_end_spc + floor((loglambda_end_ref - loglambda_end_spc)/m_logGridStep)*m_logGridStep; 
        lambdaRange_ref = TFloat64Range(exp(loglambda_start_ref), exp(loglambda_end_ref));// this is the effective lambda range, to be used in InferTemplateRebinningSetup

        // no need to compute m_lambdaRange_spc (used only for resampling input spectra)

    }else{
        m_logGridStep = zInputStep;

        // compute reference lambda range (the effective lambda range of rebinned spectrum when initial spectrum overlaps lambdaRange)
        Int32 loglambda_count_ref;
        Float64 loglambda_start_ref = log(lambdaRange.GetBegin());
        Float64 loglambda_end_ref = log(lambdaRange.GetEnd());
        loglambda_count_ref =  Int32( floor(( loglambda_end_ref - loglambda_start_ref)/m_logGridStep) );// we should sample inside range hence floor
        loglambda_end_ref = loglambda_start_ref + loglambda_count_ref*m_logGridStep;
        lambdaRange_ref = TFloat64Range(lambdaRange.GetBegin(), exp(loglambda_end_ref));// this is the effective lambda range, to be used in InferTemplateRebinningSetup

        // compute rebinned spectrum lambda range m_lambdaRange_spc (ie clamp on reference grid)
        // to be passed to computeTargetLogSpectralAxis in LogLambdaRebinSpectrum
        spectrum.GetSpectralAxis().ClampLambdaRange(lambdaRange_ref, m_lambdaRange_spc);
        Float64 loglambda_start_spc = log(m_lambdaRange_spc.GetBegin());
        Float64 loglambda_end_spc = log(m_lambdaRange_spc.GetEnd());
        if (m_lambdaRange_spc.GetBegin() > lambdaRange_ref.GetBegin()) 
        {
            loglambda_start_spc = loglambda_start_ref + ceil((loglambda_start_spc - loglambda_start_ref)/m_logGridStep)*m_logGridStep; // ceil to be bigger or equal the first sample
            m_lambdaRange_spc.SetBegin(exp(loglambda_start_spc));
        }
        if (m_lambdaRange_spc.GetEnd() < lambdaRange_ref.GetEnd())
        {
            loglambda_end_spc = loglambda_end_ref - ceil((loglambda_end_ref - loglambda_end_spc)/m_logGridStep)*m_logGridStep; // ceil to be less or equal the last sample
            m_lambdaRange_spc.SetEnd(exp(loglambda_end_spc));
        }
        Float64 count_ = (loglambda_end_spc - loglambda_start_spc)/m_logGridStep; // should integer at numerical precision...
        m_loglambda_count_spc = round(count_) + 1; 

        if(m_loglambda_count_spc<2){
            Log.LogError("   Operator-TemplateFittingLog: logGridCount = %d <2", m_loglambda_count_spc);
            throw runtime_error("   Operator-TemplateFittingLog: logGridCount <2. Abort");
        }
    }
    Log.LogDetail("  Log-Rebin: logGridStep = %f", m_logGridStep);

    // compute the effective zrange of the new redshift grid
    // set the min to the initial min
    // set the max to an interger number of log(z+1) steps 
    Float64 zmin_new = zInputRange.GetBegin(),
            zmax_new = zInputRange.GetEnd();
    {
        Float64 log_zmin_new_p1 = log(zmin_new + 1.);
        Float64 log_zmax_new_p1 = log(zmax_new + 1.);
        Int32 nb_z = Int32( ceil((log_zmax_new_p1 - log_zmin_new_p1)/m_logGridStep/SubSamplingRatio) );
        zmax_new = exp(log_zmin_new_p1 + nb_z*m_logGridStep*SubSamplingRatio) - 1.;
    }
    m_zrange = TFloat64Range(zmin_new, zmax_new); //updating zrange based on the new logstep 

    InferTemplateRebinningSetup(lambdaRange_ref);

    return;
}

/**
*  Rebin the spectrum with the calculated logGridStep if spectrum not already rebinned:
*  step1: construct the spectralAxis
*  step2: do the rebin
*/
std::shared_ptr< CSpectrum> CSpectrumLogRebinning::LoglambdaRebinSpectrum( std::shared_ptr<const CSpectrum> spectrum, 
                                                                                std::string errorRebinMethod) const
{ 
    Log.LogInfo("  Operator-TemplateFittingLog: Log-regular lambda resampling START");

    //prepare return rebinned vector  
    auto spectrumRebinedLog = make_shared<CSpectrum>(spectrum->GetName());
    CMask mskRebinedLog;

    CSpectrumSpectralAxis targetSpectralAxis =  computeTargetLogSpectralAxis(m_lambdaRange_spc, m_loglambda_count_spc);

    TFloat64Range spcLbdaRange(targetSpectralAxis[0] - 0.5 * m_logGridStep,
                            targetSpectralAxis[m_loglambda_count_spc-1] + 0.5 * m_logGridStep);

    // rebin the spectrum
    Bool ret = spectrum->Rebin(spcLbdaRange, targetSpectralAxis,
                    *spectrumRebinedLog,
                    mskRebinedLog, m_rebinMethod, errorRebinMethod);
    if(!ret){
        throw runtime_error("Cant rebin spectrum");
    }
    //spectrumRebinedLog->GetSpectralAxis().SetLogScale(); //we need to do a convertScale before setting scale
    spectrumRebinedLog->GetSpectralAxis().IsLogSampled(m_logGridStep);//double make sure that sampling is well done

    //the rebinned spectrum and we change logscale
    Log.LogDetail("  Operator-TemplateFittingLog: Log-regular lambda resampling FINISHED");

    return spectrumRebinedLog;  
}

/*
    Aims at computing the template lambda range once for all template catalog
*/
void CSpectrumLogRebinning::InferTemplateRebinningSetup(const TFloat64Range & lambdaRange_ref)
{
    Float64 loglbdamin = log( lambdaRange_ref.GetBegin()/(1.0 + m_zrange.GetEnd()));
    Float64 loglbdamax = log( lambdaRange_ref.GetEnd()/(1.0 + m_zrange.GetBegin()));
    UInt32 _round = std::round((loglbdamax - loglbdamin)/m_logGridStep)+1;
    Float64 _neat = (loglbdamax - loglbdamin)/m_logGridStep + 1;//we expect to get an int value with no need to any rounding
    if( std::abs(_round - _neat)>1E-8)
    {
        Log.LogError("Problem in logrebinning setup");
        throw runtime_error("Problem in logrebinning setup!");
    }
    m_loglambda_count_tpl = std::round((loglbdamax - loglbdamin)/m_logGridStep) + 1;

    Float64 tgt_loglbdamax = loglbdamax;
    Float64 tgt_loglbdamin = loglbdamax - (m_loglambda_count_tpl-1)* m_logGridStep;
    m_lambdaRange_tpl = TFloat64Range(exp(tgt_loglbdamin), exp(tgt_loglbdamax));
    Log.LogDetail("  Operator-TemplateFittingLog: Log-Rebin: tpl raw loglbdamin=%f : raw loglbdamax=%f", loglbdamin, loglbdamax);
	Log.LogDetail("  Operator-TemplateFittingLog: zmin_new = %f, tpl->lbdamax = %f", m_zrange.GetBegin(), exp(loglbdamax));
	Log.LogDetail("  Operator-TemplateFittingLog: zmax_new = %f, tpl->lbdamin = %f", m_zrange.GetEnd(), exp(loglbdamin)); 
   return;
}
/** 
 *                          
 * Brief: Log Rebin the template spectral axis
 * Important: below function considered that we have already constructed the new log spectralAxis of our spectrum
 * Step1: align the rebined spectrum max lambda at min redshift. 
 * Step2: get grid count for target template and update size accordingly
 * Step3: construct target loglambda axis for the template and check borders
 * Step4: rebin the template --> rebinned flux is saved in templateRebinedLog
**/
std::shared_ptr<CTemplate> CSpectrumLogRebinning::LoglambdaRebinTemplate(std::shared_ptr<const CTemplate> tpl) const
{
    Log.LogInfo("  Operator-TemplateFittingLog: Log-regular lambda resampling START for template %s", tpl->GetName().c_str());
    // check template coverage is enough for zrange and spectrum coverage
    Bool overlapFull = true;
    if (m_lambdaRange_tpl.GetBegin() < tpl->GetSpectralAxis()[0])
        overlapFull = false;
    if (m_lambdaRange_tpl.GetEnd() > tpl->GetSpectralAxis()[tpl->GetSampleCount() - 1])
        overlapFull = false;
    if (!overlapFull)
    {
        Log.LogError(" CSpectrumLogRebinning::LoglambdaRebinTemplate: overlap found to be lower than 1.0 for the template %s", tpl->GetName().c_str());
        throw std::runtime_error("CInputContext::RebinInputs: overlap found to be lower than 1.0");
    }

    CSpectrumSpectralAxis targetSpectralAxis =  computeTargetLogSpectralAxis(m_lambdaRange_tpl,
                                                                            m_loglambda_count_tpl);

    if(targetSpectralAxis[m_loglambda_count_tpl-1]< targetSpectralAxis[0])
        throw runtime_error(" Last elements of the target spectral axis are not valid. Template count is not well computed due to exp/conversions");

    auto templateRebinedLog = make_shared<CTemplate>(tpl->GetName(), tpl->GetCategory());

    CMask mskRebinedLog;
    
    TFloat64Range tplLbdaRange(targetSpectralAxis[0] - 0.5 * m_logGridStep,
                        targetSpectralAxis[m_loglambda_count_tpl-1] + 0.5 * m_logGridStep);   
                                                               
    Bool ret = tpl->Rebin(  tplLbdaRange,
                targetSpectralAxis,
                *templateRebinedLog,
                mskRebinedLog);
    if(!ret){
        throw runtime_error("Cant rebin template");
    }

    templateRebinedLog->GetSpectralAxis().IsLogSampled(m_logGridStep);

    return templateRebinedLog;
}

CSpectrumSpectralAxis CSpectrumLogRebinning::computeTargetLogSpectralAxis(const TFloat64Range & lambdarange, UInt32 count) const
{//spreadoverlog expects m_Begin to be non-log value
    TFloat64List axis = lambdarange.SpreadOverLog(m_logGridStep);
    if(axis.size()!= count)
    {
        Log.LogError("  CSpectrumLogRebinning::computeTargetLogSpectralAxis: computed axis has not expected samples number");
        throw runtime_error("  CSpectrumLogRebinning::computeTargetLogSpectralAxis: computed axis has not expected samples number");
    }
    CSpectrumSpectralAxis targetSpectralAxis(std::move(axis));
    return targetSpectralAxis;
}

Bool CSpectrumLogRebinning::CheckTemplateAlignment(const std::shared_ptr<const CTemplate> &tpl) const 
{
    const TAxisSampleList &w = tpl->GetSpectralAxis().GetSamplesVector();
    const Float64 &lstart = m_lambdaRange_tpl.GetBegin();
    Int32 idx=CIndexing<Float64>::getCloserIndex(w, lstart);
    return (lstart-w[idx])/w[idx]<=2E-7;
}
