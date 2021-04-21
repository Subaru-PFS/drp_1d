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

    TFloat64Range   redshiftRange = inputContext.m_ParameterStore->Get<TFloat64Range>("galaxy.redshiftrange");
    Float64         redshiftStep = inputContext.m_ParameterStore->Get<Float64>( "galaxy.redshiftstep" );

    SetupRebinning(*inputContext.m_Spectrum, 
                   inputContext.m_lambdaRange, 
                   redshiftStep, 
                   redshiftRange);
    if(inputContext.m_Spectrum->GetSpectralAxis().IsLogSampled()){
        inputContext.m_rebinnedSpectrum = inputContext.m_Spectrum;
    }else
        inputContext.m_rebinnedSpectrum = LoglambdaRebinSpectrum(inputContext.m_Spectrum, errorRebinMethod);        
    
    inputContext.m_redshiftRangeFFT = m_zrange;
    inputContext.m_redshiftStepFFT = m_logGridStep;
    //rebin templates using previously identified parameters,
    //TODO: rebin only if parameters to use are different from previously used params
    for(std::string s : inputContext.m_TemplateCatalog->GetCategoryList()) //should retstrict to galaxy templates for now... (else depends on the fftprocessing by object type)
    { 
        // check existence of already  & correctly logsampled templates
        inputContext.m_TemplateCatalog->m_logsampling = true;
        TTemplateRefList  TplList = inputContext.m_TemplateCatalog->GetTemplate(TStringList{s});
        inputContext.m_TemplateCatalog->m_logsampling = false;
        if (!TplList.empty())
        {            
            for (auto tpl : TplList)
            {   
                Bool needrebinning = false;
                if( tpl->GetSpectralAxis().IsLogSampled(m_logGridStep) )
                    needrebinning = true;
                if (m_lambdaRange_tpl.GetBegin() < tpl->GetSpectralAxis()[0])
                    needrebinning = true;
                if (m_lambdaRange_tpl.GetEnd() > tpl->GetSpectralAxis()[tpl->GetSampleCount() - 1])
                    needrebinning = true;
                
                if (needrebinning)
                {
                    Log.LogDetail(" CInputContext::RebinInputs: need to rebin again the template: %s", tpl->GetName().c_str());
                    std::shared_ptr<const CTemplate> input_tpl = inputContext.m_TemplateCatalog->GetTemplateByName(TStringList{s}, tpl->GetName());  
                    tpl = LoglambdaRebinTemplate(input_tpl); // assigin the tpl pointer to a new rebined template
                }
            } 
            break; // next category
        }

        // no rebined templates in the category: rebin all templates
        TplList = inputContext.m_TemplateCatalog->GetTemplate(TStringList{s});
        inputContext.m_TemplateCatalog->m_logsampling = true;
        for (auto tpl : TplList)
        {   
            std::shared_ptr<CTemplate> rebinnedTpl = LoglambdaRebinTemplate(tpl);
            inputContext.m_TemplateCatalog->Add(rebinnedTpl);
        }
        inputContext.m_TemplateCatalog->m_logsampling = false;
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
                                            const Float64 zInputStep,
                                            const TFloat64Range zInputRange)
{   
 
    if(spectrum.GetSpectralAxis().IsLogSampled() )
    { 
        spectrum.GetSpectralAxis().ClampLambdaRange(lambdaRange,m_lambdaRange_spc);
        m_lambdaRange_tpl = m_lambdaRange_spc;
        m_logGridStep = spectrum.GetSpectralAxis().GetlogGridStep();
    }else{
        m_logGridStep = zInputStep;

        // compute reference lambda range (the effective lambda range of rebinned spectrum when initial spectrum overlaps lambdaRange)
        Float64 loglambda_start_tpl = log(lambdaRange.GetBegin());
        Float64 loglambda_end_tpl = log(lambdaRange.GetEnd());
        m_loglambda_count_tpl =  Int32( floor(( loglambda_end_tpl - loglambda_start_tpl)/m_logGridStep) );// we should sample inside range hence floor
        loglambda_end_tpl = loglambda_start_tpl + m_loglambda_count_tpl*m_logGridStep;
        m_lambdaRange_tpl = TFloat64Range(lambdaRange.GetBegin(), exp(loglambda_end_tpl));// this is the effective lambda range, to be used in InferTemplateRebinningSetup

        // compute rebinned spectrum lambda range (ie clamp on reference grid)
        // to be passed to computeTargetLogSpectralAxis in LogLambdaRebinSpectrum
        spectrum.GetSpectralAxis().ClampLambdaRange(m_lambdaRange_tpl, m_lambdaRange_spc);
        Float64 loglambda_start_spc = log(m_lambdaRange_spc.GetBegin());
        Float64 loglambda_end_spc = log(m_lambdaRange_spc.GetEnd());
        if (m_lambdaRange_spc.GetBegin() > m_lambdaRange_tpl.GetBegin()) 
        {
            loglambda_start_spc = loglambda_start_tpl + ceil((loglambda_start_spc - loglambda_start_tpl)/m_logGridStep)*m_logGridStep; // ceil to be bigger or equal the first sample
            m_lambdaRange_spc.SetBegin(exp(loglambda_start_spc));
        }
        if (m_lambdaRange_spc.GetEnd() < m_lambdaRange_tpl.GetEnd())
        {
            loglambda_end_spc = loglambda_end_tpl - ceil((loglambda_end_tpl - loglambda_end_spc)/m_logGridStep)*m_logGridStep; // ceil to be less or equal the last sample
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
        Int32 nb_z = Int32( ceil((log_zmax_new_p1 - log_zmin_new_p1)/m_logGridStep) );
        zmax_new = exp(log_zmin_new_p1 + nb_z*m_logGridStep) - 1.;
    }
    m_zrange = TFloat64Range(zmin_new, zmax_new); //updating zrange based on the new logstep 

    InferTemplateRebinningSetup();

    return;
}

/**
*  Rebin the spectrum with the calculated logGridStep if spectrum not already rebinned:
*  step1: construct the spectralAxis
*  step2: do the rebin
*/
std::shared_ptr< CSpectrum> CSpectrumLogRebinning::LoglambdaRebinSpectrum( std::shared_ptr<const CSpectrum> spectrum, 
                                                                                std::string errorRebinMethod)
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
    bool verbose = true;  
    if(verbose){
        // save rebinned data
        FILE *f = fopen("loglbda_rebinlog_spclogrebin_dbg.txt", "w+");
        const CSpectrumSpectralAxis & w = spectrumRebinedLog->GetSpectralAxis();
        const CSpectrumFluxAxis & F = spectrumRebinedLog->GetFluxAxis();
        for (Int32 t = 0; t < spectrumRebinedLog->GetSampleCount(); t++)
        {
            fprintf(f, "%f\t%f\n", w[t], F[t]*1e16);
        }
        fclose(f);
    }
    return spectrumRebinedLog;  
}

/*
    Aims at computing the template lambda range once for all template catalog
*/
void CSpectrumLogRebinning::InferTemplateRebinningSetup()
{
    Float64 loglbdamin = log( m_lambdaRange_tpl.GetBegin()/(1.0 + m_zrange.GetEnd()));
    Float64 loglbdamax = log( m_lambdaRange_tpl.GetEnd()/(1.0 + m_zrange.GetBegin()));
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
std::shared_ptr<CTemplate> CSpectrumLogRebinning::LoglambdaRebinTemplate(std::shared_ptr<const CTemplate> tpl)
{
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


    m_loglambda_count_tpl = targetSpectralAxis.GetSamplesCount();

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
//to remove
if(tpl->GetName() == "ssp_5Gyr_z008.dat"){
        // save rebinned data
        FILE *f = fopen("loglbda_rebinlog_templatelogrebin_dbg6004.txt", "w+");
        const CSpectrumSpectralAxis & w = templateRebinedLog->GetSpectralAxis();
        const CSpectrumFluxAxis & F = templateRebinedLog->GetFluxAxis();
        for (Int32 t = 0; t < templateRebinedLog->GetSampleCount(); t++)
        {
            fprintf(f, "%f\t%e\n", w[t], F[t]);
        }
        fclose(f);

}

    return templateRebinedLog;
}

CSpectrumSpectralAxis CSpectrumLogRebinning::computeTargetLogSpectralAxis(TFloat64Range lambdarange, UInt32 count)
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