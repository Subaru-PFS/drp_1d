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
#include "RedshiftLibrary/spectrum/logrebinning.h"
#include "RedshiftLibrary/log/log.h"
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

void CSpectrumLogRebinning::RebinInputs(CInputContext& inputContext, std::string category, Float64 redshiftStep)
{
    m_logGridStep = redshiftStep;
    std::string errorRebinMethod = "rebinVariance";//rebin error axis as well

    TFloat64Range   redshiftRange = inputContext.GetParameterStore()->Get<TFloat64Range>(category + ".redshiftrange");
    UInt32          SSratio = 1;
    if(inputContext.GetParameterStore()->Has<UInt32>( category + ".linemodelsolve.linemodel.firstpass.largegridstepratio"))
        SSratio = inputContext.GetParameterStore()->Get<UInt32>( category + ".linemodelsolve.linemodel.firstpass.largegridstepratio");

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
                   redshiftRange,
                   SSratio);
                   
    if(inputContext.GetSpectrum()->GetSpectralAxis().IsLogSampled()){
        inputContext.SetRebinnedSpectrum(spc);
    }else
        inputContext.SetRebinnedSpectrum(LoglambdaRebinSpectrum(inputContext.GetSpectrum(), errorRebinMethod));        
    
    //rebin templates using previously identified parameters,
    // rebin only if rebinning parameters are different from previously used ones
    std::shared_ptr<CTemplateCatalog> tplcat = inputContext.GetTemplateCatalog();
    const TStringList categoryList = tplcat->GetCategoryList();
    tplcat->m_logsampling = true;

        // check existence of already  & correctly logsampled templates
        const UInt32 ntpl = tplcat->GetTemplateCount(category);
        if (ntpl>0)
        {        
            for (UInt32 i=0; i<ntpl; i++)
            {
                std::shared_ptr<const CTemplate> tpl = tplcat->GetTemplate(category, i);
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
                    tplcat->m_logsampling = false;
                    std::shared_ptr<const CTemplate> input_tpl = tplcat->GetTemplateByName(TStringList{category}, tpl->GetName());
                    tplcat->m_logsampling = true;
                    tplcat->SetTemplate(LoglambdaRebinTemplate(input_tpl), i);// assigin the tpl pointer to a new rebined template
                }
            }            
        } else {
            // no rebined templates in the category: rebin all templates
            tplcat->m_logsampling = false;
            const TTemplateConstRefList TplList = std::const_pointer_cast<const CTemplateCatalog>(tplcat)->GetTemplateList(TStringList{category});
            tplcat->m_logsampling = true;
            for (auto tpl : TplList)
            {   
                std::shared_ptr<CTemplate> rebinnedTpl = LoglambdaRebinTemplate(tpl);
                tplcat->Add(rebinnedTpl);
            }
        }

    tplcat->m_logsampling = false;
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
                                            const TFloat64Range & zInputRange,
                                            UInt32 SubSamplingRatio)
{   
    TFloat64Range lambdaRange_ref;
    if(spectrum.GetSpectralAxis().IsLogSampled() )
    { 
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
    Log.LogInfo("  Log-regular lambda resampling START");

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
    Log.LogDetail("  Log-regular lambda resampling FINISHED");

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
