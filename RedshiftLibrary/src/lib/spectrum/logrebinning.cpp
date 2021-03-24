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

/**
 * Brief: Get loglambdastep and update the zrange accordingly
 * Below code relies on the fact that both loglambda grid and the log(Redshift+1) grid follows the same arithmetic progession with a common step
 * if spectrum is already rebinned, then it imposes the rebinning and the creation of zGrid 
 * Otherwise, it's the input zrange that decides on the rebinning param. 
 * Note : If we want the log step to be log(1.0+redshiftstep) then spreadOverlog should 
  be modified to use 1+redshiftstep for the common ratio (or construct the grid using arithmetic log progression).
*/
void CSpectrumLogRebinning::Computelogstep( CSpectrum &spectrum,
                                            const TFloat64Range &lambdaRange, 
                                            const Float64 zInputStep,
                                            const TFloat64Range zInputRange)
{   

    if(spectrum.GetSpectralAxis().IsLogSampled() )
    { 
        m_loglambdaRange = TLambdaRange(spectrum.GetSpectralAxis().GetSamplesVector());
        m_loglambdaRange_ref = m_loglambdaRange;
        m_logGridStep = spectrum.GetSpectralAxis()[1] - spectrum.GetSpectralAxis()[0];  
    }else{
        m_logGridStep = zInputStep;

        // compute reference lambda range (the effective lambda range of rebinned spectrum when initial spectrum overlaps lambdaRange)
        Float64 loglambda_start_ref = log(lambdaRange.GetBegin());
        Float64 loglambda_end_ref = log(lambdaRange.GetEnd());
        m_loglambda_count_ref =  Int32( floor(( loglambda_end_ref - loglambda_start_ref)/m_logGridStep) );// we should sample inside range hence floor
        loglambda_end_ref = loglambda_start_ref + m_loglambda_count_ref*m_logGridStep;
        m_loglambdaRange_ref = TFloat64Range(loglambda_start_ref, loglambda_end_ref);// this is the effective lambda range, to be used in InferTemplateRebinningSetup

        // compute rebinned spectrum lambda range 
        const TLambdaRange spclambdaRange = spectrum.GetSpectralAxis().GetLambdaRange();
        Float64 loglambda_start = loglambda_start_ref + ceil((log(spclambdaRange.GetBegin()) - loglambda_start_ref )/m_logGridStep)*m_logGridStep; // ceil to be bigger or equal the first sample
        m_loglambda_count =  Int32( floor(( log(spclambdaRange.GetEnd()) - loglambda_start)/m_logGridStep) ); // we should sample inside available samples hence floor
        Float64 loglambda_end = loglambda_start + m_loglambda_count*m_logGridStep;
        m_loglambdaRange = TLambdaRange(loglambda_start, loglambda_end); //to be passed to computeTargetLogSpectralAxis in LogLambdaRebinSpectrum
        if(m_loglambda_count<2){
            Log.LogError("   Operator-TemplateFittingLog: logGridCount = %d <2", m_loglambda_count);
            throw runtime_error("   Operator-TemplateFittingLog: logGridCount <2. Abort");
        }
    }
    Log.LogDetail("  Log-Rebin: logGridStep = %f, loglambdaReference = %f", m_logGridStep, m_loglambdaRange_ref);

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
    
    InferTemplateRebinningSetup(m_loglambdaRange_ref);

    return;
}

/**
*  Rebin the spectrum with the calculated logGridStep if spectrum not already rebinned:
*  step1: construct the spectralAxis
*  step2: do the rebin
*/
std::shared_ptr<const CSpectrum> CSpectrumLogRebinning::LoglambdaRebinSpectrum( std::shared_ptr<CSpectrum> spectrum, 
                                                                                const TFloat64Range& lambdaRange, 
                                                                                std::string errorRebinMethod)
{    
    if(spectrum->GetSpectralAxis().IsLogSampled()){
        return spectrum;
    }
    Log.LogDetail("  Operator-TemplateFittingLog: Log-regular lambda resampling START");

    //prepare return rebinned vector  
    auto spectrumRebinedLog = make_shared<CSpectrum>(spectrum->GetName());
    spectrumRebinedLog->GetSpectralAxis().SetSize(m_loglambda_count);
    spectrumRebinedLog->GetFluxAxis().SetSize(m_loglambda_count);   
    CMask mskRebinedLog(m_loglambda_count);

    CSpectrumSpectralAxis targetSpectralAxis =  computeTargetLogSpectralAxis(m_loglambdaRange, m_loglambda_count);

    TFloat64Range spcLbdaRange(exp(targetSpectralAxis[0] - 0.5 * m_logGridStep),
                            exp(targetSpectralAxis[m_loglambda_count-1] + 0.5 * m_logGridStep));
    // rebin the spectrum
    spectrum->Rebin(spcLbdaRange, targetSpectralAxis,
                    *spectrumRebinedLog,
                    mskRebinedLog, m_rebinMethod, errorRebinMethod);

    spectrumRebinedLog->GetSpectralAxis().SetLogScale(); //doesnt work before ::Rebin, cause we targetSpectralAxis into
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
void CSpectrumLogRebinning::InferTemplateRebinningSetup(TFloat64Range lambdaRange)
{
    Float64 loglbdamin = log( exp(lambdaRange.GetBegin())/(1.0 + m_zrange.GetEnd()));
    Float64 loglbdamax = log( exp(lambdaRange.GetEnd())/ (1.0 + m_zrange.GetBegin()));
    m_loglambda_count_ref = Int32(ceil((loglbdamax - loglbdamin)/m_logGridStep)) + 1;

    Float64 tgt_loglbdamax = loglbdamax;
    Float64 tgt_loglbdamin = loglbdamax - (m_loglambda_count_ref-1)* m_logGridStep;
    m_loglambdaRange_ref = TFloat64Range(exp(tgt_loglbdamin - 0.5 * m_logGridStep), exp(tgt_loglbdamax + 0.5 * m_logGridStep));
	Log.LogDetail("  Operator-TemplateFittingLog: Log-Rebin: tpl raw loglbdamin=%f : raw loglbdamax=%f", loglbdamin, loglbdamax);
	Log.LogDetail("  Operator-TemplateFittingLog: zmin_new = %f, tpl.lbdamax = %f", m_zrange.GetBegin(), exp(loglbdamax));
	Log.LogDetail("  Operator-TemplateFittingLog: zmax_new = %f, tpl.lbdamin = %f", m_zrange.GetEnd(), exp(loglbdamin));
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
std::shared_ptr<CTemplate> CSpectrumLogRebinning::LoglambdaRebinTemplate(const CTemplate &tpl)
{
    //InferTemplateRebinningSetup(); 
    auto templateRebinedLog = make_shared<CTemplate>(tpl.GetName(), tpl.GetCategory());
    templateRebinedLog->GetSpectralAxis().SetSize(m_loglambda_count_ref);
    templateRebinedLog->GetFluxAxis().SetSize(m_loglambda_count_ref);
    templateRebinedLog->m_ismCorrectionCalzetti = tpl.m_ismCorrectionCalzetti;
    templateRebinedLog->m_igmCorrectionMeiksin = tpl.m_igmCorrectionMeiksin;

    CMask mskRebinedLog(m_loglambda_count_ref);
	// Recheck the template coverage is larger,
	//  by construction only the min has to be re-checked,
	//  since the max is aligned to the max of the rebined spectrum at zmin which is
	// smaller to the max input lamdba range at min already checked 
	if (exp(m_loglambdaRange_ref.GetBegin()) < tpl.GetSpectralAxis()[0] )
	{
	    Log.LogError("  Operator-TemplateFittingLog: overlap found to be lower than 1.0 for this redshift range");
        Log.LogError("  Operator-TemplateFittingLog: for zmax=%f, tpl.lbdamin is %f (should be <%f)",
			        m_zrange.GetEnd(), tpl.GetSpectralAxis()[0], exp(m_loglambdaRange_ref.GetBegin()));
	    throw std::runtime_error("  Operator-TemplateFittingLog: overlap found to be lower than 1.0 for this redshift range");
	}

    CSpectrumSpectralAxis tpl_targetSpectralAxis =  computeTargetLogSpectralAxis(m_loglambdaRange_ref,
                                                                                m_loglambda_count_ref);
    tpl.Rebin(  m_loglambdaRange_ref,
                tpl_targetSpectralAxis,
                *templateRebinedLog,
                mskRebinedLog);
    
    templateRebinedLog->GetSpectralAxis().SetLogScale();
//to remove
if(tpl.GetName() == "ssp_5Gyr_z008.dat"){
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
/**
 * the target log spectral axis should contain the log(m_loglambdaRange_ref)
 * Thus the construction of the new axis should be based on including this value
*/
CSpectrumSpectralAxis CSpectrumLogRebinning::computeTargetLogSpectralAxis(TFloat64Range lambdarange, UInt32 count)
{
    CSpectrumSpectralAxis targetSpectralAxis(lambdarange.SpreadOverLog(m_logGridStep), count);
    return targetSpectralAxis;
}
