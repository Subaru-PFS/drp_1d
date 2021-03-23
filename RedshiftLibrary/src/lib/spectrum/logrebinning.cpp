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
TFloat64Range& CSpectrumLogRebinning::Computelogstep(const TFloat64Range &lambdaRange, 
                                                    const Float64 zInputStep,
                                                    const TFloat64Range zInputRange)
{   

    if(spectrum.GetSpectralAxis().IsLogSampled() )
    { 
        m_loglambdaReference = spectrum.GetSpectralAxis()[0];
        m_logGridStep = spectrum.GetSpectralAxis()[1] - spectrum.GetSpectralAxis()[0];  
    }else{
        m_logGridStep = zInputStep;
        m_loglambdaReference = log(lambdaRange.GetBegin()); 
    }
    Log.LogDetail("  Log-Rebin: logGridStep = %f and loglambdaReference = %f", m_logGridStep, m_loglambdaReference);
    //computing the starting loglambda value considering that m_loglambdaReference should be included in the rebinned spectrum
    m_log_lambda_start = m_loglambdaReference + ceil((log(lambdaRange.GetBegin()) - m_loglambdaReference)/m_logGridStep)*m_logGridStep;

    //spectrum grid count using m_log_lambda_start
    m_spectrumLogGridCount = (Int32) floor((log(lambdaRange.GetEnd()) -  log_lambda_start) / m_logGridStep + 1);
    if(m_spectrumLogGridCount<2){
        Log.LogError("   Operator-TemplateFittingLog: logGridCount = %d <2", m_spectrumLogGridCount);
        throw runtime_error("   Operator-TemplateFittingLog: logGridCount <2. Abort");
    }
    Log.LogDetail("  Operator-TemplateFittingLog: Log-Rebin: logGridCount = %d", m_spectrumLogGridCount);

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
    InferTemplateRebinningSetup(lambdaRange);
    return m_zrange; //returned to be passed to processFlow
}

/**
*  Rebin the spectrum with the calculated logGridStep if spectrum not already rebinned:
*  step1: construct the spectralAxis
*  step2: do the rebin
*/
std::shared_ptr<const CSpectrum> CSpectrumLogRebinning::LoglambdaRebinSpectrum( CSpectrum& spectrum, 
                                                                                const TFloat64Range& lambdaRange, 
                                                                                std::string errorRebinMethod)
{    
    if(spectrum.GetSpectralAxis().IsLogSampled()){
        auto spectrumRebinedLog = make_shared<CSpectrum>(spectrum);
        m_rebinnedspcLambdaRange = TFloat64Range(spectrumRebinedLog->GetSpectralAxis().GetSamplesVector());
        return spectrumRebinedLog;
    }

    Float64 loglbdamin = log(lambdaRange.GetBegin()),
            loglbdamax = log(lambdaRange.GetEnd());

    TFloat64Range spcLbdaRange(exp(m_loglambdaReference - 0.5 * m_logGridStep),
                            exp(loglbdamax + 0.5 * m_logGridStep));
    Log.LogDetail("  Operator-TemplateFittingLog: Log-regular lambda resampling START");

    //prepare return rebinned vector
       
    auto spectrumRebinedLog = make_shared<CSpectrum>(spectrum.GetName());
    spectrumRebinedLog->GetSpectralAxis().SetSize(m_spectrumLogGridCount);
    spectrumRebinedLog->GetFluxAxis().SetSize(m_spectrumLogGridCount);   
    CMask mskRebinedLog(m_spectrumLogGridCount);

    CSpectrumSpectralAxis targetSpectralAxis =  computeTargetLogSpectralAxis(spectrum.GetSpectralAxis(),
                                                                            m_log_lambda_start,
                                                                            spectrumRebinedLog->GetSampleCount());
    // rebin the spectrum
    spectrum.Rebin(spcLbdaRange, targetSpectralAxis,
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
    m_rebinnedspcLambdaRange = TFloat64Range(spectrumRebinedLog->GetSpectralAxis().GetSamplesVector());
    return spectrumRebinedLog;  
}

/*
    Aims at computing the template lambda range once for all template catalog
*/
void CSpectrumLogRebinning::InferTemplateRebinningSetup(TFloat64Range lambdaRange)
{
    Float64 loglbdamin, loglbdamax;
    if(spectrum.GetSpectralAxis().IsLogSampled() )
        //below is valid if input spectrum is in log
        loglbdamin = log( m_rebinnedspcLambdaRange.GetBegin()/(1.0 + m_zrange.GetEnd()));
        loglbdamax = log( m_rebinnedspcLambdaRange.GetEnd()/ (1.0 + m_zrange.GetBegin()));
        m_templateLogGridCount = Int32(ceil((loglbdamax - loglbdamin)/m_logGridStep)) + 1;
    }else{
        //case of non-log spectrum
        //consider a full lambda range and not the rebinnedspcLambdaRange
        //we should get the rebinnedspcLambdaRange considering that the spectrum has a spectral axis corresponding or higher than lambdarange   
        loglbdamin = log(lambdaRange.GetBegin()/(1.0 + m_zrange.GetEnd()));//not sure about this
        m_templateLogGridCount = Int32(ceil((lambdaRange.GetEnd() - m_loglambdaReference)/m_logGridStep);
        loglbdamax = m_loglambdaReference + m_templateLogGridCount*m_logGridStep;
    }
    
    Float64 tgt_loglbdamax = loglbdamax;
    Float64 tgt_loglbdamin = loglbdamax - (loglbdamin-1)* m_logGridStep;
    m_tplLambdaRange = TFloat64Range(exp(tgt_loglbdamin - 0.5 * m_logGridStep), exp(tgt_loglbdamax + 0.5 * m_logGridStep));
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
    templateRebinedLog->GetSpectralAxis().SetSize(m_templateLogGridCount);
    templateRebinedLog->GetFluxAxis().SetSize(m_templateLogGridCount);
    templateRebinedLog->m_ismCorrectionCalzetti = tpl.m_ismCorrectionCalzetti;
    templateRebinedLog->m_igmCorrectionMeiksin = tpl.m_igmCorrectionMeiksin;

    CMask mskRebinedLog(m_templateLogGridCount);
	// Recheck the template coverage is larger,
	//  by construction only the min has to be re-checked,
	//  since the max is aligned to the max of the rebined spectrum at zmin which is
	// smaller to the max input lamdba range at min already checked 
	if (exp(m_tplLambdaRange.GetBegin()) < tpl.GetSpectralAxis()[0] )
	{
	    Log.LogError("  Operator-TemplateFittingLog: overlap found to be lower than 1.0 for this redshift range");
        Log.LogError("  Operator-TemplateFittingLog: for zmax=%f, tpl.lbdamin is %f (should be <%f)",
			        m_zrange.GetEnd(), tpl.GetSpectralAxis()[0], exp(m_tplLambdaRange.GetBegin()));
	    throw std::runtime_error("  Operator-TemplateFittingLog: overlap found to be lower than 1.0 for this redshift range");
	}

    CSpectrumSpectralAxis tpl_targetSpectralAxis =  computeTargetLogSpectralAxis(tpl.GetSpectralAxis(),
                                                                                m_tplLambdaRange.GetBegin(),
                                                                                m_templateLogGridCount);
    tpl.Rebin(  m_tplLambdaRange,
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
 * the target log spectral axis should contain the log(m_loglambdaReference)
 * Thus the construction of the new axis should be based on including this value
*/
CSpectrumSpectralAxis CSpectrumLogRebinning::computeTargetLogSpectralAxis(const CSpectrumSpectralAxis &ref_axis, 
                                                          Float64 tgt_loglbdamin,
                                                          Float64 logGridCount)
{
    TFloat64Range range(tgt_loglbdamin, tgt_loglbdamin + m_logGridStep*logGridCount);
    CSpectrumSpectralAxis targetSpectralAxis(range.SpreadOverLog(m_logGridStep), tgt_loglbdamin);
    return targetSpectralAxis;
}
