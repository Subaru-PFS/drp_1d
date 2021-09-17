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
#include "RedshiftLibrary/processflow/inputcontext.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include <float.h>
using namespace NSEpic;

CInputContext::CInputContext(std::shared_ptr<CSpectrum> spc,
                             std::shared_ptr<CTemplateCatalog> tmplCatalog,
                             std::shared_ptr<CRayCatalog> gal_rayCatalog,
                             std::shared_ptr<CRayCatalog> qso_rayCatalog,
                             std::shared_ptr<CParameterStore> paramStore):

  m_Spectrum(std::move(spc)),
  m_TemplateCatalog(std::move(tmplCatalog)),
  m_gal_RayCatalog(std::move(gal_rayCatalog)),
  m_qso_RayCatalog(std::move(qso_rayCatalog)),
  m_ParameterStore(std::move(paramStore))
{
    m_Spectrum->InitSpectrum(*m_ParameterStore);
    //non clamped lambdaRange: to be clamped depending on used spectra
    m_lambdaRange = m_ParameterStore->Get<TFloat64Range>("lambdarange");
    
    RebinInputs();

    // Calzetti ISM & Meiksin IGM initialization, for both rebinned and original templates
    std::string calibrationPath =  m_ParameterStore->Get<std::string>( "calibrationDir");  
    m_TemplateCatalog->InitIsmIgm(calibrationPath, m_ParameterStore, m_Spectrum->GetLSF());

    std::string enableInputSpcCorrectStr = m_ParameterStore->Get<std::string>( "autocorrectinput");
    Bool enableInputSpcCorrect = enableInputSpcCorrectStr == "yes";
    //we should replace spectrum with m_inputContext->
    validateSpectrum(m_Spectrum, m_lambdaRange, enableInputSpcCorrect);
    if(m_use_LogLambaSpectrum)
    {
        validateSpectrum(m_rebinnedSpectrum, m_lambdaRange, enableInputSpcCorrect);
        m_rebinnedSpectrum->SetLSF(m_Spectrum->GetLSF());
    }
    //orthog only if rebinning happens
    OrthogonalizeTemplates(calibrationPath);
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
    std::vector<std::string> categories{"galaxy", "qso", "star"};
    Bool fft_processing_gal = m_ParameterStore->HasFFTProcessing(categories[0]); 
    Bool fft_processing_qso = m_ParameterStore->HasFFTProcessing(categories[1]);
    Bool fft_processing_star = m_ParameterStore->HasFFTProcessing(categories[2]);
    
    if(fft_processing_star)
    {
        Log.LogError("FFT processing is not yet supported for stars");
        throw std::runtime_error("FFT processing is not yet supported for stars");
    }

    m_use_LogLambaSpectrum = fft_processing_gal || fft_processing_qso || fft_processing_star;

    if(!m_use_LogLambaSpectrum) return;

    Float64 logGridStep;
    if(m_Spectrum->GetSpectralAxis().IsLogSampled())
    {
        m_rebinnedSpectrum = std::make_shared<CSpectrum>(m_Spectrum->GetName());
        CSpectrumSpectralAxis  spcWav = m_Spectrum->GetSpectralAxis();
        spcWav.RecomputePreciseLoglambda(); // in case input spectral values have been rounded
        //save into the rebinnedSpectrum
        m_rebinnedSpectrum->SetSpectralAndFluxAxes(std::move(spcWav), m_Spectrum->GetFluxAxis());
        logGridStep = m_rebinnedSpectrum->GetSpectralAxis().GetlogGridStep();
    }else
    {
      Float64 zInputStep_gal = fft_processing_gal?m_ParameterStore->Get<Float64>( categories[0]+".redshiftstep" ):DBL_MAX;
      Float64 zInputStep_qso = fft_processing_qso?m_ParameterStore->Get<Float64>( categories[1]+".redshiftstep" ):DBL_MAX;        
      logGridStep = (zInputStep_gal>zInputStep_qso)?zInputStep_qso:zInputStep_gal;
    }
    std::string category;
    std::string errorRebinMethod = "rebinVariance";
    CSpectrumLogRebinning logReb(*this, logGridStep);

    if(!m_Spectrum->GetSpectralAxis().IsLogSampled())
      m_rebinnedSpectrum = logReb.LoglambdaRebinSpectrum(m_Spectrum, errorRebinMethod);

    TFloat64Range zrange;
    if(fft_processing_gal){
      zrange = logReb.LogRebinTemplateCatalog(categories[0]);
      m_logRebin.insert({categories[0], std::move(SRebinResults{zrange, logGridStep})});
    }
    if(fft_processing_qso){
      zrange = logReb.LogRebinTemplateCatalog(categories[1]);
      m_logRebin.insert({categories[1], std::move(SRebinResults{zrange, logGridStep})});
    }

    return;
}


void CInputContext::validateSpectrum(std::shared_ptr<CSpectrum> spectrum, 
                                                TFloat64Range lambdaRange, 
                                                Bool enableInputSpcCorrect)
{
  TFloat64Range clampedlambdaRange;
  spectrum->GetSpectralAxis().ClampLambdaRange(lambdaRange, clampedlambdaRange);
  Log.LogInfo( "Validate spectrum: (CLambdaRange: %f-%f:%f)",
               clampedlambdaRange.GetBegin(),
               clampedlambdaRange.GetEnd(),
               spectrum->GetResolution());

  Float64 lmin = clampedlambdaRange.GetBegin();
  Float64 lmax = clampedlambdaRange.GetEnd();

  if(enableInputSpcCorrect)
  {
      //Check if the Spectrum is valid on the lambdarange
      //correctInputSpectrum(ctx.GetInputContext()->m_lambdaRange);

      if( spectrum->correctSpectrum( lmin,lmax ))
        Log.LogInfo( "Successfully corrected noise on wavelength range (%.1f ; %.1f)",  lmin, lmax );
  }

   if( !spectrum->IsFluxValid( lmin, lmax ) ){
      Log.LogError("Failed to validate spectrum flux on wavelength range (%.1f ; %.1f)",
                   lmin, lmax );
      throw std::runtime_error("Failed to validate spectrum flux");
    }else{
      Log.LogDetail( "Successfully validated spectrum flux, on wavelength range (%.1f ; %.1f)", lmin, lmax );
    }
	//Check if the noise is valid in the clampedlambdaRange
    if( !spectrum->IsNoiseValid( lmin, lmax ) ){
      Log.LogError("Failed to validate noise on wavelength range (%.1f ; %.1f)",
                   lmin, lmax );
      throw std::runtime_error("Failed to validate noise from spectrum");
    }else{
      Log.LogDetail( "Successfully validated noise on wavelength range (%.1f ; %.1f)", lmin, lmax );
    }
}

void CInputContext::OrthogonalizeTemplates(const std::string& calibrationPath)
{
    Bool orthog_gal = m_ParameterStore->HasToOrthogonalizeTemplates("galaxy"); 
    Bool orthog_qso = m_ParameterStore->HasToOrthogonalizeTemplates("qso");

    std::shared_ptr<const CLSF> lsf = m_Spectrum->GetLSF(); //to be changed in #6680
    
    if(orthog_gal)
    {
      CTemplatesOrthogonalization tplOrtho;
      tplOrtho.Orthogonalize(*this,"galaxy",calibrationPath, lsf);
    }
    if(orthog_qso)
    {
      CTemplatesOrthogonalization tplOrtho_;//tplOrtho could be reused..TBC
      tplOrtho_.Orthogonalize(*this,"qso",calibrationPath, lsf);
    }
    return;
}