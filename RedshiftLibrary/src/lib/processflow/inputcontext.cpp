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
    Bool enableInputSpcCorrect = m_ParameterStore->Get<bool>( "autocorrectinput");
    //non clamped lambdaRange: to be clamped depending on used spectra
    m_lambdaRange = m_ParameterStore->Get<TFloat64Range>("lambdarange");

    m_Spectrum->ValidateSpectrum(m_lambdaRange, enableInputSpcCorrect);
    m_Spectrum->InitSpectrum(*m_ParameterStore); 
    
    // set template continuum removal parameters
    m_TemplateCatalog->InitContinuumRemoval(m_ParameterStore);
  
    // Calzetti ISM & Meiksin IGM initialization, for only original templates, 
    //only when lsf changes notably when LSFType is fromspectrumdata
    //or the first time InitIsmIgm is called
    m_TemplateCatalog->m_logsampling = 0; m_TemplateCatalog->m_orthogonal = 0; 
    if(m_TemplateCatalog->GetTemplate(m_TemplateCatalog->GetCategoryList()[0], 0)->CalzettiInitFailed())
    {
        m_TemplateCatalog->InitIsmIgm(m_ParameterStore, m_Spectrum->GetLSF());
    }
    else
      {
      if(m_ParameterStore->Get<std::string>("LSF.LSFType") == "FROMSPECTRUMDATA") //redo the convolution
      {
        m_TemplateCatalog->GetTemplate(m_TemplateCatalog->GetCategoryList()[0], 0)->m_igmCorrectionMeiksin->ConvolveAll(m_Spectrum->GetLSF());
      }
    }

    RebinInputs();

    if(m_use_LogLambaSpectrum)
    {
        m_rebinnedSpectrum->ValidateSpectrum(m_lambdaRange, enableInputSpcCorrect);
        m_rebinnedSpectrum->SetLSF(m_Spectrum->GetLSF());
    }

    OrthogonalizeTemplates();
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
    Bool fft_processing_gal = m_ParameterStore->HasFFTProcessing(m_categories[0]); 
    Bool fft_processing_qso = m_ParameterStore->HasFFTProcessing(m_categories[1]);
    Bool fft_processing_star = m_ParameterStore->HasFFTProcessing(m_categories[2]);
    
    if(fft_processing_star)
    {
        throw GlobalException(INTERNAL_ERROR,"FFT processing is not yet supported for stars");
    }

    m_use_LogLambaSpectrum = fft_processing_gal || fft_processing_qso || fft_processing_star;

    if(!m_use_LogLambaSpectrum) return;

    if(m_Spectrum->GetSpectralAxis().IsLogSampled())
    {
        m_rebinnedSpectrum = std::make_shared<CSpectrum>(m_Spectrum->GetName());
        CSpectrumSpectralAxis  spcWav = m_Spectrum->GetSpectralAxis();
        spcWav.RecomputePreciseLoglambda(); // in case input spectral values have been rounded

        //Intersect with input lambdaRange and get indexes
        Int32 kstart = -1; Int32 kend = -1;
        bool ret = m_lambdaRange.getClosedIntervalIndices(spcWav.GetSamplesVector(), kstart, kend);
        if(!ret) throw GlobalException(INTERNAL_ERROR, "LambdaRange borders are outside the spectralAxis range");
        //save into the rebinnedSpectrum
        m_rebinnedSpectrum->SetSpectralAndFluxAxes(spcWav.extract(kstart,kend), m_Spectrum->GetFluxAxis().extract(kstart,kend));
        m_logGridStep = m_rebinnedSpectrum->GetSpectralAxis().GetlogGridStep();
    }else
    {
      Float64 zInputStep_gal = fft_processing_gal?m_ParameterStore->Get<Float64>( m_categories[0]+".redshiftstep" ):DBL_MAX;
      Float64 zInputStep_qso = fft_processing_qso?m_ParameterStore->Get<Float64>( m_categories[1]+".redshiftstep" ):DBL_MAX;        
      m_logGridStep = (zInputStep_gal>zInputStep_qso)?zInputStep_qso:zInputStep_gal;
    }
    std::string category;
    std::string errorRebinMethod = "rebinVariance";
    CSpectrumLogRebinning logReb(*this);

    if(!m_Spectrum->GetSpectralAxis().IsLogSampled())
      m_rebinnedSpectrum = logReb.LoglambdaRebinSpectrum(m_Spectrum, errorRebinMethod);

    TFloat64Range zrange;
    if(fft_processing_gal){
      zrange = logReb.LogRebinTemplateCatalog(m_categories[0]);
      m_logRebin.insert({m_categories[0], SRebinResults{zrange}});
    }
    if(fft_processing_qso){
      zrange = logReb.LogRebinTemplateCatalog(m_categories[1]);
      m_logRebin.insert({m_categories[1], SRebinResults{zrange}});
    }

    return;
}

void CInputContext::OrthogonalizeTemplates()
{
    Bool orthog_gal = m_ParameterStore->HasToOrthogonalizeTemplates( m_categories[0]); 
    Bool orthog_qso = m_ParameterStore->HasToOrthogonalizeTemplates( m_categories[1]);

    Float64 lambda = (m_lambdaRange.GetBegin() + m_lambdaRange.GetEnd())/2;
    std::shared_ptr<TLSFArguments> args = std::make_shared<TLSFGaussianConstantWidthArgs>(m_Spectrum->GetLSF()->GetWidth(lambda));
    std::shared_ptr<const CLSF> lsf = LSFFactory.Create("GaussianConstantWidth", args);

    if(orthog_gal)
    {
      CTemplatesOrthogonalization tplOrtho;
      tplOrtho.Orthogonalize(*this, m_categories[0],lsf);
    }
    if(orthog_qso)
    {
      CTemplatesOrthogonalization tplOrtho_;//tplOrtho could be reused..TBC
      tplOrtho_.Orthogonalize(*this, m_categories[1],lsf);
    }
    return;
}
