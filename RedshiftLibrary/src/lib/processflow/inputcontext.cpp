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
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include <float.h>
using namespace NSEpic;

CInputContext::CInputContext(std::shared_ptr<CParameterStore> paramStore)
    : m_ParameterStore(std::move(paramStore)) {}
/*
Two cases exist:
1. input spectrum is linear-sampled
2. input spectrum is log-sampled

C1 _Case1:
1. To log-sample it:
** determine log-step based on zgrid step (from parameter.json)
** lambda reference = lambdaRange.GetBegin() , thus independent from spectrum
spectral axis cause anyway we are gona rebin the spectrum
** use this log-step and lambda reference for template rebin as well

C2 _Case2:
1. Check that the spectrum is well rebinned, using
@logReb::CheckLoglambdaRebinSpectrum@
** determine logLambda step and log(Z+1) step and the reference lambda value
** use these params to rebin templates.
** in the case of running bunches, check that all spectra are rebinned the same
way, otherwise we have to rebin
** the template catalog differently for each spectrum.

Template rebinning:
* for the same reference lambda value and logLambda step, rebin ONCE the
template catalogs
* If any param changes, redo the template rebinning:
** Add getters for templates as follows: @GetTemplate(..., "log", logstep, ref
lambda)@
** the getter: 1)returns rebinnedTemplates OR 2) redo the template rebinning

Rebinning parameters for _Case2 should be extracted from m_Spectrum object, thus
the client has the responsibility to add these info to each spectrum
*/
void CInputContext::RebinInputs() {

  std::map<std::string, bool> fft_processing;
  m_use_LogLambaSpectrum = false;
  for (std::string cat : m_categories) {
    fft_processing[cat] = m_ParameterStore->HasFFTProcessing(cat);
    if (fft_processing[cat])
      m_use_LogLambaSpectrum = true;
  }
  if (fft_processing.find("star") != fft_processing.end() &&
      fft_processing["star"]) {
    throw GlobalException(INTERNAL_ERROR,
                          "FFT processing is not yet supported for stars");
  }

  if (!m_use_LogLambaSpectrum)
    return;

  if (m_Spectrum->GetSpectralAxis().IsLogSampled()) {
    m_rebinnedSpectrum = std::make_shared<CSpectrum>(m_Spectrum->GetName());
    CSpectrumSpectralAxis spcWav = m_Spectrum->GetSpectralAxis();
    spcWav.RecomputePreciseLoglambda(); // in case input spectral values have
                                        // been rounded

    // Intersect with input lambdaRange and get indexes
    Int32 kstart = -1;
    Int32 kend = -1;
    bool ret = m_lambdaRange.getClosedIntervalIndices(spcWav.GetSamplesVector(),
                                                      kstart, kend);
    if (!ret)
      throw GlobalException(
          INTERNAL_ERROR,
          "LambdaRange borders are outside the spectralAxis range");
    // save into the rebinnedSpectrum
    m_rebinnedSpectrum->SetSpectralAndFluxAxes(
        spcWav.extract(kstart, kend),
        m_Spectrum->GetFluxAxis().extract(kstart, kend));
    m_logGridStep = m_rebinnedSpectrum->GetSpectralAxis().GetlogGridStep();
  } else {
    m_logGridStep = DBL_MAX;
    for (std::string cat : m_categories) {
      if (fft_processing[cat]) {
        Float64 redshift_step =
            m_ParameterStore->Get<Float64>(cat + ".redshiftstep");

        if (redshift_step < m_logGridStep) {
          m_logGridStep = redshift_step;
        }
      }
    }
  }
  Log.LogInfo(Formatter() << "loggrid step=" << m_logGridStep);
  std::string category;
  std::string errorRebinMethod = "rebinVariance";
  CSpectrumLogRebinning logReb(*this);

  if (!m_Spectrum->GetSpectralAxis().IsLogSampled())
    m_rebinnedSpectrum =
        logReb.LoglambdaRebinSpectrum(m_Spectrum, errorRebinMethod);

  TFloat64Range zrange;
  for (std::string cat : m_categories) {
    if (fft_processing[cat]) {
      zrange = logReb.LogRebinTemplateCatalog(cat);
      m_logRebin.insert({cat, SRebinResults{zrange}});
    }
  }

  return;
}

void CInputContext::OrthogonalizeTemplates() {
  Float64 lambda = (m_lambdaRange.GetBegin() + m_lambdaRange.GetEnd()) / 2;
  Float64 resolution = CLSFGaussianConstantResolution::computeResolution(
      lambda, m_Spectrum->GetLSF()->GetWidth(lambda));
  std::shared_ptr<TLSFArguments> args =
      std::make_shared<TLSFGaussianConstantResolutionArgs>(resolution);
  std::shared_ptr<const CLSF> lsf =
      LSFFactory.Create("GaussianConstantResolution", args);

  for (std::string cat : m_categories) {
    if (m_ParameterStore->HasToOrthogonalizeTemplates(cat)) {
      CTemplatesOrthogonalization tplOrtho;
      tplOrtho.Orthogonalize(*this, cat, lsf);
    }
  }
}

void CInputContext::setLineCatalog(
    const std::string &objectType,
    const std::shared_ptr<CLineCatalog> &catalog) {
  m_lineCatalogs[objectType] = catalog;
}

void CInputContext::setLineRatioCatalogCatalog(
    const std::string &objectType,
    const std::shared_ptr<CLineCatalogsTplShape> &catalog) {
  m_lineRatioCatalogCatalogs[objectType] = catalog;
}

void CInputContext::Init() {
  m_categories = m_ParameterStore->GetList<std::string>("objects");
  for (std::string cat : m_categories)
    Log.LogInfo(cat);

  bool enableInputSpcCorrect = m_ParameterStore->Get<bool>("autocorrectinput");
  // non clamped lambdaRange: to be clamped depending on used spectra
  m_lambdaRange = m_ParameterStore->Get<TFloat64Range>("lambdarange");

  m_Spectrum->ValidateSpectrum(m_lambdaRange, enableInputSpcCorrect);
  m_Spectrum->InitSpectrum(*m_ParameterStore);

  // set template continuum removal parameters
  m_TemplateCatalog->InitContinuumRemoval(m_ParameterStore);

  // convolve IGM by LSF
  if (!m_igmcorrectionMeiksin->isConvolved() ||
      m_ParameterStore->Get<std::string>("LSF.LSFType") == "FROMSPECTRUMDATA")
    m_igmcorrectionMeiksin->convolveByLSF(m_Spectrum->GetLSF(), m_lambdaRange);

  // Calzetti ISM & Meiksin IGM initialization, for only original templates,
  // only when lsf changes notably when LSFType is fromspectrumdata
  // or the first time InitIsmIgm is called
  m_TemplateCatalog->m_logsampling = 0;
  m_TemplateCatalog->m_orthogonal = 0;
  if (m_TemplateCatalog->GetTemplate(m_TemplateCatalog->GetCategoryList()[0], 0)
          ->CalzettiInitFailed()) {
    m_TemplateCatalog->InitIsmIgm(m_igmcorrectionMeiksin,
                                  m_ismcorrectionCalzetti);
  }

  RebinInputs();

  if (m_use_LogLambaSpectrum) {
    m_rebinnedSpectrum->ValidateSpectrum(m_lambdaRange, enableInputSpcCorrect);
    m_rebinnedSpectrum->SetLSF(m_Spectrum->GetLSF());
  }

  OrthogonalizeTemplates();
}

void CInputContext::resetSpectrumSpecific() {
  m_Spectrum.reset();
  m_rebinnedSpectrum.reset();
  // not always spectrum specific
  m_TemplateCatalog.reset();
  // those one should not be here, they stay until api modification (only load
  // spectrum specific data in Context::run)
  m_lineCatalogs.clear();
  m_lineRatioCatalogCatalogs.clear();
  m_photBandCatalog.reset();
}
