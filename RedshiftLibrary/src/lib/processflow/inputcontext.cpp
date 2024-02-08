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
#include <climits>

#include <boost/range/combine.hpp>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/tuple_like_boost_tuple.h"
#include "RedshiftLibrary/processflow/inputcontext.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/spectrum/template/template.h"

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
  m_use_LogLambaSpectrum =
      m_ParameterStore->hasToLogRebin(m_categories, fft_processing);
  if (!m_use_LogLambaSpectrum)
    return;
  for (auto const &[spectrum_ptr, lambdaRange_ptr] :
       boost::combine(m_spectra, m_lambdaRanges)) {
    if (spectrum_ptr->GetSpectralAxis().IsLogSampled()) {
      addRebinSpectrum(std::make_shared<CSpectrum>(spectrum_ptr->GetName()));
      CSpectrumSpectralAxis spcWav = spectrum_ptr->GetSpectralAxis();
      spcWav.RecomputePreciseLoglambda(); // in case input spectral values have
                                          // been rounded

      // Intersect with input lambdaRange and get indexes
      Int32 kstart = -1;
      Int32 kend = -1;
      bool ret = lambdaRange_ptr->getClosedIntervalIndices(
          spcWav.GetSamplesVector(), kstart, kend);
      if (!ret)
        THROWG(INTERNAL_ERROR,
               "LambdaRange borders are outside the spectralAxis range");
      // save into the rebinnedSpectrum
      m_rebinnedSpectra.back()->SetSpectralAndFluxAxes(
          spcWav.extract(kstart, kend),
          spectrum_ptr->GetFluxAxis().extract(kstart, kend));
      m_logGridStep =
          m_rebinnedSpectra.back()->GetSpectralAxis().GetlogGridStep();
    } else {
      m_logGridStep =
          m_ParameterStore->getMinZStepForFFTProcessing(fft_processing);
    }
  }
  Log.LogInfo(Formatter() << "loggrid step=" << m_logGridStep);
  std::string const errorRebinMethod = "rebinVariance";
  CSpectrumLogRebinning logReb(*this);

  for (auto const &[spectrum_ptr, lambdaRange_ptr,
                    rebinnedClampedLambdaRange_ptr] :
       boost::combine(m_spectra, m_lambdaRanges,
                      m_rebinnedClampedLambdaRanges)) {
    if (!spectrum_ptr->GetSpectralAxis().IsLogSampled())
      addRebinSpectrum(
          logReb.loglambdaRebinSpectrum(*spectrum_ptr, errorRebinMethod));

    TFloat64Range zrange;
    for (std::string cat : m_categories) {
      if (fft_processing[cat]) {
        zrange = logReb.logRebinTemplateCatalog(cat);
        m_logRebin.insert({cat, SRebinResults{zrange}});
      }
    }
    // Initialize rebinned clamped lambda range
    m_rebinnedSpectra.back()->GetSpectralAxis().ClampLambdaRange(
        *lambdaRange_ptr, *rebinnedClampedLambdaRange_ptr);
  }
  return;
}

void CInputContext::OrthogonalizeTemplates() {
  Float64 lambda =
      (m_lambdaRanges[0]->GetBegin() + m_lambdaRanges[0]->GetEnd()) / 2;
  if (m_spectra[0]->GetLSF() == nullptr)
    THROWG(INTERNAL_ERROR, "No defined lsf");
  Float64 resolution = CLSFGaussianConstantResolution::computeResolution(
      lambda, m_spectra[0]->GetLSF()->GetWidth(lambda));
  std::shared_ptr<TLSFArguments> args =
      std::make_shared<TLSFGaussianConstantResolutionArgs>(resolution);
  std::shared_ptr<const CLSF> lsf =
      LSFFactory.Create("gaussianConstantResolution", args);

  for (std::string cat : m_categories) {
    if (m_ParameterStore->HasToOrthogonalizeTemplates(cat)) {
      CTemplatesOrthogonalization tplOrtho;
      tplOrtho.Orthogonalize(*this, cat, lsf);
    }
  }
}

void CInputContext::setLineCatalog(
    const std::string &spectrumModel, const std::string &method,
    const std::shared_ptr<CLineCatalog> &catalog) {
  m_lineCatalogs[spectrumModel][method] = catalog;
}

void CInputContext::setLineRatioCatalogCatalog(
    const std::string &spectrumModel,
    const std::shared_ptr<CLineCatalogsTplRatio> &catalog) {
  m_lineRatioCatalogCatalogs[spectrumModel] = catalog;
}

void CInputContext::Init() {
  m_categories = m_ParameterStore->GetList<std::string>("spectrumModels");
  for (std::string cat : m_categories)
    Log.LogInfo(cat);

  // set template continuum removal parameters
  m_TemplateCatalog->InitContinuumRemoval(m_ParameterStore);

  bool enableInputSpcCorrect = m_ParameterStore->Get<bool>("autoCorrectInput");

  // non clamped lambdaRange: to be clamped depending on used spectra
  for (auto spectrum : m_spectra) {
    TFloat64Range lambdaRange;
    if (m_spectra.size() > 1)
      lambdaRange = m_ParameterStore->Get<TFloat64Range>(
          Formatter() << "lambdaRange." << spectrum->getObsID());
    else
      lambdaRange = m_ParameterStore->Get<TFloat64Range>("lambdaRange");
    m_lambdaRanges.push_back(std::make_shared<TFloat64Range>(lambdaRange));
    m_clampedLambdaRanges.emplace_back(new TFloat64Range());
    m_rebinnedClampedLambdaRanges.emplace_back(new TFloat64Range());
    m_constClampedLambdaRanges.push_back(m_clampedLambdaRanges.back());
    m_constRebinnedClampedLambdaRanges.push_back(
        m_rebinnedClampedLambdaRanges.back());
  }

  // clamp spectral axis at lambdaRange
  for (auto const &[spectrum_ptr, lambdaRange_ptr, clampedLambdaRange_ptr] :
       boost::combine(m_spectra, m_lambdaRanges, m_clampedLambdaRanges)) {
    spectrum_ptr->GetSpectralAxis().ClampLambdaRange(*lambdaRange_ptr,
                                                     *clampedLambdaRange_ptr);
  }

  // validate spectra and initialize spectra continuum removal
  for (auto const &[spectrum_ptr, lambdaRange_ptr] :
       boost::combine(m_spectra, m_lambdaRanges)) {
    spectrum_ptr->ValidateSpectrum(*lambdaRange_ptr, enableInputSpcCorrect);
    spectrum_ptr->InitSpectrumContinuum(*m_ParameterStore);
    // convolve IGM by LSF

    if (!m_igmcorrectionMeiksin->isConvolved() ||
        m_ParameterStore->Get<std::string>("lsf.lsfType") == "fromSpectrumData")
      m_igmcorrectionMeiksin->convolveByLSF(spectrum_ptr->GetLSF(),
                                            *lambdaRange_ptr);
  }

  // insert extinction correction objects if needed
  m_TemplateCatalog->m_logsampling = 0;
  m_TemplateCatalog->m_orthogonal = 0;
  m_TemplateCatalog->SetIsmIgmCorrection(
      m_ParameterStore, m_igmcorrectionMeiksin, m_ismcorrectionCalzetti);

  // log-lambda resampling if needed
  RebinInputs();

  // validate log-lambda resampled spectra
  if (m_use_LogLambaSpectrum) {
    for (auto const &[spectrum_ptr, lambdaRange_ptr, rebinedSpectrum_ptr] :
         boost::combine(m_spectra, m_lambdaRanges, m_rebinnedSpectra)) {
      rebinedSpectrum_ptr->ValidateSpectrum(*lambdaRange_ptr,
                                            enableInputSpcCorrect);
      rebinedSpectrum_ptr->SetLSF(spectrum_ptr->GetLSF());
    }
  }
  // template orthogonalisation with linemodel
  OrthogonalizeTemplates();
}

void CInputContext::resetSpectrumSpecific() {
  m_spectra.clear();
  m_constSpectra.clear();
  m_rebinnedSpectra.clear();
  m_constRebinnedSpectra.clear();
  m_lambdaRanges.clear();
  m_rebinnedClampedLambdaRanges.clear();
  m_clampedLambdaRanges.clear();
  m_constLambdaRanges.clear();
  m_constClampedLambdaRanges.clear();
  m_constRebinnedClampedLambdaRanges.clear();
  // not always spectrum specific
  m_TemplateCatalog.reset();
  // those one should not be here, they stay until api modification (only load
  // spectrum specific data in Context::run)
  m_lineCatalogs.clear();
  m_lineRatioCatalogCatalogs.clear();
  m_photBandCatalog.reset();
}
