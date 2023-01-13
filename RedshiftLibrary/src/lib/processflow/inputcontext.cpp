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
#include <climits>
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

  auto lambdaRange = m_lambdaRanges.begin();
  for (auto spectrum_it = m_spectra.begin(); spectrum_it != m_spectra.end();
       ++spectrum_it, ++lambdaRange) {

    auto spectrum = *spectrum_it;
    if (spectrum->GetSpectralAxis().IsLogSampled()) {
      addRebinSpectrum(std::make_shared<CSpectrum>(spectrum->GetName()));
      CSpectrumSpectralAxis spcWav = spectrum->GetSpectralAxis();
      spcWav.RecomputePreciseLoglambda(); // in case input spectral values have
      // been rounded

      // Intersect with input lambdaRange and get indexes
      Int32 kstart = -1;
      Int32 kend = -1;
      bool ret = (*lambdaRange)
                     ->getClosedIntervalIndices(spcWav.GetSamplesVector(),
                                                kstart, kend);
      if (!ret)
        THROWG(INTERNAL_ERROR,
               "LambdaRange borders are outside the spectralAxis range");
      // save into the rebinnedSpectrum
      m_rebinnedSpectra.back()->SetSpectralAndFluxAxes(
          spcWav.extract(kstart, kend),
          spectrum->GetFluxAxis().extract(kstart, kend));
      m_logGridStep =
          m_rebinnedSpectra.back()->GetSpectralAxis().GetlogGridStep();
    } else {
      m_logGridStep =
          m_ParameterStore->getMinZStepForFFTProcessing(fft_processing);
    }
  }
  Log.LogInfo(Formatter() << "loggrid step=" << m_logGridStep);
  std::string category;
  std::string errorRebinMethod = "rebinVariance";
  CSpectrumLogRebinning logReb(*this);

  for (auto it = std::make_tuple(m_spectra.begin(), m_lambdaRanges.begin(),
                                 m_rebinnedClampedLambdaRanges.begin());
       std::get<0>(it) != m_spectra.end();
       ++std::get<0>(it), ++std::get<1>(it), ++std::get<2>(it)) {
    auto spectrum = *std::get<0>(it);
    auto lambdaRange = *std::get<1>(it);
    auto rebinnedClampedLambdaRange = *std::get<2>(it);
    if (!spectrum->GetSpectralAxis().IsLogSampled())
      addRebinSpectrum(
          logReb.loglambdaRebinSpectrum(spectrum, errorRebinMethod));

    TFloat64Range zrange;
    for (std::string cat : m_categories) {
      if (fft_processing[cat]) {
        zrange = logReb.logRebinTemplateCatalog(cat);
        m_logRebin.insert({cat, SRebinResults{zrange}});
      }
    }
    // Initialize rebinned clamped lambda range
    m_rebinnedSpectra.back()->GetSpectralAxis().ClampLambdaRange(
        *(lambdaRange), *(rebinnedClampedLambdaRange));
  }
  return;
}

void CInputContext::OrthogonalizeTemplates() {
  Float64 lambda =
      (m_lambdaRanges[0]->GetBegin() + m_lambdaRanges[0]->GetEnd()) / 2;
  Float64 resolution = CLSFGaussianConstantResolution::computeResolution(
      lambda, m_spectra[0]->GetLSF()->GetWidth(lambda));
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
    const std::string &objectType, const std::string &method,
    const std::shared_ptr<CLineCatalog> &catalog) {
  m_lineCatalogs[objectType][method] = catalog;
}

void CInputContext::setLineRatioCatalogCatalog(
    const std::string &objectType,
    const std::shared_ptr<CLineCatalogsTplRatio> &catalog) {
  m_lineRatioCatalogCatalogs[objectType] = catalog;
}

void CInputContext::Init() {
  m_categories = m_ParameterStore->GetList<std::string>("objects");
  for (std::string cat : m_categories)
    Log.LogInfo(cat);

  // set template continuum removal parameters
  m_TemplateCatalog->InitContinuumRemoval(m_ParameterStore);

  bool enableInputSpcCorrect = m_ParameterStore->Get<bool>("autocorrectinput");
  // non clamped lambdaRange: to be clamped depending on used spectra

  for (auto spectrum : m_spectra) {
    TFloat64Range lambdaRange;
    if (m_spectra.size() > 1)
      lambdaRange = m_ParameterStore->Get<TFloat64Range>(
          Formatter() << "lambdarange." << spectrum->getObsID());
    else
      lambdaRange = m_ParameterStore->Get<TFloat64Range>("lambdarange");
    m_lambdaRanges.push_back(std::make_shared<TFloat64Range>(lambdaRange));
    std::shared_ptr<TFloat64Range> clr =
        std::make_shared<TFloat64Range>(TFloat64Range());
    std::shared_ptr<TFloat64Range> rclr =
        std::make_shared<TFloat64Range>(TFloat64Range());
    m_clampedLambdaRanges.push_back(clr);
    m_rebinnedClampedLambdaRanges.push_back(rclr);
    m_constClampedLambdaRanges.push_back(clr);
    m_constRebinnedClampedLambdaRanges.push_back(rclr);
  }
  for (auto it = std::make_tuple(m_spectra.begin(), m_lambdaRanges.begin(),
                                 m_rebinnedClampedLambdaRanges.begin());
       std::get<0>(it) != m_spectra.end();
       ++std::get<0>(it), ++std::get<1>(it), ++std::get<2>(it)) {
    auto spectrum = *std::get<0>(it);
    auto lambdaRange = *std::get<1>(it);
    spectrum->ValidateSpectrum(*(lambdaRange), enableInputSpcCorrect);
    spectrum->InitSpectrum(*m_ParameterStore);
    // convolve IGM by LSF

    if (!m_igmcorrectionMeiksin->isConvolved() ||
        m_ParameterStore->Get<std::string>("LSF.LSFType") == "FROMSPECTRUMDATA")
      m_igmcorrectionMeiksin->convolveByLSF(spectrum->GetLSF(), *(lambdaRange));
  }
  // insert extinction correction objects if needed
  m_TemplateCatalog->m_logsampling = 0;
  m_TemplateCatalog->m_orthogonal = 0;
  m_TemplateCatalog->SetIsmIgmCorrection(
      m_ParameterStore, m_igmcorrectionMeiksin, m_ismcorrectionCalzetti);

  // log-lambda resampling if needed
  RebinInputs();

  if (m_use_LogLambaSpectrum) {
    for (auto it = std::make_tuple(m_spectra.begin(), m_lambdaRanges.begin(),
                                   m_rebinnedSpectra.begin());
         std::get<0>(it) != m_spectra.end();
         ++std::get<0>(it), ++std::get<1>(it), ++std::get<2>(it)) {
      auto spectrum = *std::get<0>(it);
      auto lambdaRange = *std::get<1>(it);

      auto rspectrum = *std::get<2>(it);
      rspectrum->ValidateSpectrum(*(lambdaRange), enableInputSpcCorrect);
      rspectrum->SetLSF(spectrum->GetLSF());
    }
  }
  OrthogonalizeTemplates();

  for (auto it = std::make_tuple(m_spectra.begin(), m_lambdaRanges.begin(),
                                 m_clampedLambdaRanges.begin());
       std::get<0>(it) != m_spectra.end();
       ++std::get<0>(it), ++std::get<1>(it), ++std::get<2>(it)) {
    auto spectrum = *std::get<0>(it);
    auto lambdaRange = *std::get<1>(it);
    auto clampedLambdaRange = *std::get<2>(it);

    spectrum->GetSpectralAxis().ClampLambdaRange(*(lambdaRange),
                                                 *(clampedLambdaRange));
  }
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
