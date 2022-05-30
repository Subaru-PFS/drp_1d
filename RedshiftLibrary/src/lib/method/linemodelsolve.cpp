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
#include "RedshiftLibrary/method/linemodelsolve.h"

#include "RedshiftLibrary/log/log.h"

#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"

#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/statistics/deltaz.h"
#include "RedshiftLibrary/statistics/pdfcandidateszresult.h"
#include "RedshiftLibrary/statistics/zprior.h"

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <fstream>
#include <iostream>
#include <string>

using namespace NSEpic;
using namespace std;
using namespace boost;
// throw GlobalException\(\s*\n*ErrorCode::([A-Z_]+),((?:\s*".+"\n*)+)\);
// throw
// GlobalException\(\s*\n*ErrorCode::([A-Z_]+),\n*\s+Formatter\(\)\s*\n*((?:\s*<<\s*".+"\n*)+)\);
// throw
// GlobalException\(\s*\n*ErrorCode::([A-Z_]+),\n*\s+Formatter\(\)\s*\n*<<((?:\s*\s*".+"\n*)+)\);
// throw
// GlobalException\(\s*\n*ErrorCode::([A-Z_]+),\n*\s+Formatter\(\)\s*\n*((?:\s*<<\s*"?.+"?\n*)+)\);
/**
 * \brief Empty constructor.
 **/
CLineModelSolve::CLineModelSolve(TScopeStack &scope, string objectType)
    : CObjectSolve("LineModelSolve", scope, objectType) {}

/**
 * \brief
 * Populates the method parameters from the dataStore into the class members
 * Returns true if successful, false otherwise
 **/
bool CLineModelSolve::PopulateParameters(
    std::shared_ptr<const CParameterStore> parameterStore) {
  m_redshiftSeparation = parameterStore->Get<Float64>(
      "extremaredshiftseparation"); //,m_redshiftSeparation,2e-3);
  m_opt_linetypefilter =
      parameterStore->GetScoped<std::string>("linemodel.linetypefilter");
  m_opt_lineforcefilter =
      parameterStore->GetScoped<std::string>("linemodel.lineforcefilter");
  m_opt_fittingmethod =
      parameterStore->GetScoped<std::string>("linemodel.fittingmethod");

  m_opt_secondpasslcfittingmethod = parameterStore->GetScoped<std::string>(
      "linemodel.secondpasslcfittingmethod");
  m_opt_skipsecondpass =
      parameterStore->GetScoped<bool>("linemodel.skipsecondpass");
  m_opt_secondpass_continuumfit = parameterStore->GetScoped<std::string>(
      "linemodel.secondpass.continuumfit");
  m_opt_secondpass_halfwindowsize =
      parameterStore->GetScoped<Float64>("linemodel.secondpass.halfwindowsize");

  m_opt_firstpass_fittingmethod = parameterStore->GetScoped<std::string>(
      "linemodel.firstpass.fittingmethod");
  m_opt_firstpass_largegridstepRatio = parameterStore->GetScoped<Int32>(
      "linemodel.firstpass.largegridstepratio");
  m_opt_firstpass_tplratio_ismfit =
      parameterStore->GetScoped<bool>("linemodel.firstpass.tplratio_ismfit");
  m_opt_firstpass_disablemultiplecontinuumfit = parameterStore->GetScoped<bool>(
      "linemodel.firstpass.multiplecontinuumfit_disable");

  m_opt_firstpass_largegridsampling = m_redshiftSampling;
  Log.LogDetail("    firstpass - largegridsampling (auto set from "
                "redshiftsampling param.): %s",
                m_opt_firstpass_largegridsampling.c_str());

  m_opt_continuumcomponent =
      parameterStore->GetScoped<std::string>("linemodel.continuumcomponent");
  if (m_opt_continuumcomponent == "tplfit" ||
      m_opt_continuumcomponent == "tplfitauto") {
    m_opt_tplfit_fftprocessing =
        parameterStore->GetScoped<bool>("linemodel.continuumfit.fftprocessing");
    m_opt_tplfit_use_photometry = false;
    if (parameterStore->HasScoped<bool>("linemodel.enablephotometry")) {
      m_opt_tplfit_use_photometry =
          parameterStore->GetScoped<bool>("linemodel.enablephotometry");
      if (m_opt_tplfit_use_photometry)
        m_opt_tplfit_photo_weight =
            parameterStore->GetScoped<Float64>("linemodel.photometry.weight");
    }
    if (m_opt_tplfit_fftprocessing && m_opt_tplfit_use_photometry)
      THROWG(INTERNAL_ERROR,
             "CLineModelSolve::PopulateParameters: fftprocessing "
             "not implemented with photometry enabled");
    m_opt_tplfit_fftprocessing_secondpass = m_opt_tplfit_fftprocessing;
    m_opt_tplfit_dustfit =
        parameterStore->GetScoped<bool>("linemodel.continuumfit.ismfit");
    m_opt_tplfit_igmfit =
        parameterStore->GetScoped<bool>("linemodel.continuumfit.igmfit");
    m_opt_continuumfitcount =
        parameterStore->GetScoped<Int32>("linemodel.continuumfit.count");
    m_opt_continuum_neg_amp_threshold = parameterStore->GetScoped<Float64>(
        "linemodel.continuumfit.negativethreshold");
    m_opt_continuum_null_amp_threshold = parameterStore->GetScoped<Float64>(
        "linemodel.continuumfit.nullthreshold");
    m_opt_tplfit_ignoreLinesSupport = parameterStore->GetScoped<bool>(
        "linemodel.continuumfit.ignorelinesupport");
    m_opt_tplfit_continuumprior_betaA = parameterStore->GetScoped<Float64>(
        "linemodel.continuumfit.priors.betaA");
    m_opt_tplfit_continuumprior_betaTE = parameterStore->GetScoped<Float64>(
        "linemodel.continuumfit.priors.betaTE");
    m_opt_tplfit_continuumprior_betaZ = parameterStore->GetScoped<Float64>(
        "linemodel.continuumfit.priors.betaZ");
    m_opt_tplfit_continuumprior_dirpath =
        parameterStore->GetScoped<std::string>(
            "linemodel.continuumfit.priors.catalog_dirpath"); // no priors by
                                                              // default
  }

  /*if(m_opt_igmfit!=m_opt_tplfit_igmfit)
      throw GlobalException(INCOHERENT_INPUTPARAMETERS, "igmfit should be
     activated for both continuum and linemodel");*/

  m_opt_rigidity = parameterStore->GetScoped<std::string>("linemodel.rigidity");

  if (m_opt_rigidity == "tplshape") {
    m_opt_tplratio_reldirpath =
        parameterStore->GetScoped<std::string>("linemodel.tplratio_catalog");
    m_opt_tplratio_ismfit =
        parameterStore->GetScoped<bool>("linemodel.tplratio_ismfit");

    m_opt_tplratio_prior_betaA =
        parameterStore->GetScoped<Float64>("linemodel.tplratio.priors.betaA");
    m_opt_tplratio_prior_betaTE =
        parameterStore->GetScoped<Float64>("linemodel.tplratio.priors.betaTE");
    m_opt_tplratio_prior_betaZ =
        parameterStore->GetScoped<Float64>("linemodel.tplratio.priors.betaZ");
    m_opt_tplratio_prior_dirpath = parameterStore->GetScoped<std::string>(
        "linemodel.tplratio.priors.catalog_dirpath"); // no priors by default
  } else if (m_opt_rigidity == "rules") {
    m_opt_enableImproveBalmerFit =
        parameterStore->GetScoped<bool>("linemodel.improveBalmerFit");
  }

  m_opt_lineWidthType =
      parameterStore->GetScoped<std::string>("linemodel.linewidthtype");
  m_opt_nsigmasupport =
      parameterStore->GetScoped<Float64>("linemodel.nsigmasupport");
  m_opt_velocity_emission =
      parameterStore->GetScoped<Float64>("linemodel.velocityemission");
  m_opt_velocity_absorption =
      parameterStore->GetScoped<Float64>("linemodel.velocityabsorption");
  m_opt_velocityfit = parameterStore->GetScoped<bool>("linemodel.velocityfit");
  if (m_opt_velocityfit) {
    m_opt_em_velocity_fit_min =
        parameterStore->GetScoped<Float64>("linemodel.emvelocityfitmin");
    m_opt_em_velocity_fit_max =
        parameterStore->GetScoped<Float64>("linemodel.emvelocityfitmax");
    m_opt_em_velocity_fit_step =
        parameterStore->GetScoped<Float64>("linemodel.emvelocityfitstep");
    m_opt_abs_velocity_fit_min =
        parameterStore->GetScoped<Float64>("linemodel.absvelocityfitmin");
    m_opt_abs_velocity_fit_max =
        parameterStore->GetScoped<Float64>("linemodel.absvelocityfitmax");
    m_opt_abs_velocity_fit_step =
        parameterStore->GetScoped<Float64>("linemodel.absvelocityfitstep");
  }
  m_opt_lya_forcefit = parameterStore->GetScoped<bool>("linemodel.lyaforcefit");
  m_opt_lya_forcedisablefit =
      parameterStore->GetScoped<bool>("linemodel.lyaforcedisablefit");
  m_opt_lya_fit_asym_min =
      parameterStore->GetScoped<Float64>("linemodel.lyafit.asymfitmin");
  m_opt_lya_fit_asym_max =
      parameterStore->GetScoped<Float64>("linemodel.lyafit.asymfitmax");
  m_opt_lya_fit_asym_step =
      parameterStore->GetScoped<Float64>("linemodel.lyafit.asymfitstep");
  m_opt_lya_fit_width_min =
      parameterStore->GetScoped<Float64>("linemodel.lyafit.widthfitmin");
  m_opt_lya_fit_width_max =
      parameterStore->GetScoped<Float64>("linemodel.lyafit.widthfitmax");
  m_opt_lya_fit_width_step =
      parameterStore->GetScoped<Float64>("linemodel.lyafit.widthfitstep");
  m_opt_lya_fit_delta_min =
      parameterStore->GetScoped<Float64>("linemodel.lyafit.deltafitmin");
  m_opt_lya_fit_delta_max =
      parameterStore->GetScoped<Float64>("linemodel.lyafit.deltafitmax");
  m_opt_lya_fit_delta_step =
      parameterStore->GetScoped<Float64>("linemodel.lyafit.deltafitstep");

  m_opt_continuumreest =
      parameterStore->GetScoped<std::string>("linemodel.continuumreestimation");
  m_opt_rules = parameterStore->GetScoped<std::string>("linemodel.rules");
  m_opt_extremacount =
      parameterStore->GetScoped<Int32>("linemodel.extremacount");
  m_opt_extremacountB =
      parameterStore->GetScoped<Int32>("linemodel.extremacountB");
  m_opt_candidatesLogprobaCutThreshold =
      parameterStore->GetScoped<Float64>("linemodel.extremacutprobathreshold");
  m_opt_stronglinesprior =
      parameterStore->GetScoped<Float64>("linemodel.stronglinesprior");
  m_opt_haPrior = parameterStore->GetScoped<Float64>("linemodel.haprior");
  m_opt_euclidNHaEmittersPriorStrength =
      parameterStore->GetScoped<Float64>("linemodel.euclidnhaemittersStrength");
  m_opt_pdfcombination =
      parameterStore->GetScoped<std::string>("linemodel.pdfcombination");
  m_opt_pdf_margAmpCorrection =
      parameterStore->GetScoped<bool>("linemodel.pdf.margampcorr");

  // Auto-correct fitting method
  std::string forcefittingmethod = "individual";
  if (m_opt_rigidity == "tplshape" && m_opt_fittingmethod == "hybrid") {
    THROWG(BAD_PARAMETER_VALUE, "rigidity = tplshape and fitting_method=hybrid "
                                "imply fittingmethod=individual");
    /*
    m_opt_fittingmethod = forcefittingmethod;
      parameterStore->SetScopedParam("linemodel.fittingmethod",
    m_opt_fittingmethod); Log.LogInfo( "Linemodel fitting method auto-correct
    due to tplshape rigidity");
    */
  }
  if (m_opt_rigidity == "tplshape" &&
      m_opt_firstpass_fittingmethod == "hybrid") {
    THROWG(BAD_PARAMETER_VALUE,
           "rigidity = tplshape and firstpass_fitting_method=hybrid imply "
           "fittingmethod=individual");

    //   m_opt_firstpass_fittingmethod = forcefittingmethod;
    //  parameterStore.SetScopedParam("linemodel.firstpass.fittingmethod",
    //  m_opt_firstpass_fittingmethod);
    // Log.LogInfo( "Linemodel first pass fitting method auto-correct due to
    // tplshape rigidity");
  }

  Log.LogInfo("Linemodel parameters:");
  Log.LogInfo("    -linetypefilter: %s", m_opt_linetypefilter.c_str());
  Log.LogInfo("    -lineforcefilter: %s", m_opt_lineforcefilter.c_str());
  Log.LogInfo("    -fittingmethod: %s", m_opt_fittingmethod.c_str());
  Log.LogInfo("    -linewidthtype: %s", m_opt_lineWidthType.c_str());
  if (m_opt_lineWidthType == "combined" ||
      m_opt_lineWidthType == "velocitydriven") {
    Log.LogInfo("    -velocity emission: %.2f", m_opt_velocity_emission);
    Log.LogInfo("    -velocity absorption: %.2f", m_opt_velocity_absorption);
    Log.LogInfo("    -velocity fit: %s", m_opt_velocityfit ? "true" : "false");
  }

  Log.LogInfo("    -nsigmasupport: %.1f", m_opt_nsigmasupport);

  if (m_opt_velocityfit) {
    Log.LogInfo("    -em velocity fit min : %.1f", m_opt_em_velocity_fit_min);
    Log.LogInfo("    -em velocity fit max : %.1f", m_opt_em_velocity_fit_max);
    Log.LogInfo("    -em velocity fit step : %.1f", m_opt_em_velocity_fit_step);
    Log.LogInfo("    -abs velocity fit min : %.1f", m_opt_abs_velocity_fit_min);
    Log.LogInfo("    -abs velocity fit max : %.1f", m_opt_abs_velocity_fit_max);
    Log.LogInfo("    -abs velocity fit step : %.1f",
                m_opt_abs_velocity_fit_step);
  }

  Log.LogInfo("    -rigidity: %s", m_opt_rigidity.c_str());
  if (m_opt_rigidity == "rules") {
    Log.LogInfo("      -rules: %s", m_opt_rules.c_str());
    if (m_opt_fittingmethod == "hybrid") {
      Log.LogInfo("      -Improve Balmer Fit: %s",
                  m_opt_enableImproveBalmerFit ? "true" : "false");
    }
  } else if (m_opt_rigidity == "tplshape") {
    Log.LogInfo("      -tplratio_catalog: %s",
                m_opt_tplratio_reldirpath.c_str());
    Log.LogInfo("      -tplratio_ismfit: %s",
                m_opt_tplratio_ismfit ? "true" : "false");
    Log.LogInfo("      -tplfit_priors_dirpath: %s",
                m_opt_tplratio_prior_dirpath.c_str());
    Log.LogInfo("      -tplratio_priors_betaA:  %f",
                m_opt_tplratio_prior_betaA);
    Log.LogInfo("      -tplratio_priors_betaTE:  %f",
                m_opt_tplratio_prior_betaTE);
    Log.LogInfo("      -tplratio_priors_betaZ:  %f",
                m_opt_tplratio_prior_betaZ);
  }

  Log.LogInfo("    -continuumcomponent: %s", m_opt_continuumcomponent.c_str());
  if (m_opt_continuumcomponent == "tplfit" ||
      m_opt_continuumcomponent == "tplfitauto") {
    Log.LogInfo("      -tplfit_fftprocessing: %d", m_opt_tplfit_fftprocessing);
    Log.LogInfo("      -tplfit_ismfit: %s",
                m_opt_tplfit_dustfit ? "true" : "false");
    Log.LogInfo("      -tplfit_igmfit: %s",
                m_opt_tplfit_igmfit ? "true" : "false");
    Log.LogInfo("      -continuum fit count:  %.0f", m_opt_continuumfitcount);
    Log.LogInfo("      -tplfit_ignorelinesupport: %s",
                m_opt_tplfit_ignoreLinesSupport ? "true" : "false");
    Log.LogInfo("      -tplfit_secondpass-LC-fitting-method: %s",
                m_opt_secondpasslcfittingmethod.c_str());
    Log.LogInfo("      -tplfit_priors_dirpath: %s",
                m_opt_tplfit_continuumprior_dirpath.c_str());
    Log.LogInfo("      -tplfit_priors_betaA:  %f",
                m_opt_tplfit_continuumprior_betaA);
    Log.LogInfo("      -tplfit_priors_betaTE:  %f",
                m_opt_tplfit_continuumprior_betaTE);
    Log.LogInfo("      -tplfit_priors_betaZ:  %f",
                m_opt_tplfit_continuumprior_betaZ);
  }
  Log.LogInfo("    -continuumreestimation: %s", m_opt_continuumreest.c_str());
  Log.LogInfo("    -extremacount: %i", m_opt_extremacount);
  Log.LogInfo("    -extremacount-firstpass B: %i", m_opt_extremacountB);
  Log.LogInfo("    -extrema cut proba-threshold: %.0f",
              m_opt_candidatesLogprobaCutThreshold);
  Log.LogInfo("    -first pass:");
  Log.LogInfo("      -largegridstepratio: %d",
              m_opt_firstpass_largegridstepRatio);
  Log.LogInfo("      -fittingmethod: %s",
              m_opt_firstpass_fittingmethod.c_str());
  Log.LogInfo("      -tplratio_ismfit: %s",
              m_opt_firstpass_tplratio_ismfit ? "true" : "false");
  Log.LogInfo("      -multiplecontinuumfit_disable: %s",
              m_opt_firstpass_disablemultiplecontinuumfit ? "true" : "false");

  Log.LogInfo("    -second pass:");
  Log.LogInfo("      -skip second pass: %s",
              m_opt_skipsecondpass ? "true" : "false");
  Log.LogInfo("      -continuum fit method: %s",
              m_opt_secondpass_continuumfit.c_str());

  Log.LogInfo("    -pdf-stronglinesprior: %e", m_opt_stronglinesprior);
  Log.LogInfo("    -pdf-hapriorstrength: %e", m_opt_haPrior);
  Log.LogInfo("    -pdf-euclidNHaEmittersPriorStrength: %e",
              m_opt_euclidNHaEmittersPriorStrength);
  Log.LogInfo("    -pdf-combination: %s",
              m_opt_pdfcombination
                  .c_str()); // "marg";    // "bestchi2";    // "bestproba";
  Log.LogInfo("    -pdf-margAmpCorrection: %s",
              m_opt_pdf_margAmpCorrection ? "true" : "false");

  return true;
}

/**
 * \brief Calls the Solve method and returns a new "result" object.
 * Call Solve.
 * Return a pointer to an empty CLineModelSolveResult. (The results for
 *Linemodel will reside in the linemodel.linemodel result).
 **/
std::shared_ptr<CSolveResult>
CLineModelSolve::compute(std::shared_ptr<const CInputContext> inputContext,
                         std::shared_ptr<COperatorResultStore> resultStore,
                         TScopeStack &scope) {
  const CSpectrum &rebinnedSpc = *(inputContext->GetRebinnedSpectrum());
  const CSpectrum &spc = *(inputContext->GetSpectrum());
  const CTemplateCatalog &tplCatalog = *(inputContext->GetTemplateCatalog());
  const CLineCatalog &restlinecatalog =
      *(inputContext->GetLineCatalog(m_objectType, m_name));
  const auto &photBandCat = inputContext->GetPhotBandCatalog();
  const CLineCatalogsTplShape &tplRatioCatalog =
      *(inputContext->GetTemplateRatioCatalog(m_objectType));

  PopulateParameters(inputContext->GetParameterStore());

  // useloglambdasampling param is relevant only if linemodel.continuumfit is
  // set to use fftprocessing below we explicit this check on this condition
  bool useloglambdasampling =
      inputContext->GetParameterStore()->GetScoped<bool>(
          "linemodel.useloglambdasampling");
  useloglambdasampling &= inputContext->GetParameterStore()->GetScoped<bool>(
      "linemodel.continuumfit.fftprocessing");

  CLineCatalog::TLineVector restLineList = restlinecatalog.GetFilteredList(
      m_opt_linetypefilter, m_opt_lineforcefilter);
  Log.LogDebug("restLineList.size() = %d", restLineList.size());

  bool retSolve =
      Solve(resultStore, useloglambdasampling ? rebinnedSpc : spc, rebinnedSpc,
            tplCatalog, std::move(restLineList), tplRatioCatalog, m_lambdaRange,
            m_redshifts, photBandCat, m_opt_tplfit_use_photometry);

  if (!retSolve) {
    return NULL;
  }

  auto results = resultStore->GetScopedGlobalResult("linemodel");
  if (results.expired()) {
    Log.LogError("linemodelsolve: Unable to retrieve linemodel results");
    return NULL;
  }
  std::shared_ptr<const CLineModelResult> result =
      std::dynamic_pointer_cast<const CLineModelResult>(results.lock());

  //  suggestion : CSolve::GetCurrentScopeName(TScopeStack)
  //  prepare the linemodel chisquares and prior results for pdf computation
  ChisquareArray chisquares = BuildChisquareArray(result);

  /*
  Log.LogDetail("    linemodelsolve: Storing priors (size=%d)",
  zpriorResult->Redshifts.size());

  resultStore->StoreScopedGlobalResult( "priorpdf", zpriorResult); //TODO review
  name if uncommented
  */

  /* ------------------------  COMPUTE POSTERIOR PDF  --------------------------
   */

  COperatorPdfz pdfz(
      m_opt_pdfcombination,
      0.0,                // no peak Separation in 2nd pass
      0.0,                // cut threshold
      m_opt_extremacount, // max nb of final (2nd pass) candidates
      "SPE",              // Id_prefix
      false,              // do not allow extrema at border
      1,                  // one peak/window only
      m_linemodel.m_secondpass_parameters_extremaResult.ExtendedRedshifts,
      m_linemodel.m_secondpass_parameters_extremaResult.m_ranked_candidates);

  std::shared_ptr<PdfCandidatesZResult> candidateResult =
      pdfz.Compute(chisquares);

  // store PDF results
  Log.LogInfo("%s: Storing PDF results", __func__);

  resultStore->StoreScopedGlobalResult("pdf", pdfz.m_postmargZResult);

  TFloat64Range clampedlambdaRange;
  spc.GetSpectralAxis().ClampLambdaRange(m_lambdaRange, clampedlambdaRange);
  // Get linemodel results at extrema (recompute spectrum model etc.)
  std::shared_ptr<LineModelExtremaResult> ExtremaResult =
      m_linemodel.buildExtremaResults(
          spc, clampedlambdaRange, candidateResult->m_ranked_candidates,
          m_opt_continuumreest); // maybe its better to pass
                                 // resultStore->GetGlobalResult so that we
                                 // constuct extremaResult all at once in
                                 // linemodel operator
  // store extrema results
  storeExtremaResults(resultStore, ExtremaResult);

  // SaveContinuumPDF(dataStore, result);
  //  TBD

  // create the solveresult
  std::shared_ptr<CLineModelSolveResult> lmsolveresult =
      std::make_shared<CLineModelSolveResult>(
          ExtremaResult->m_ranked_candidates[0].second, m_opt_pdfcombination,
          pdfz.m_postmargZResult->valMargEvidenceLog);

  return lmsolveresult;
}

void CLineModelSolve::GetZpriorsOptions(
    bool &zPriorStrongLinePresence, bool &zPriorHaStrongestLine,
    bool &zPriorNLineSNR, Float64 &opt_nlines_snr_penalization_factor,
    bool &zPriorEuclidNHa) const {
  zPriorStrongLinePresence = (m_opt_stronglinesprior > 0.0);
  if (zPriorStrongLinePresence) {
    Log.LogDetail("%s: StrongLinePresence prior enabled: factor=%e", __func__,
                  m_opt_stronglinesprior);
  } else {
    Log.LogDetail("%s: StrongLinePresence prior disabled", __func__);
  }

  zPriorHaStrongestLine = (m_opt_haPrior > 0.0);
  if (zPriorHaStrongestLine) {
    Log.LogDetail("%s: Ha strongest line prior enabled: factor=%e", __func__,
                  m_opt_haPrior);
  } else {
    Log.LogDetail("%s: Ha strongest line prior disabled", __func__);
  }

  opt_nlines_snr_penalization_factor = -1;
  zPriorNLineSNR = (opt_nlines_snr_penalization_factor > 0.0);
  if (zPriorNLineSNR) {
    Log.LogDetail("%s: N lines snr>cut prior enabled: factor=%e", __func__,
                  opt_nlines_snr_penalization_factor);
  } else {
    Log.LogDetail("%s: N lines snr>cut prior disabled", __func__);
  }

  // hardcoded Euclid-NHaZprior parameter
  zPriorEuclidNHa = false;
  if (m_opt_euclidNHaEmittersPriorStrength > 0.0) {
    zPriorEuclidNHa = true;
    Log.LogDetail("%s: EuclidNHa prior enabled, with strength-coeff: %e",
                  __func__, m_opt_euclidNHaEmittersPriorStrength);
  } else {
    Log.LogDetail("%s: EuclidNHa prior disabled", __func__);
  }
}

TFloat64List CLineModelSolve::BuildZpriors(
    const std::shared_ptr<const CLineModelResult> &result,
    Int32 kTplShape) const {
  TFloat64List zpriors;

  CZPrior zpriorhelper;

  bool zPriorStrongLinePresence, zPriorHaStrongestLine, zPriorEuclidNHa,
      zPriorNLineSNR;
  Float64 opt_nlines_snr_penalization_factor;
  GetZpriorsOptions(zPriorStrongLinePresence, zPriorHaStrongestLine,
                    zPriorNLineSNR, opt_nlines_snr_penalization_factor,
                    zPriorEuclidNHa);

  if (zPriorStrongLinePresence) {
    if (kTplShape == -1) {
      Int32 lineTypeFilter = 1; // for emission lines only
      const TBoolList strongLinePresence = result->getStrongLinesPresence(
          lineTypeFilter, result->LineModelSolutions);
      zpriors = zpriorhelper.GetStrongLinePresenceLogZPrior(
          strongLinePresence, m_opt_stronglinesprior);
    } else {
      const TBoolList &strongLinePresence =
          result->StrongELPresentTplshapes[kTplShape];
      zpriors = zpriorhelper.GetStrongLinePresenceLogZPrior(
          strongLinePresence, m_opt_stronglinesprior);
    }

  } else {
    zpriors = zpriorhelper.GetConstantLogZPrior(result->Redshifts.size());
  }

  if (zPriorHaStrongestLine) {
    TFloat64List zlogPriorHaStrongest;
    if (kTplShape == -1) {
      const TBoolList wHaStronglinePresence =
          result->getStrongestLineIsHa(result->LineModelSolutions);
      zlogPriorHaStrongest = zpriorhelper.GetStrongLinePresenceLogZPrior(
          wHaStronglinePresence, m_opt_haPrior);
    } else {
      const TBoolList &wHaStronglinePresence =
          result->StrongHalphaELPresentTplshapes[kTplShape];
      zlogPriorHaStrongest = zpriorhelper.GetStrongLinePresenceLogZPrior(
          wHaStronglinePresence, m_opt_haPrior);
    }
    zpriors = zpriorhelper.CombineLogZPrior(zpriors, zlogPriorHaStrongest);
  }

  if (zPriorEuclidNHa) {
    TFloat64List zlogPriorNHa = zpriorhelper.GetEuclidNhaLogZPrior(
        result->Redshifts, m_opt_euclidNHaEmittersPriorStrength);
    zpriors = zpriorhelper.CombineLogZPrior(zpriors, zlogPriorNHa);
  }

  if (zPriorNLineSNR) {
    TFloat64List zlogPriorNLinesAboveSNR;
    TInt32List n_lines_above_snr;
    if (kTplShape == -1) {
      const TInt32List n_lines_above_snr =
          result->getNLinesAboveSnrcut(result->LineModelSolutions);
      zlogPriorNLinesAboveSNR = zpriorhelper.GetNLinesSNRAboveCutLogZPrior(
          n_lines_above_snr, opt_nlines_snr_penalization_factor);
    } else {
      const TInt32List &n_lines_above_snr =
          result->NLinesAboveSNRTplshapes[kTplShape];
      zlogPriorNLinesAboveSNR = zpriorhelper.GetNLinesSNRAboveCutLogZPrior(
          n_lines_above_snr, opt_nlines_snr_penalization_factor);
    }
    zpriors = zpriorhelper.CombineLogZPrior(zpriors, zlogPriorNLinesAboveSNR);
  }

  return zpriors;
}

ChisquareArray CLineModelSolve::BuildChisquareArray(
    const std::shared_ptr<const CLineModelResult> &result) const {
  Log.LogDetail("LinemodelSolve: building chisquare array");

  ChisquareArray chisquarearray;
  std::vector<TFloat64List> &chisquares = chisquarearray.chisquares;
  std::vector<TFloat64List> &zpriors = chisquarearray.zpriors;

  chisquarearray.cstLog = result->cstLog;
  Log.LogDetail("%s: using cstLog = %f", __func__, chisquarearray.cstLog);

  chisquarearray.redshifts = result->Redshifts;

  const Int32 zsize = result->Redshifts.size();

  if (m_opt_pdfcombination == "bestchi2") {
    zpriors.push_back(BuildZpriors(result));
    chisquares.push_back(result->ChiSquare);

  } else if (m_opt_pdfcombination == "bestproba" ||
             m_opt_pdfcombination == "marg") {

    if (m_opt_rigidity != "tplshape") {
      zpriors.push_back(BuildZpriors(result));
      chisquares.push_back(result->ChiSquare);
    } else {

      const Int32 ntplshapes = result->ChiSquareTplshapes.size();

      bool zPriorLines = false;
      Log.LogDetail("%s: PriorLinesTplshapes.size()=%d", __func__,
                    result->PriorLinesTplshapes.size());
      if (boost::filesystem::exists(m_opt_tplratio_prior_dirpath) &&
          result->PriorLinesTplshapes.size() == ntplshapes) {
        zPriorLines = true;
        Log.LogDetail("%s: Lines Prior enabled", __func__);
      } else {
        Log.LogDetail("%s: Lines Prior disabled", __func__);
      }

      chisquares.reserve(ntplshapes);
      zpriors.reserve(ntplshapes);

      for (Int32 k = 0; k < ntplshapes; ++k) {
        zpriors.push_back(BuildZpriors(result, k));
        chisquares.push_back(result->ChiSquareTplshapes[k]);

        // correct chi2 if necessary
        TFloat64List &logLikelihoodCorrected = chisquares.back();
        /*
        if(m_opt_pdf_margAmpCorrection) //nb: this is experimental.
        {
            //find max scalemargcorr
            Float64 maxscalemargcorr=-DBL_MAX;
            for ( Int32 kz=0; kz<zsize; kz++ )
                if(maxscalemargcorr <
        result->ScaleMargCorrectionTplshapes[k][kz]) maxscalemargcorr =
        result->ScaleMargCorrectionTplshapes[k][kz];

            Log.LogDetail("%s: maxscalemargcorr= %e", __func__,
        maxscalemargcorr); for ( Int32 kz=0; kz<zsize; kz++ )
                if(result->ScaleMargCorrectionTplshapes[k][kz]!=0) //warning,
        this is experimental. logLikelihoodCorrected[kz] +=
        result->ScaleMargCorrectionTplshapes[k][kz] - maxscalemargcorr;

            // need to add maxscalemargcorr ?
        }*/

        if (zPriorLines && result->PriorLinesTplshapes[k].size() == zsize)
          for (Int32 kz = 0; kz < zsize; kz++)
            logLikelihoodCorrected[kz] += result->PriorLinesTplshapes[k][kz];
      }

      chisquarearray.modelpriors = result->PriorTplshapes;
    }

    if (!result->ChiSquareTplContinuum.empty()) {

      // Fullmodel (ie with continuum template fitting): store all continuum tpl
      // fitting chisquares (ChiSquareTplContinuum size will be null if not
      // tplfit)

      //  note: the continuum xi2 are computed on continuum ~orthogonal to the
      //  linemodel, ie corresponding
      //        to a fullmodel with perfect linemodel fit (null Xi2 under the
      //        lines). They are thus better (smaller) than the fullmodel Xi2
      //        with which they will be summed in the marginalization. To
      //        mitigate this, for each continuum Xi2, we have to add the Xi2
      //        part of the linemodel alone. The latter can be obtained by
      //        subtracting the corresponding continuum Xi2 of the fullmodel.

      const auto newsize =
          chisquares.size() * result->ChiSquareTplContinuum.size();
      chisquares.reserve(newsize);
      zpriors.reserve(newsize);
      chisquarearray.modelpriors.reserve(newsize);

      // divide model prior by the number of continuum templates
      // TODO: need to add/handle tpl continuum priors
      // (note: in the case of ATEZ, which is exclusive of zpriors and
      // Tplshapes, priors are included already in the chisquares of both
      // tplcontinuum chi2 and tplshape chi2)
      for (auto &prior : chisquarearray.modelpriors)
        prior /= result->ChiSquareTplContinuum.size();

      // loop on all tplshapes (or 1 linemodel free)
      // TFloat64List linemodel_alone_chi2(zsize);
      for (Int32 k = 0, ke = chisquares.size(); k < ke; ++k) {

        // loop on all continuum templates, skiping 1st one, already used with
        // linemodel
        for (auto it = result->ChiSquareTplContinuum.cbegin() + 1,
                  itend = result->ChiSquareTplContinuum.cend();
             it != itend; ++it) {

          zpriors.push_back(zpriors[k]); // duplicate zpriors
          chisquares.emplace_back(zsize);
          TFloat64List &fullmodel_chi2 = chisquares.back();
          for (Int32 iz = 0; iz < zsize;
               ++iz) { // estimate chi2 of fullmodel with other continuum
            fullmodel_chi2[iz] =
                chisquares[k][iz] +
                std::max(0., (*it)[iz] - result->ChiSquareTplContinuum[0][iz]);
          }
          if (!chisquarearray.modelpriors.empty())
            chisquarearray.modelpriors.push_back(
                chisquarearray.modelpriors[k]); // duplicate modelpriors
        }
      }
    }

  } else {
    Log.LogError("Linemodel: Unable to parse pdf combination method option");
  }

  return chisquarearray;
}

///
/// \brief COperatorLineModel::storeGlobalModelResults
/// stores the linemodel results as global results in the datastore
///

void CLineModelSolve::storeExtremaResults(
    std::shared_ptr<COperatorResultStore> resultStore,
    std::shared_ptr<const LineModelExtremaResult> ExtremaResult) const {
  resultStore->StoreScopedGlobalResult("extrema_results", ExtremaResult);

  Int32 nResults = ExtremaResult->size();
}

void CLineModelSolve::StoreChisquareTplShapeResults(
    std::shared_ptr<COperatorResultStore> resultStore,
    std::shared_ptr<const CLineModelResult> result) const {
  for (Int32 km = 0; km < result->ChiSquareTplshapes.size(); km++) {
    std::shared_ptr<CLineModelResult> result_chisquaretplshape =
        std::shared_ptr<CLineModelResult>(new CLineModelResult());
    result_chisquaretplshape->Init(result->Redshifts, result->restLineList, 0,
                                   0, TFloat64List());
    for (Int32 kz = 0; kz < result->Redshifts.size(); kz++) {
      result_chisquaretplshape->ChiSquare[kz] =
          result->ChiSquareTplshapes[km][kz];
    }

    std::string resname =
        (boost::format(
             "linemodel_chisquaretplshape/linemodel_chisquaretplshape_%d") %
         km)
            .str();
    resultStore->StoreScopedGlobalResult(resname.c_str(),
                                         result_chisquaretplshape);
  }

  // Save scaleMargCorrTplshape results
  for (Int32 km = 0; km < result->ScaleMargCorrectionTplshapes.size(); km++) {
    std::shared_ptr<CLineModelResult> result_chisquaretplshape =
        std::shared_ptr<CLineModelResult>(new CLineModelResult());
    result_chisquaretplshape->Init(result->Redshifts, result->restLineList, 0,
                                   0, TFloat64List());
    for (Int32 kz = 0; kz < result->Redshifts.size(); kz++) {
      result_chisquaretplshape->ChiSquare[kz] =
          result->ScaleMargCorrectionTplshapes[km][kz];
    }

    std::string resname =
        (boost::format(
             "linemodel_chisquaretplshape/linemodel_scalemargcorrtplshape_%d") %
         km)
            .str();
    resultStore->StoreScopedGlobalResult(resname.c_str(),
                                         result_chisquaretplshape);
  }

  // Save PriorLinesTplshapes results
  for (Int32 km = 0; km < result->PriorLinesTplshapes.size(); km++) {
    std::shared_ptr<CLineModelResult> result_chisquaretplshape =
        std::shared_ptr<CLineModelResult>(new CLineModelResult());
    result_chisquaretplshape->Init(result->Redshifts, result->restLineList, 0,
                                   0, TFloat64List());
    for (Int32 kz = 0; kz < result->Redshifts.size(); kz++) {
      result_chisquaretplshape->ChiSquare[kz] =
          result->PriorLinesTplshapes[km][kz];
    }

    std::string resname =
        (boost::format(
             "linemodel_chisquaretplshape/linemodel_priorlinestplshape_%d") %
         km)
            .str();
    resultStore->StoreScopedGlobalResult(resname.c_str(),
                                         result_chisquaretplshape);
  }

  // Save PriorContinuumTplshapes results
  std::shared_ptr<CLineModelResult> result_chisquaretplshape =
      std::shared_ptr<CLineModelResult>(new CLineModelResult());
  result_chisquaretplshape->Init(result->Redshifts, result->restLineList, 0, 0,
                                 TFloat64List());
  for (Int32 kz = 0; kz < result->Redshifts.size(); kz++) {
    result_chisquaretplshape->ChiSquare[kz] =
        result->ContinuumModelSolutions[kz].tplLogPrior;
  }

  std::string resname =
      (boost::format(
           "linemodel_chisquaretplshape/linemodel_priorcontinuumtplshape"))
          .str();
  resultStore->StoreScopedGlobalResult(resname.c_str(),
                                       result_chisquaretplshape);
}

/**
 * \brief
 * Retrieve the true-velocities from a hardcoded ref file path
 * nb: this is a hack for development purposes
 **/
Int32 getVelocitiesFromRefFile(const char *filePath, std::string spcid,
                               Float64 &elv, Float64 &alv) {
  std::ifstream file;

  file.open(filePath, std::ifstream::in);
  if (file.rdstate() & ios_base::failbit)
    return false;

  string line;

  // Read file line by line
  while (getline(file, line)) {
    // remove comments
    if (line.compare(0, 1, "#", 1) == 0) {
      continue;
    }
    char_separator<char> sep(" \t");

    // Tokenize each line
    typedef tokenizer<char_separator<char>> ttokenizer;
    ttokenizer tok(line, sep);

    // Check if it's not a comment
    ttokenizer::iterator it = tok.begin();
    if (it != tok.end() && *it != "#") {
      string name;
      if (it != tok.end()) {
        name = *it;
      }
      std::size_t foundstr = name.find(spcid.c_str());
      if (foundstr == std::string::npos) {
        continue;
      }

      // Found the correct spectrum ID: now read the ref values
      Int32 nskip = 7;
      for (Int32 i = 0; i < nskip; i++) {
        ++it;
      }
      if (it != tok.end()) {

        elv = 0.0;
        try {
          elv = lexical_cast<double>(*it);
        } catch (bad_lexical_cast &) {
          elv = 0.0;
          return false;
        }
      }
      ++it;
      if (it != tok.end()) {
        alv = 0.0;
        try {
          alv = lexical_cast<double>(*it);
        } catch (bad_lexical_cast &) {
          alv = 0.0;
          return false;
        }
      }
    }
  }
  file.close();
  return true;
}

/**
 * \brief
 * Create a continuum object by subtracting spcWithoutContinuum from the spc.
 * Configure the opt_XXX variables from the dataStore scope parameters.
 * LogInfo the opt_XXX values.
 * Create a COperatorLineModel, call its Compute method.
 * If that returned true, store results.
 **/

bool CLineModelSolve::Solve(
    std::shared_ptr<COperatorResultStore> resultStore, const CSpectrum &spc,
    const CSpectrum &rebinnedSpc, const CTemplateCatalog &tplCatalog,
    const CLineCatalog::TLineVector &restLineList,
    const CLineCatalogsTplShape &tplRatioCatalog,
    const TFloat64Range &lambdaRange, const TFloat64List &redshifts,
    const std::shared_ptr<const CPhotBandCatalog> &photBandCat,
    const Float64 photo_weight) {
  std::string scopeStr = "linemodel";
  //    //Hack: load the simulated true-velocities
  //    if(false)
  //    {
  //        Float64 elv = 0.0;
  //        Float64 alv = 0.0;
  //        namespace fs = boost::filesystem;
  //        fs::path
  //        refFilePath("/home/aschmitt/data/simu_linemodel/simulm_20160513/simulation_pfswlinemodel_20160513_10spcperbin/refz.txt");
  //        if ( fs::exists(refFilePath) )
  //        {
  //            std::string spcSubStringId = spc.GetName().substr(0, 20);
  //            getVelocitiesFromRefFile( refFilePath.c_str(), spcSubStringId,
  //            elv, alv);
  //        }
  //        Float64 offsetv = 0.0;
  //        m_opt_velocity_emission = elv+offsetv;
  //        m_opt_velocity_absorption = alv+offsetv;
  //        Log.LogInfo( "Linemodel - hack - Loaded velocities for spc %s :
  //        elv=%4.1f, alv=%4.1f", spc.GetName().c_str(), elv, alv);
  //    }

  // Compute with linemodel operator
  Int32 retInit =
      m_linemodel.Init(spc, redshifts, restLineList, m_categoryList,
                       m_opt_continuumcomponent, m_opt_nsigmasupport,
                       m_opt_secondpass_halfwindowsize, m_redshiftSeparation);
  if (retInit != 0)
    THROWG(INTERNAL_ERROR, "Linemodel, init failed. Aborting");

  m_linemodel.m_opt_firstpass_fittingmethod = m_opt_firstpass_fittingmethod;
  //
  if (m_opt_continuumcomponent == "tplfit" ||
      m_opt_continuumcomponent == "tplfitauto") {
    Log.LogDetail(
        "  method Linemodel wit tplfit: fitcontinuum_maxN set to %.0f",
        m_opt_continuumfitcount);

    m_linemodel.m_opt_tplfit_fftprocessing = m_opt_tplfit_fftprocessing;
    m_linemodel.m_opt_tplfit_fftprocessing_secondpass =
        m_opt_tplfit_fftprocessing_secondpass;
    m_linemodel.m_opt_tplfit_use_photometry = m_opt_tplfit_use_photometry;
    m_linemodel.m_opt_tplfit_dustFit = Int32(m_opt_tplfit_dustfit);
    m_linemodel.m_opt_tplfit_extinction = Int32(m_opt_tplfit_igmfit);
    m_linemodel.m_opt_fitcontinuum_maxN = m_opt_continuumfitcount;
    m_linemodel.m_opt_tplfit_ignoreLinesSupport =
        Int32(m_opt_tplfit_ignoreLinesSupport);
    m_linemodel.m_opt_secondpasslcfittingmethod =
        m_opt_secondpasslcfittingmethod;
    m_linemodel.m_opt_tplfit_continuumprior_dirpath =
        m_opt_tplfit_continuumprior_dirpath;
    m_linemodel.m_opt_tplfit_continuumprior_betaA =
        m_opt_tplfit_continuumprior_betaA;
    m_linemodel.m_opt_tplfit_continuumprior_betaTE =
        m_opt_tplfit_continuumprior_betaTE;
    m_linemodel.m_opt_tplfit_continuumprior_betaZ =
        m_opt_tplfit_continuumprior_betaZ;
    m_linemodel.m_opt_continuum_neg_amp_threshold =
        m_opt_continuum_neg_amp_threshold;
    m_linemodel.m_opt_continuum_null_amp_threshold =
        m_opt_continuum_null_amp_threshold;
  }

  m_linemodel.m_opt_lya_forcefit = m_opt_lya_forcefit;
  m_linemodel.m_opt_lya_forcedisablefit = m_opt_lya_forcedisablefit;
  m_linemodel.m_opt_lya_fit_asym_min = m_opt_lya_fit_asym_min;
  m_linemodel.m_opt_lya_fit_asym_max = m_opt_lya_fit_asym_max;
  m_linemodel.m_opt_lya_fit_asym_step = m_opt_lya_fit_asym_step;
  m_linemodel.m_opt_lya_fit_width_min = m_opt_lya_fit_width_min;
  m_linemodel.m_opt_lya_fit_width_max = m_opt_lya_fit_width_max;
  m_linemodel.m_opt_lya_fit_width_step = m_opt_lya_fit_width_step;
  m_linemodel.m_opt_lya_fit_delta_min = m_opt_lya_fit_delta_min;
  m_linemodel.m_opt_lya_fit_delta_max = m_opt_lya_fit_delta_max;
  m_linemodel.m_opt_lya_fit_delta_step = m_opt_lya_fit_delta_step;

  if (m_opt_rigidity == "tplshape") {
    m_linemodel.m_opt_tplratio_ismFit = Int32(m_opt_tplratio_ismfit);
    m_linemodel.m_opt_firstpass_tplratio_ismFit =
        Int32(m_opt_firstpass_tplratio_ismfit);

    m_linemodel.m_opt_tplratio_prior_dirpath = m_opt_tplratio_prior_dirpath;
    m_linemodel.m_opt_tplratio_prior_betaA = m_opt_tplratio_prior_betaA;
    m_linemodel.m_opt_tplratio_prior_betaTE = m_opt_tplratio_prior_betaTE;
    m_linemodel.m_opt_tplratio_prior_betaZ = m_opt_tplratio_prior_betaZ;
  }

  if (m_opt_rigidity == "rules") {
    m_linemodel.m_opt_enableImproveBalmerFit = m_opt_enableImproveBalmerFit;
  }
  // logstep from redshift

  //**************************************************
  // FIRST PASS
  //**************************************************
  Int32 retFirstPass = m_linemodel.ComputeFirstPass(
      spc, rebinnedSpc, tplCatalog, tplRatioCatalog, lambdaRange, photBandCat,
      photo_weight, m_opt_fittingmethod, m_opt_lineWidthType,
      m_opt_velocity_emission, m_opt_velocity_absorption, m_opt_continuumreest,
      m_opt_rules, m_opt_velocityfit, m_opt_firstpass_largegridstepRatio,
      m_opt_firstpass_largegridsampling, m_opt_rigidity, m_opt_haPrior);
  if (retFirstPass != 0) {
    THROWG(INTERNAL_ERROR, "Linemodel, first pass failed. Aborting");
    return false;
  }

  //**************************************************
  // Compute z-candidates
  //**************************************************
  std::shared_ptr<const CLineModelResult> lmresult =
      std::dynamic_pointer_cast<const CLineModelResult>(
          m_linemodel.getResult());

  ChisquareArray chisquares = BuildChisquareArray(lmresult);

  // TODO deal with the case lmresult->Redshifts=1
  Int32 extremacount = 5;
  COperatorPdfz pdfz(m_opt_pdfcombination,
                     2 * m_opt_secondpass_halfwindowsize, // peak separation
                     m_opt_candidatesLogprobaCutThreshold, extremacount, "FPE");

  std::shared_ptr<PdfCandidatesZResult> candResult_fp =
      pdfz.Compute(chisquares, false);
  m_linemodel.SetFirstPassCandidates(candResult_fp->m_ranked_candidates);

  // save firstpass pdf
  // do the necessary to pass a pdf with 1E-3 precision only
  const std::shared_ptr<const CPdfMargZLogResult> coarsePDFZ =
      pdfz.compressFirstpassPDF(m_opt_firstpass_largegridstepRatio);
  resultStore->StoreScopedGlobalResult("firstpass_pdf", coarsePDFZ);

  //**************************************************
  // FIRST PASS + CANDIDATES - B
  //**************************************************
  bool enableFirstpass_B = (m_opt_extremacountB > 0) &&
                           (m_opt_continuumcomponent == "tplfit" ||
                            m_opt_continuumcomponent == "tplfitauto") &&
                           (m_opt_extremacountB > 1);
  COperatorLineModel linemodel_fpb;
  std::string fpb_opt_continuumcomponent =
      "fromspectrum"; // Note: this is hardocoded! given that condition for FPB
                      // relies on having "tplfit"
  Int32 retInitB =
      linemodel_fpb.Init(spc, redshifts, restLineList, m_categoryList,
                         fpb_opt_continuumcomponent, m_opt_nsigmasupport,
                         m_opt_secondpass_halfwindowsize, m_redshiftSeparation);
  if (retInitB != 0) {
    Log.LogError("Linemodel fpB, init failed. Aborting");
    return false;
  }
  if (enableFirstpass_B) {
    Log.LogInfo("Linemodel FIRST PASS B enabled. Computing now.");

    linemodel_fpb.m_opt_firstpass_fittingmethod = m_opt_firstpass_fittingmethod;

    if (fpb_opt_continuumcomponent == "tplfit" ||
        fpb_opt_continuumcomponent == "tplfitauto") {
      linemodel_fpb.m_opt_tplfit_dustFit = Int32(m_opt_tplfit_dustfit);
      linemodel_fpb.m_opt_tplfit_extinction = Int32(m_opt_tplfit_igmfit);
      Log.LogDetail("  method tplfit Linemodel: fitcontinuum_maxN set to %d",
                    m_opt_continuumfitcount);
      linemodel_fpb.m_opt_fitcontinuum_maxN = m_opt_continuumfitcount;
      linemodel_fpb.m_opt_tplfit_ignoreLinesSupport =
          Int32(m_opt_tplfit_ignoreLinesSupport);
      linemodel_fpb.m_opt_secondpasslcfittingmethod =
          m_opt_secondpasslcfittingmethod;
    }

    if (m_opt_rigidity == "tplshape") {
      linemodel_fpb.m_opt_tplratio_ismFit = Int32(m_opt_tplratio_ismfit);
      linemodel_fpb.m_opt_firstpass_tplratio_ismFit =
          Int32(m_opt_firstpass_tplratio_ismfit);
    }

    //**************************************************
    // FIRST PASS B
    //**************************************************
    Int32 retFirstPass = linemodel_fpb.ComputeFirstPass(
        spc, rebinnedSpc, tplCatalog, tplRatioCatalog, lambdaRange, photBandCat,
        photo_weight, m_opt_fittingmethod, m_opt_lineWidthType,
        m_opt_velocity_emission, m_opt_velocity_absorption,
        m_opt_continuumreest, m_opt_rules, m_opt_velocityfit,
        m_opt_firstpass_largegridstepRatio, m_opt_firstpass_largegridsampling,
        m_opt_rigidity);
    if (retFirstPass != 0) {
      Log.LogError("Linemodel, first pass failed. Aborting");
      return false;
    }

    //**************************************************
    // Compute z-candidates B
    //**************************************************
    std::shared_ptr<const CLineModelResult> lmresult =
        std::dynamic_pointer_cast<const CLineModelResult>(
            linemodel_fpb.getResult());

    ChisquareArray chisquares = BuildChisquareArray(lmresult);

    // TODO deal with the case lmresult->Redshifts=1
    Int32 extremacount = 5;
    COperatorPdfz pdfz(m_opt_pdfcombination,
                       2 * m_opt_secondpass_halfwindowsize, // peak separation
                       m_opt_candidatesLogprobaCutThreshold, extremacount,
                       "FPB");

    std::shared_ptr<PdfCandidatesZResult> candResult =
        pdfz.Compute(chisquares, false);

    linemodel_fpb.SetFirstPassCandidates(candResult->m_ranked_candidates);
    // resultStore->StoreScopedGlobalResult( "firstpassb_pdf",
    // pdfz.m_postmargZResult);

    //**************************************************
    // COMBINE CANDIDATES
    //**************************************************
    m_linemodel.Combine_firstpass_candidates(
        linemodel_fpb.m_firstpass_extremaResult);
  }

  std::shared_ptr<const LineModelExtremaResult> fpExtremaResult =
      m_linemodel.buildFirstPassExtremaResults(
          m_linemodel.m_firstpass_extremaResult->m_ranked_candidates);

  // save linemodel firstpass extrema results
  std::string firstpassExtremaResultsStr = scopeStr;
  firstpassExtremaResultsStr.append("_firstpass_extrema");
  resultStore->StoreScopedGlobalResult(firstpassExtremaResultsStr.c_str(),
                                       fpExtremaResult);

  //**************************************************
  // SECOND PASS
  //**************************************************
  if (!m_opt_skipsecondpass) {
    Int32 retSecondPass = m_linemodel.ComputeSecondPass(
        spc, rebinnedSpc, tplCatalog, lambdaRange, photBandCat, fpExtremaResult,
        photo_weight, m_opt_fittingmethod, m_opt_lineWidthType,
        m_opt_velocity_emission, m_opt_velocity_absorption,
        m_opt_continuumreest, m_opt_rules, m_opt_velocityfit, m_opt_rigidity,
        m_opt_em_velocity_fit_min, m_opt_em_velocity_fit_max,
        m_opt_em_velocity_fit_step, m_opt_abs_velocity_fit_min,
        m_opt_abs_velocity_fit_max, m_opt_abs_velocity_fit_step,
        m_opt_secondpass_continuumfit);
    if (retSecondPass != 0) {
      Log.LogError("Linemodel, second pass failed. Aborting");
      return false;
    }
  } else {
    m_linemodel.m_secondpass_parameters_extremaResult =
        *m_linemodel.m_firstpass_extremaResult;
  }

  // read it as constant to save it
  std::shared_ptr<const CLineModelResult> result =
      std::dynamic_pointer_cast<const CLineModelResult>(
          m_linemodel.getResult());

  if (!result)
    THROWG(INTERNAL_ERROR, "Failed to get linemodel result");

  // save linemodel chisquare results
  resultStore->StoreScopedGlobalResult(scopeStr.c_str(), result);

  // don't save linemodel extrema results, since will change with pdf
  // computation

  return true;
}
