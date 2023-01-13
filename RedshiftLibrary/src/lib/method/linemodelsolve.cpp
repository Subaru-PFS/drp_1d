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
#include "RedshiftLibrary/method/linemodelsolveresult.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"

#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/statistics/pdfcandidateszresult.h"
#include "RedshiftLibrary/statistics/zprior.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace NSEpic;
using namespace std;
using namespace boost;
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

  m_opt_lineratiotype =
      parameterStore->GetScoped<std::string>("linemodel.lineRatioType");
  m_opt_continuumreest =
      parameterStore->GetScoped<std::string>("linemodel.continuumreestimation");
  m_opt_continuumcomponent =
      parameterStore->GetScoped<std::string>("linemodel.continuumcomponent");

  m_opt_pdfcombination =
      parameterStore->GetScoped<std::string>("linemodel.pdfcombination");
  m_opt_extremacount =
      parameterStore->GetScoped<Int32>("linemodel.extremacount");
  m_opt_extremacountB =
      parameterStore->GetScoped<Int32>("linemodel.extremacountB");

  m_opt_stronglinesprior =
      parameterStore->GetScoped<Float64>("linemodel.stronglinesprior");
  m_opt_haPrior = parameterStore->GetScoped<Float64>("linemodel.haprior");
  m_opt_euclidNHaEmittersPriorStrength =
      parameterStore->GetScoped<Float64>("linemodel.euclidnhaemittersStrength");

  m_opt_secondpass_halfwindowsize =
      parameterStore->GetScoped<Float64>("linemodel.secondpass.halfwindowsize");

  m_opt_candidatesLogprobaCutThreshold =
      parameterStore->GetScoped<Float64>("linemodel.extremacutprobathreshold");

  m_opt_firstpass_largegridstepRatio = parameterStore->GetScoped<Int32>(
      "linemodel.firstpass.largegridstepratio");

  m_opt_skipsecondpass =
      parameterStore->GetScoped<bool>("linemodel.skipsecondpass");

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
  const CSpectrum &spc = *(inputContext->GetSpectrum(false));
  PopulateParameters(inputContext->GetParameterStore());

  Solve();

  auto results = resultStore->GetScopedGlobalResult("linemodel");
  if (results.expired())
    THROWG(INTERNAL_ERROR,
           "linemodelsolve: Unable to retrieve linemodel results");

  std::shared_ptr<const CLineModelResult> result =
      std::dynamic_pointer_cast<const CLineModelResult>(results.lock());

  //  suggestion : CSolve::GetCurrentScopeName(TScopeStack)
  //  prepare the linemodel chisquares and prior results for pdf computation
  ChisquareArray chisquares =
      BuildChisquareArray(result, m_linemodel.getSPZGridParams());

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
      m_redshiftSampling == "log",
      "SPE", // Id_prefix
      false, // do not allow extrema at border
      1,     // one peak/window only
      true   // integrate under peaks
  );

  TFloat64RangeList ExtendedRedshiftsRange(
      m_linemodel.m_secondpass_parameters_extremaResult.ExtendedRedshifts
          .cbegin(),
      m_linemodel.m_secondpass_parameters_extremaResult.ExtendedRedshifts
          .cend());

  std::shared_ptr<PdfCandidatesZResult> candidateResult = pdfz.Compute(
      chisquares, ExtendedRedshiftsRange,
      m_linemodel.m_secondpass_parameters_extremaResult.m_ranked_candidates);

  // store PDF results
  Log.LogInfo("%s: Storing PDF results", __func__);
  resultStore->StoreScopedGlobalResult("pdf", pdfz.m_postmargZResult);
  // TODO: clean below line once #7646 is solved
  resultStore->StoreScopedGlobalResult("pdf_params", pdfz.m_postmargZResult);

  // Get linemodel results at extrema (recompute spectrum model etc.)
  std::shared_ptr<LineModelExtremaResult> ExtremaResult =
      m_linemodel.buildExtremaResults(
          spc, *(Context.GetClampedLambdaRange(false)),
          candidateResult->m_ranked_candidates,
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
    Int32 kTplRatio) const {
  TFloat64List zpriors;

  CZPrior zpriorhelper;

  bool zPriorStrongLinePresence, zPriorHaStrongestLine, zPriorEuclidNHa,
      zPriorNLineSNR;
  Float64 opt_nlines_snr_penalization_factor;
  GetZpriorsOptions(zPriorStrongLinePresence, zPriorHaStrongestLine,
                    zPriorNLineSNR, opt_nlines_snr_penalization_factor,
                    zPriorEuclidNHa);

  if (zPriorStrongLinePresence) {
    if (kTplRatio == -1) {
      Int32 lineTypeFilter = 1; // for emission lines only
      const TBoolList strongLinePresence = result->getStrongLinesPresence(
          lineTypeFilter, result->LineModelSolutions);
      zpriors = zpriorhelper.GetStrongLinePresenceLogZPrior(
          strongLinePresence, m_opt_stronglinesprior);
    } else {
      const TBoolList &strongLinePresence =
          result->StrongELPresentTplratios[kTplRatio];
      zpriors = zpriorhelper.GetStrongLinePresenceLogZPrior(
          strongLinePresence, m_opt_stronglinesprior);
    }

  } else {
    zpriors = zpriorhelper.GetConstantLogZPrior(result->Redshifts.size());
  }

  if (zPriorHaStrongestLine) {
    TFloat64List zlogPriorHaStrongest;
    if (kTplRatio == -1) {
      const TBoolList wHaStronglinePresence =
          result->getStrongestLineIsHa(result->LineModelSolutions);
      zlogPriorHaStrongest = zpriorhelper.GetStrongLinePresenceLogZPrior(
          wHaStronglinePresence, m_opt_haPrior);
    } else {
      const TBoolList &wHaStronglinePresence =
          result->StrongHalphaELPresentTplratios[kTplRatio];
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
    if (kTplRatio == -1) {
      const TInt32List n_lines_above_snr =
          result->getNLinesAboveSnrcut(result->LineModelSolutions);
      zlogPriorNLinesAboveSNR = zpriorhelper.GetNLinesSNRAboveCutLogZPrior(
          n_lines_above_snr, opt_nlines_snr_penalization_factor);
    } else {
      const TInt32List &n_lines_above_snr =
          result->NLinesAboveSNRTplratios[kTplRatio];
      zlogPriorNLinesAboveSNR = zpriorhelper.GetNLinesSNRAboveCutLogZPrior(
          n_lines_above_snr, opt_nlines_snr_penalization_factor);
    }
    zpriors = zpriorhelper.CombineLogZPrior(zpriors, zlogPriorNLinesAboveSNR);
  }

  return zpriors;
}

ChisquareArray CLineModelSolve::BuildChisquareArray(
    const std::shared_ptr<const CLineModelResult> &result,
    const TZGridListParams &spZgridParams) const {
  Log.LogDetail("LinemodelSolve: building chisquare array");

  if (m_opt_pdfcombination != "bestchi2" &&
      m_opt_pdfcombination != "bestproba" && m_opt_pdfcombination != "marg")
    THROWG(BAD_PARAMETER_VALUE,
           "PdfCombination can only be {bestchi2, bestproba, marg");

  ChisquareArray chisquarearray;
  std::vector<TFloat64List> &chisquares = chisquarearray.chisquares;
  std::vector<TFloat64List> &zpriors = chisquarearray.zpriors;
  chisquarearray.zstep = m_coarseRedshiftStep;
  chisquarearray.zgridParams = spZgridParams;

  chisquarearray.cstLog = result->cstLog;
  Log.LogDetail("%s: using cstLog = %f", __func__, chisquarearray.cstLog);

  chisquarearray.redshifts = result->Redshifts;

  const Int32 zsize = result->Redshifts.size();

  if (m_opt_pdfcombination == "bestchi2") {
    zpriors.push_back(BuildZpriors(result));
    chisquares.push_back(result->ChiSquare);
    return chisquarearray;
  }

  if (m_opt_lineratiotype != "tplratio") {
    zpriors.push_back(BuildZpriors(result));
    chisquares.push_back(result->ChiSquare);
  } else {
    fillChisquareArrayForTplRatio(result, chisquarearray);
  }

  if (result->ChiSquareTplContinuum.empty())
    return chisquarearray;

  // Fullmodel (ie with continuum template fitting): store all continuum
  // tpl fitting chisquares (ChiSquareTplContinuum size will be null if
  // not tplfit)

  //  note: the continuum xi2 are computed on continuum ~orthogonal to the
  //  linemodel, ie corresponding
  //        to a fullmodel with perfect linemodel fit (null Xi2 under the
  //        lines). They are thus better (smaller) than the fullmodel Xi2
  //        with which they will be summed in the marginalization. To
  //        mitigate this, for each continuum Xi2, we have to add the Xi2
  //        part of the linemodel alone. The latter can be obtained by
  //        subtracting the corresponding continuum Xi2 of the fullmodel.

  const auto newsize = chisquares.size() * result->ChiSquareTplContinuum.size();
  chisquares.reserve(newsize);
  zpriors.reserve(newsize);
  chisquarearray.modelpriors.reserve(newsize);

  // divide model prior by the number of continuum templates
  // TODO: need to add/handle tpl continuum priors
  // (note: in the case of ATEZ, which is exclusive of zpriors and
  // Tplratios, priors are included already in the chisquares of both
  // tplcontinuum chi2 and tplratio chi2)
  Int32 n = result->ChiSquareTplContinuum.size();
  std::transform(chisquarearray.modelpriors.begin(),
                 chisquarearray.modelpriors.end(),
                 chisquarearray.modelpriors.begin(),
                 [n](Float64 prior) { return prior / n; });

  // loop on all tplratios (or 1 linemodel free)
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
      for (Int32 iz = 0; iz < zsize; ++iz) {
        // estimate chi2 of fullmodel with other continuum
        fullmodel_chi2[iz] =
            chisquares[k][iz] +
            std::max(0., (*it)[iz] - result->ChiSquareTplContinuum[0][iz]);
      }
      if (!chisquarearray.modelpriors.empty())
        chisquarearray.modelpriors.push_back(
            chisquarearray.modelpriors[k]); // duplicate modelpriors
    }
  }

  return chisquarearray;
}

void CLineModelSolve::fillChisquareArrayForTplRatio(
    const std::shared_ptr<const CLineModelResult> &result,
    ChisquareArray &chisquarearray) const {

  std::vector<TFloat64List> &chisquares = chisquarearray.chisquares;
  std::vector<TFloat64List> &zpriors = chisquarearray.zpriors;

  const Int32 zsize = result->Redshifts.size();

  const Int32 ntplratios = result->ChiSquareTplratios.size();

  bool zPriorLines = false;
  Log.LogDetail("%s: PriorLinesTplratios.size()=%d", __func__,
                result->PriorLinesTplratios.size());
  if (result->PriorLinesTplratios.size() == ntplratios) {
    zPriorLines = true;
    Log.LogDetail("%s: Lines Prior enabled", __func__);
  } else {
    Log.LogDetail("%s: Lines Prior disabled", __func__);
  }

  chisquares.reserve(ntplratios);
  zpriors.reserve(ntplratios);

  for (Int32 k = 0; k < ntplratios; ++k) {
    zpriors.push_back(BuildZpriors(result, k));
    chisquares.push_back(result->ChiSquareTplratios[k]);

    // correct chi2 if necessary
    TFloat64List &logLikelihoodCorrected = chisquares.back();
    /*
    if(m_opt_pdf_margAmpCorrection) //nb: this is experimental.
    {
        //find max scalemargcorr
        Float64 maxscalemargcorr=-DBL_MAX;
        for ( Int32 kz=0; kz<zsize; kz++ )
            if(maxscalemargcorr <
    result->ScaleMargCorrectionTplratios[k][kz]) maxscalemargcorr =
    result->ScaleMargCorrectionTplratios[k][kz];

        Log.LogDetail("%s: maxscalemargcorr= %e", __func__,
    maxscalemargcorr); for ( Int32 kz=0; kz<zsize; kz++ )
            if(result->ScaleMargCorrectionTplratios[k][kz]!=0) //warning,
    this is experimental. logLikelihoodCorrected[kz] +=
    result->ScaleMargCorrectionTplratios[k][kz] - maxscalemargcorr;

        // need to add maxscalemargcorr ?
    }*/
    if (!zPriorLines || result->PriorLinesTplratios[k].size() != zsize)
      continue;

    for (Int32 kz = 0; kz < zsize; kz++)
      logLikelihoodCorrected[kz] += result->PriorLinesTplratios[k][kz];
  }

  chisquarearray.modelpriors = result->PriorTplratios;
  return;
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

void CLineModelSolve::StoreChisquareTplRatioResults(
    std::shared_ptr<COperatorResultStore> resultStore,
    std::shared_ptr<const CLineModelResult> result) const {
  for (Int32 km = 0; km < result->ChiSquareTplratios.size(); km++) {
    std::shared_ptr<CLineModelResult> result_chisquaretplratio =
        std::shared_ptr<CLineModelResult>(new CLineModelResult());
    result_chisquaretplratio->Init(result->Redshifts, result->restLineList, 0,
                                   0, TFloat64List());
    for (Int32 kz = 0; kz < result->Redshifts.size(); kz++) {
      result_chisquaretplratio->ChiSquare[kz] =
          result->ChiSquareTplratios[km][kz];
    }

    std::string resname =
        (boost::format(
             "linemodel_chisquaretplratio/linemodel_chisquaretplratio_%d") %
         km)
            .str();
    resultStore->StoreScopedGlobalResult(resname.c_str(),
                                         result_chisquaretplratio);
  }

  // Save scaleMargCorrTplratio results
  for (Int32 km = 0; km < result->ScaleMargCorrectionTplratios.size(); km++) {
    std::shared_ptr<CLineModelResult> result_chisquaretplratio =
        std::shared_ptr<CLineModelResult>(new CLineModelResult());
    result_chisquaretplratio->Init(result->Redshifts, result->restLineList, 0,
                                   0, TFloat64List());
    for (Int32 kz = 0; kz < result->Redshifts.size(); kz++) {
      result_chisquaretplratio->ChiSquare[kz] =
          result->ScaleMargCorrectionTplratios[km][kz];
    }

    std::string resname = (boost::format("linemodel_chisquaretplratio/"
                                         "linemodel_scalemargcorrtplratio_%d") %
                           km)
                              .str();
    resultStore->StoreScopedGlobalResult(resname.c_str(),
                                         result_chisquaretplratio);
  }

  // Save PriorLinesTplratios results
  for (Int32 km = 0; km < result->PriorLinesTplratios.size(); km++) {
    std::shared_ptr<CLineModelResult> result_chisquaretplratio =
        std::shared_ptr<CLineModelResult>(new CLineModelResult());
    result_chisquaretplratio->Init(result->Redshifts, result->restLineList, 0,
                                   0, TFloat64List());
    for (Int32 kz = 0; kz < result->Redshifts.size(); kz++) {
      result_chisquaretplratio->ChiSquare[kz] =
          result->PriorLinesTplratios[km][kz];
    }

    std::string resname =
        (boost::format(
             "linemodel_chisquaretplratio/linemodel_priorlinestplratio_%d") %
         km)
            .str();
    resultStore->StoreScopedGlobalResult(resname.c_str(),
                                         result_chisquaretplratio);
  }

  // Save PriorContinuumTplratios results
  std::shared_ptr<CLineModelResult> result_chisquaretplratio =
      std::shared_ptr<CLineModelResult>(new CLineModelResult());
  result_chisquaretplratio->Init(result->Redshifts, result->restLineList, 0, 0,
                                 TFloat64List());
  for (Int32 kz = 0; kz < result->Redshifts.size(); kz++) {
    result_chisquaretplratio->ChiSquare[kz] =
        result->ContinuumModelSolutions[kz].tplLogPrior;
  }

  std::string resname =
      (boost::format(
           "linemodel_chisquaretplratio/linemodel_priorcontinuumtplshape"))
          .str();
  resultStore->StoreScopedGlobalResult(resname.c_str(),
                                       result_chisquaretplratio);
}

/**
 * \brief
 * Create a continuum object by subtracting spcWithoutContinuum from the spc.
 * Configure the opt_XXX variables from the dataStore scope parameters.
 * LogInfo the opt_XXX values.
 * Create a COperatorLineModel, call its Compute method.
 * If that returned true, store results.
 **/

void CLineModelSolve::Solve() {
  std::string scopeStr = "linemodel";

  std::shared_ptr<COperatorResultStore> resultStore = Context.GetResultStore();

  // Compute with linemodel operator
  m_linemodel.Init(m_redshifts, m_redshiftStep, m_redshiftSampling);

  // logstep from redshift

  //**************************************************
  // FIRST PASS
  //**************************************************
  m_linemodel.ComputeFirstPass();

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
                     m_opt_candidatesLogprobaCutThreshold, extremacount,
                     m_redshiftSampling == "log", "FPE", true, 0, false);

  std::shared_ptr<PdfCandidatesZResult> candResult_fp =
      pdfz.Compute(chisquares);
  m_linemodel.SetFirstPassCandidates(candResult_fp->m_ranked_candidates);

  resultStore->StoreScopedGlobalResult("firstpass_pdf", pdfz.m_postmargZResult);
  // had to duplicate it to allow access from hdf5
  resultStore->StoreScopedGlobalResult("firstpass_pdf_params",
                                       pdfz.m_postmargZResult);

  //**************************************************
  // FIRST PASS + CANDIDATES - B
  //**************************************************
  bool enableFirstpass_B = (m_opt_extremacountB > 0) &&
                           (m_opt_continuumcomponent == "tplfit" ||
                            m_opt_continuumcomponent == "tplfitauto") &&
                           (m_opt_extremacountB > 1);
  COperatorLineModel linemodel_fpb;
  std::string fpb_opt_continuumcomponent =
      "fromspectrum"; // Note: this is hardocoded! given that condition for
                      // FPB relies on having "tplfit"
  linemodel_fpb.Init(m_redshifts, m_redshiftStep, m_redshiftSampling);

  if (enableFirstpass_B) {
    Log.LogInfo("Linemodel FIRST PASS B enabled. Computing now.");

    //**************************************************
    // FIRST PASS B
    //**************************************************
    linemodel_fpb.ComputeFirstPass();

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
                       m_redshiftSampling == "log", "FPB", true, 0, false);

    std::shared_ptr<PdfCandidatesZResult> candResult = pdfz.Compute(chisquares);

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
    m_linemodel.ComputeSecondPass(fpExtremaResult);
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
}

void CLineModelSolve::createRedshiftGrid(
    const std::shared_ptr<const CInputContext> &inputContext,
    const TFloat64Range &redshiftRange) {

  Int32 opt_twosteplargegridstep_ratio =
      inputContext->GetParameterStore()->GetScoped<Int32>(
          "LineModelSolve.linemodel.firstpass.largegridstepratio");

  m_coarseRedshiftStep = m_redshiftStep * opt_twosteplargegridstep_ratio;

  CZGridParam zp(redshiftRange, m_coarseRedshiftStep);
  m_redshifts = zp.getZGrid(m_redshiftSampling == "log");

  if (m_redshifts.size() < MIN_GRID_COUNT) {
    m_coarseRedshiftStep = m_redshiftStep;
    CObjectSolve::createRedshiftGrid(
        inputContext, redshiftRange); // fall back to creating fine grid
    Log.LogInfo("Operator-Linemodel: 1st pass coarse zgrid auto disabled: "
                "raw %d redshifts will be calculated",
                m_redshifts.size());
  } else {
    Log.LogInfo(
        "Operator-Linemodel: 1st pass coarse zgrid enabled: %d redshifts "
        "will be calculated on the coarse grid",
        m_redshifts.size());
  }
}
