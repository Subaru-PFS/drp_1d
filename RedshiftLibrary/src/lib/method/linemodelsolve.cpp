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
#include <fstream>
#include <iostream>
#include <string>

#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/method/linemodelsolve.h"
#include "RedshiftLibrary/method/linemodelsolveresult.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/statistics/pdfcandidateszresult.h"
#include "RedshiftLibrary/statistics/zprior.h"

using namespace NSEpic;
using namespace std;
using namespace boost;
/**
 * \brief Empty constructor.
 **/
CLineModelSolve::CLineModelSolve() : CTwoPassSolve("lineModelSolve") {
  m_opt_skipsecondpass = false;
}

/**
 * \brief
 * Populates the method parameters from the dataStore into the class members
 * Returns true if successful, false otherwise
 **/
bool CLineModelSolve::PopulateParameters(
    std::shared_ptr<const CParameterStore> parameterStore) {

  m_opt_lineratiotype =
      parameterStore->GetScoped<std::string>("lineModel.lineRatioType");

  m_opt_continuumreest =
      parameterStore->GetScoped<std::string>("lineModel.continuumReestimation");
  m_opt_continuumcomponent = TContinuumComponent(
      parameterStore->GetScoped<std::string>("lineModel.continuumComponent"));

  m_opt_pdfcombination =
      parameterStore->GetScoped<std::string>("lineModel.pdfCombination");
  m_opt_extremacount =
      parameterStore->GetScoped<Int32>("lineModel.extremaCount");
  m_opt_extremacountB =
      parameterStore->GetScoped<Int32>("lineModel.extremaCountB");
  m_opt_maxCandidate =
      parameterStore->GetScoped<Int32>("lineModel.firstPass.extremaCount");

  m_opt_stronglinesprior =
      parameterStore->GetScoped<Float64>("lineModel.strongLinesPrior");
  m_opt_haPrior = parameterStore->GetScoped<Float64>("lineModel.hAlphaPrior");
  m_opt_euclidNHaEmittersPriorStrength =
      parameterStore->GetScoped<Float64>("lineModel.nOfZPriorStrength");

  m_opt_secondpass_halfwindowsize =
      parameterStore->GetScoped<Float64>("lineModel.secondPass.halfWindowSize");

  m_opt_candidatesLogprobaCutThreshold =
      parameterStore->GetScoped<Float64>("lineModel.extremaCutProbaThreshold");

  m_opt_firstpass_largegridstepRatio = parameterStore->GetScoped<Int32>(
      "lineModel.firstPass.largeGridStepRatio");

  m_opt_skipsecondpass =
      parameterStore->GetScoped<bool>("lineModel.skipSecondPass");

  m_useloglambdasampling =
      parameterStore->GetScoped<bool>("lineModel.useLogLambdaSampling");
  return true;
}

/**
 * \brief Calls the Solve method and returns a new "result" object.
 * Call Solve.
 * Return a pointer to an empty CLineModelSolveResult. (The results for
 *Linemodel will reside in the linemodel.linemodel result).
 **/
std::shared_ptr<CSolveResult> CLineModelSolve::compute() {
  auto const &inputContext = Context.GetInputContext();
  auto const &resultStore = Context.GetResultStore();

  PopulateParameters(inputContext->GetParameterStore());

  Solve();

  auto result = resultStore->GetScopedGlobalResult("lineModel");
  if (result.expired())
    THROWG(ErrorCode::INTERNAL_ERROR,
           "lineModelsolve: Unable to retrieve linemodel results");
  std::shared_ptr<const CLineModelResult> lmresult =
      std::dynamic_pointer_cast<const CLineModelResult>(result.lock());

  //  suggestion : CSolve::GetCurrentScopeName(CScopeStack)
  //  prepare the linemodel chisquares and prior results for pdf computation
  ChisquareArray chisquares = BuildChisquareArray(
      lmresult, m_linemodel.getSPZGridParams(),
      m_linemodel.m_firstpass_extremaResult.m_ranked_candidates);

  /*
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
      m_zLogSampling,
      "SPE", // Id_prefix
      false, // do not allow extrema at border
      1      // one peak/window only
  );

  std::shared_ptr<PdfCandidatesZResult> candidateResult =
      pdfz.Compute(chisquares);

  // store PDF results
  Log.LogInfo(Formatter() << __func__ << ": Storing PDF results");
  resultStore->StoreScopedGlobalResult("pdf", pdfz.m_postmargZResult);
  // TODO: clean below line once #7646 is solved
  resultStore->StoreScopedGlobalResult("pdf_params", pdfz.m_postmargZResult);

  // Get linemodel results at extrema (recompute spectrum model etc.)
  const CSpectrum &spc = *(inputContext->GetSpectrum(m_useloglambdasampling));
  std::shared_ptr<LineModelExtremaResult> ExtremaResult =
      m_linemodel.buildExtremaResults(
          candidateResult->m_ranked_candidates,
          m_opt_continuumreest); // maybe it's better to pass
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
    Log.LogDetail(Formatter()
                  << __func__ << ": StrongLinePresence prior enabled: factor="
                  << m_opt_stronglinesprior);
  } else {
    Log.LogDetail(Formatter()
                  << __func__ << ": StrongLinePresence prior disabled");
  }

  zPriorHaStrongestLine = (m_opt_haPrior > 0.0);
  if (zPriorHaStrongestLine) {
    Log.LogDetail(Formatter()
                  << __func__ << ": Ha strongest line prior enabled: factor="
                  << m_opt_haPrior);
  } else {
    Log.LogDetail(Formatter()
                  << __func__ << ": Ha strongest line prior disabled");
  }

  opt_nlines_snr_penalization_factor = -1;
  zPriorNLineSNR = (opt_nlines_snr_penalization_factor > 0.0);
  if (zPriorNLineSNR) {
    Log.LogDetail(Formatter()
                  << __func__ << ": N lines snr>cut prior enabled: factor="
                  << opt_nlines_snr_penalization_factor);
  } else {
    Log.LogDetail(Formatter()
                  << __func__ << ": N lines snr>cut prior disabled");
  }

  // hardcoded Euclid-NHaZprior parameter
  zPriorEuclidNHa = false;
  if (m_opt_euclidNHaEmittersPriorStrength > 0.0) {
    zPriorEuclidNHa = true;
    Log.LogDetail(Formatter()
                  << __func__
                  << ": EuclidNHa prior enabled, with strength-coeff: "
                  << m_opt_euclidNHaEmittersPriorStrength);
  } else {
    Log.LogDetail(Formatter() << __func__ << "EuclidNHa prior disabled");
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
    const TZGridListParams &spZgridParams,
    const TCandidateZbyRank &parentZCand) const {
  Log.LogDetail("LinemodelSolve: building chisquare array");

  if (m_opt_pdfcombination != "bestChi2" &&
      m_opt_pdfcombination != "bestproba" && m_opt_pdfcombination != "marg")
    THROWG(ErrorCode::BAD_PARAMETER_VALUE,
           "PdfCombination can only be {bestchi2, bestproba, marg");

  ChisquareArray chisquarearray;
  std::vector<TFloat64List> &chisquares = chisquarearray.chisquares;
  std::vector<TFloat64List> &zpriors = chisquarearray.zpriors;
  chisquarearray.zstep = m_coarseRedshiftStep;
  chisquarearray.zgridParams = spZgridParams;
  if (!spZgridParams.empty())
    chisquarearray.parentCandidates = parentZCand;

  chisquarearray.cstLog = result->cstLog;
  Log.LogDetail(Formatter()
                << __func__ << ": using cstLog = " << chisquarearray.cstLog);

  chisquarearray.redshifts = result->Redshifts;

  const Int32 zsize = result->Redshifts.size();

  if (m_opt_pdfcombination == "bestChi2") {
    zpriors.push_back(BuildZpriors(result));
    chisquares.push_back(result->ChiSquare);
    return chisquarearray;
  }

  if (m_opt_lineratiotype != "tplRatio") {
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
  Log.LogDetail(Formatter() << __func__ << ": PriorLinesTplratios.size()="
                            << result->PriorLinesTplratios.size());
  if (result->PriorLinesTplratios.size() == ntplratios) {
    zPriorLines = true;
    Log.LogDetail(Formatter() << __func__ << ": Lines Prior enabled");
  } else {
    Log.LogDetail(Formatter() << __func__ << ": Lines Prior disabled");
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

/**
 * \brief
 * Create a continuum object by subtracting spcWithoutContinuum from the spc.
 * Configure the opt_XXX variables from the dataStore scope parameters.
 * LogInfo the opt_XXX values.
 * Create a COperatorLineModel, call its Compute method.
 * If that returned true, store results.
 **/

void CLineModelSolve::Solve() {
  std::string scopeStr = "lineModel";

  std::shared_ptr<COperatorResultStore> resultStore = Context.GetResultStore();

  // Compute with linemodel operator
  m_linemodel.Init(m_redshifts, m_redshiftStep, m_zLogSampling);

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
  //   Int32 extremacount = 5;
  COperatorPdfz pdfz(m_opt_pdfcombination,
                     2 * m_opt_secondpass_halfwindowsize, // peak separation
                     m_opt_candidatesLogprobaCutThreshold, m_opt_maxCandidate,
                     m_zLogSampling, "FPE", true, 0);

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
                           m_opt_continuumcomponent.isContinuumFit() &&
                           (m_opt_extremacountB > 1);
  COperatorLineModel linemodel_fpb;
  std::string fpb_opt_continuumcomponent =
      "fromSpectrum"; // Note: this is hardocoded! given that condition for
                      // FPB relies on having "tplFit"
  linemodel_fpb.Init(m_redshifts, m_redshiftStep, m_zLogSampling);

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
    // Int32 extremacount = 5;
    COperatorPdfz pdfz(m_opt_pdfcombination,
                       2 * m_opt_secondpass_halfwindowsize, // peak separation
                       m_opt_candidatesLogprobaCutThreshold, m_opt_maxCandidate,
                       m_zLogSampling, "FPB", true, 0);

    std::shared_ptr<PdfCandidatesZResult> candResult = pdfz.Compute(chisquares);

    linemodel_fpb.SetFirstPassCandidates(candResult->m_ranked_candidates);
    // resultStore->StoreScopedGlobalResult( "firstPassb_pdf",
    // pdfz.m_postmargZResult);

    //**************************************************
    // COMBINE CANDIDATES
    //**************************************************
    m_linemodel.Combine_firstpass_candidates(
        linemodel_fpb.m_firstpass_extremaResult);
  }

  std::shared_ptr<const LineModelExtremaResult> fpExtremaResult =
      m_linemodel.BuildFirstPassExtremaResults();

  // save linemodel firstpass extrema results
  std::string firstpassExtremaResultsStr = scopeStr;
  firstpassExtremaResultsStr.append("_firstpass_extrema");
  resultStore->StoreScopedGlobalResult(firstpassExtremaResultsStr.c_str(),
                                       fpExtremaResult);

  //**************************************************
  // SECOND PASS
  //**************************************************
  if (useTwoPass())
    m_linemodel.ComputeSecondPass(fpExtremaResult);

  // read it as constant to save it
  std::shared_ptr<const CLineModelResult> result =
      std::dynamic_pointer_cast<const CLineModelResult>(
          m_linemodel.getResult());

  if (!result)
    THROWG(ErrorCode::INTERNAL_ERROR, "Failed to get linemodel result");

  // save linemodel chisquare results
  resultStore->StoreScopedGlobalResult(scopeStr.c_str(), result);

  // don't save linemodel extrema results, since will change with pdf
  // computation
}

void CLineModelSolve::initTwoPassZStepFactor() {
  m_twoPassZStepFactor =
      Context.GetInputContext()->GetParameterStore()->GetScoped<Int32>(
          "lineModelSolve.lineModel.firstPass."
          "largeGridStepRatio");
};