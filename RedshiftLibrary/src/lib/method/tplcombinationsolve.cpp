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
#include <cfloat>

#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/method/tplcombinationsolve.h"
#include "RedshiftLibrary/method/tplcombinationsolveresult.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/operator/tplcombinationresult.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/statistics/zprior.h"

using namespace NSEpic;
using namespace std;

CTplCombinationSolve::CTplCombinationSolve()
    : CObjectSolve("tplCombinationSolve") {}

std::shared_ptr<CSolveResult> CTplCombinationSolve::compute() {
  auto const &inputContext = Context.GetInputContext();
  auto const &resultStore = Context.GetResultStore();

  const CSpectrum &spc = *(inputContext->GetSpectrum());
  const CTemplateCatalog &tplCatalog = *(inputContext->GetTemplateCatalog());

  bool storeResult = false;
  m_redshiftSeparation = inputContext->GetParameterStore()->Get<Float64>(
      "extremaRedshiftSeparation");
  m_opt_maxCandidate =
      inputContext->GetParameterStore()->GetScoped<int>("extremaCount");
  m_opt_pdfcombination =
      inputContext->GetParameterStore()->GetScoped<std::string>(
          "pdfCombination");
  std::string opt_interp =
      inputContext->GetParameterStore()->GetScoped<std::string>(
          "interpolation");
  bool opt_dustFit =
      inputContext->GetParameterStore()->GetScoped<bool>("ismFit");
  Float64 overlapThreshold =
      inputContext->GetParameterStore()->GetScoped<Float64>("overlapThreshold");
  bool opt_extinction =
      inputContext->GetParameterStore()->GetScoped<bool>("igmFit");

  std::string opt_spcComponent =
      inputContext->GetParameterStore()->GetScoped<std::string>(
          "spectrum.component");

  std::vector<CMask> maskList;

  std::string scopeStr = "tplcombination";

  EType _type = nType_raw;
  if (opt_spcComponent == "raw") {
    _type = nType_raw;
  } else if (opt_spcComponent == "noContinuum") {
    _type = nType_noContinuum;
    scopeStr = "tplcombination_nocontinuum";
  } else if (opt_spcComponent == "continuum") {
    _type = nType_continuumOnly;
    scopeStr = "tplcombination_continuum";
  } else if (opt_spcComponent == "all") {
    _type = nType_all;
  } else {
    THROWG(ErrorCode::INTERNAL_ERROR, "Unknown spectrum component");
  }

  // for now interp must be 'lin'. pfg not availbale for now...
  if (opt_interp != "lin") {
    THROWG(ErrorCode::INTERNAL_ERROR, "interpolation parameter must be 'lin'");
  }

  Log.LogInfo(Formatter() << "Method parameters:");
  Log.LogInfo(Formatter() << "    -interpolation: " << opt_interp);
  Log.LogInfo(Formatter() << "    -overlapThreshold: " << overlapThreshold);
  Log.LogInfo(Formatter() << "    -component: " << opt_spcComponent);
  Log.LogInfo(Formatter() << "    -IGM extinction: "
                          << (opt_extinction ? "true" : "false"));
  Log.LogInfo(Formatter() << "    -ISM dust-fit: "
                          << (opt_dustFit ? "true" : "false"));
  Log.LogInfo("");

  Solve(resultStore, spc, tplCatalog, m_lambdaRange, m_redshifts,
        overlapThreshold, maskList, _type, opt_interp, opt_extinction,
        opt_dustFit);

  COperatorPdfz pdfz(m_opt_pdfcombination, m_redshiftSeparation, 0.0,
                     m_opt_maxCandidate, m_redshiftSampling == "log");

  std::shared_ptr<PdfCandidatesZResult> candidateResult =
      pdfz.Compute(BuildChisquareArray(resultStore, scopeStr));

  // save in resultstore pdf results
  resultStore->StoreScopedGlobalResult("pdf", pdfz.m_postmargZResult);
  resultStore->StoreScopedGlobalResult("pdf_params", pdfz.m_postmargZResult);
  // save in resultstore candidates results
  resultStore->StoreScopedGlobalResult("candidatesresult", candidateResult);

  // for each candidate, get best model by reading from datastore and selecting
  // best fit
  /////////////////////////////////////////////////////////////////////////////////////
  TFloat64Range clampedLbdaRange;
  spc.GetSpectralAxis().ClampLambdaRange(m_lambdaRange, clampedLbdaRange);
  std::shared_ptr<const TplCombinationExtremaResult> extremaResult =
      buildExtremaResults(resultStore, scopeStr,
                          candidateResult->m_ranked_candidates, spc, tplCatalog,
                          clampedLbdaRange, overlapThreshold);
  // store extrema results
  StoreExtremaResults(resultStore, extremaResult);

  std::shared_ptr<CTplCombinationSolveResult> solveResult =
      std::make_shared<CTplCombinationSolveResult>(
          extremaResult->m_ranked_candidates[0].second, m_opt_pdfcombination,
          pdfz.m_postmargZResult->valMargEvidenceLog);

  return solveResult;
}

bool CTplCombinationSolve::Solve(
    std::shared_ptr<COperatorResultStore> resultStore, const CSpectrum &spc,
    const CTemplateCatalog &tplCatalog, const TFloat64Range &lambdaRange,
    const TFloat64List &redshifts, Float64 overlapThreshold,
    const std::vector<CMask> &maskList, EType spctype,
    const std::string &opt_interp, bool opt_extinction, bool opt_dustFitting) {
  std::string scopeStr = "tplcombination";
  Int32 _ntype = 1;
  CSpectrum::EType _spctype = CSpectrum::nType_raw;
  CSpectrum::EType _spctypetab[3] = {CSpectrum::nType_raw,
                                     CSpectrum::nType_noContinuum,
                                     CSpectrum::nType_continuumOnly};

  Int32 enable_extinction =
      0; // TODO: extinction should be deactivated for nocontinuum anyway ? TBD
  if (opt_extinction) {
    enable_extinction = 1;
  }
  Int32 enable_dustFitting = 0;
  if (opt_dustFitting) {
    enable_dustFitting =
        1; // here we dont distinguish between using on single ismCoeff or
           // iterating over all coeffs. Default to all!
  }

  // prepare the list of components/templates
  const TTemplateConstRefList &tplList =
      tplCatalog.GetTemplateList(TStringList{m_category});
  checkTemplates(tplList);

  // case: nType_all
  if (spctype == nType_all)
    _ntype = 3;

  const CSpectrum::EType save_spcType = spc.GetType();
  std::vector<CSpectrum::EType> save_tplTypes(tplList.size());
  std::transform(tplList.begin(), tplList.end(), save_tplTypes.begin(),
                 [](std::shared_ptr<const NSEpic::CTemplate> tpl) {
                   return tpl->GetType();
                 });

  for (Int32 i = 0; i < _ntype; i++) {
    if (spctype == nType_all)
      _spctype = _spctypetab[i];
    else
      _spctype = static_cast<CSpectrum::EType>(spctype);

    spc.SetType(_spctype);
    std::for_each(
        tplList.begin(), tplList.end(),
        [&_spctype, &opt_interp](std::shared_ptr<const NSEpic::CTemplate> tpl) {
          tpl->SetType(_spctype);
          tpl->setRebinInterpMethod(opt_interp);
        });

    scopeStr = getSpecBasedScope(_spctype);
    if (_spctype == CSpectrum::nType_noContinuum)
      enable_dustFitting = 0;

    // Compute merit function
    auto result = std::dynamic_pointer_cast<CTplCombinationResult>(
        m_tplcombinationOperator.Compute(
            spc, tplList, lambdaRange, redshifts, overlapThreshold, maskList,
            opt_interp, enable_extinction, enable_dustFitting));

    if (!result)
      return false;

    // Store results
    Log.LogDetail("tplcombinationsolver: Save tplcombination results");
    resultStore->StoreScopedGlobalResult(scopeStr.c_str(), result);
  }

  // restore component types
  spc.SetType(save_spcType);
  auto it = save_tplTypes.begin();
  std::for_each(tplList.begin(), tplList.end(),
                [&it](std::shared_ptr<const NSEpic::CTemplate> tpl) {
                  tpl->SetType(*it);
                  it++;
                });
  return true;
}

void CTplCombinationSolve::checkTemplates(
    const TTemplateConstRefList &tplList) const {
  if (tplList.empty())
    THROWG(ErrorCode::BAD_TEMPLATECATALOG,
           Formatter() << "Empty template catalog for category " << m_category);

  // check all templates have same spectralAxis
  const CSpectrumSpectralAxis &refSpcAxis = tplList[0]->GetSpectralAxis();
  Int32 axisSize = refSpcAxis.GetSamplesCount();
  for (const auto &tpl : tplList) {
    const CSpectrumSpectralAxis &currentSpcAxis = tpl->GetSpectralAxis();
    if (axisSize != tpl->GetSampleCount())
      THROWG(ErrorCode::INTERNAL_ERROR, "templates do not have same size");

    for (Int32 i = 0; i < axisSize; i++)
      if (std::abs(refSpcAxis[i] - currentSpcAxis[i]) > 1E-8)
        THROWG(ErrorCode::INTERNAL_ERROR,
               "templates do not have same spectralAxis");
  }
}

std::string CTplCombinationSolve::getSpecBasedScope(CSpectrum::EType _spctype) {
  if (_spctype == CSpectrum::nType_continuumOnly)
    // use continuum only
    return "tplcombination_continuum";

  if (_spctype == CSpectrum::nType_raw)
    // use full spectrum
    return "tplcombination";

  if (_spctype == CSpectrum::nType_noContinuum)
    // use spectrum without continuum
    return "tplcombination_nocontinuum";
  else
    THROWG(ErrorCode::INTERNAL_ERROR, "Unknown spectrum component");
}

ChisquareArray CTplCombinationSolve::BuildChisquareArray(
    std::shared_ptr<COperatorResultStore> store,
    const std::string &scopeStr) const {

  Log.LogDetail("tplcombinationsolver: build chisquare array");
  Log.LogDetail(Formatter()
                << "    tplcombinationsolver: using results in scope: "
                << store->GetScopedName(scopeStr));

  auto results = store->GetScopedGlobalResult(scopeStr.c_str());
  if (results.expired())
    THROWG(ErrorCode::INTERNAL_ERROR,
           "Unable to retrieve tplcombination results");

  std::shared_ptr<const CTplCombinationResult> result =
      std::dynamic_pointer_cast<const CTplCombinationResult>(results.lock());

  ChisquareArray chisquarearray;
  chisquarearray.cstLog = -1;
  chisquarearray.zstep = m_redshiftStep;
  Int32 retPdfz = -1;

  Int32 nISM = result->nISM;
  Int32 nIGM = result->nIGM;

  if (chisquarearray.cstLog == -1) {
    chisquarearray.cstLog = result->CstLog;
    Log.LogInfo(Formatter() << "tplcombinationsolver: using cstLog = "
                            << chisquarearray.cstLog);
  } else if (chisquarearray.cstLog != result->CstLog)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "cstLog values do not match in results: val1="
                       << chisquarearray.cstLog
                       << " != val2=" << result->CstLog);

  if (!chisquarearray.redshifts.size())
    chisquarearray.redshifts = result->Redshifts;

  // check chi2 results status for this template

  CZPrior zpriorhelper;
  for (Int32 kism = 0; kism < nISM; kism++) {
    for (Int32 kigm = 0; kigm < nIGM; kigm++) {
      chisquarearray.zpriors.push_back(
          zpriorhelper.GetConstantLogZPrior(result->Redshifts.size()));

      // correct chi2 for ampl. marg. if necessary: todo add switch,
      // currently deactivated
      chisquarearray.chisquares.emplace_back(
          result->ChiSquareIntermediate.size(), DBL_MAX);
      TFloat64List &logLikelihoodCorrected = chisquarearray.chisquares.back();
      for (Int32 kz = 0; kz < result->ChiSquareIntermediate.size(); kz++)
        logLikelihoodCorrected[kz] =
            result->ChiSquareIntermediate
                [kz][kism]
                [kigm]; // + resultXXX->ScaleMargCorrectionTplshapes[][]?;
      Log.LogDetail(
          Formatter()
          << "    tplcombinationsolver: Pdfz combine - prepared merit "
             "#"
          << chisquarearray.chisquares.size() - 1 << " for ism=" << kism
          << ", igm=" << kigm);
    }
  }
  return chisquarearray;
}

std::shared_ptr<const TplCombinationExtremaResult>
CTplCombinationSolve::buildExtremaResults(
    std::shared_ptr<const COperatorResultStore> store,
    const std::string &scopeStr, const TCandidateZbyRank &ranked_zCandidates,
    const CSpectrum &spc, const CTemplateCatalog &tplCatalog,
    const TFloat64Range &lambdaRange, Float64 overlapThreshold) {

  Log.LogDetail("CTplCombinationSolve::buildExtremaResults: building "
                "chisquare array");
  Log.LogDetail(Formatter()
                << "    tplCombinationSolve: using results in scope: "
                << store->GetScopedName(scopeStr));
  // in contrast to linemodel and TF, there is no perTemplateResults for
  // tplCombination
  auto results = store->GetScopedGlobalResult(scopeStr.c_str());
  if (results.expired()) {
    THROWG(ErrorCode::INTERNAL_ERROR,
           "Unable to retrieve tplcombination results");
  }
  auto TplFitResult =
      std::dynamic_pointer_cast<const CTplCombinationResult>(results.lock());
  const TFloat64List &redshifts = TplFitResult->Redshifts;

  bool foundRedshiftAtLeastOnce = false;

  if (TplFitResult->ChiSquare.size() != redshifts.size()) {
    THROWG(ErrorCode::INTERNAL_ERROR,
           "Size do not match among templatefitting results");
  }

  // prepare the list of components/templates
  const TTemplateConstRefList &tplList =
      tplCatalog.GetTemplateList(TStringList{m_category});

  std::shared_ptr<TplCombinationExtremaResult> extremaResult =
      make_shared<TplCombinationExtremaResult>(ranked_zCandidates);
  Int32 extremumCount = ranked_zCandidates.size();
  for (Int32 i = 0; i < extremumCount; i++) {
    Float64 z = ranked_zCandidates[i].second->Redshift;
    auto itZ = std::find(redshifts.begin(), redshifts.end(), z);
    const Int32 idx = std::distance(redshifts.begin(), itZ);

    // Fill extrema Result
    extremaResult->m_ranked_candidates[i].second->fittedTpl.merit =
        TplFitResult->ChiSquare[idx];
    extremaResult->m_ranked_candidates[i].second->fittedTpl.tplMeritPhot =
        TplFitResult->ChiSquarePhot[idx];
    extremaResult->m_ranked_candidates[i].second->fittedTpl.meiksinIdx =
        TplFitResult->FitMeiksinIdx[idx];
    extremaResult->m_ranked_candidates[i].second->fittedTpl.ebmvCoef =
        TplFitResult->FitEbmvCoeff[idx];
    extremaResult->m_ranked_candidates[i].second->FittedTplAmplitudeList =
        TplFitResult->FitAmplitude[idx];
    extremaResult->m_ranked_candidates[i].second->FittedTplAmplitudeErrorList =
        TplFitResult->FitAmplitudeError[idx];
    extremaResult->m_ranked_candidates[i].second->FittedTplAmplitudeSigmaList =
        TplFitResult->FitAmplitudeSigma[idx];
    extremaResult->m_ranked_candidates[i].second->FittedTplCovMatrix =
        TplFitResult->FitCOV[idx];
    extremaResult->m_ranked_candidates[i].second->fittedTpl.tplLogPrior = NAN;
    extremaResult->m_ranked_candidates[i].second->fittedTpl.SNR =
        TplFitResult->SNR[idx];
    // make sure tpl is non-rebinned
    bool currentSampling = tplCatalog.m_logsampling;
    tplCatalog.m_logsampling = false;
    std::shared_ptr<CModelSpectrumResult> spcmodelPtr =
        m_tplcombinationOperator.ComputeSpectrumModel(
            spc, tplList, z, TplFitResult->FitEbmvCoeff[idx],
            TplFitResult->FitMeiksinIdx[idx], TplFitResult->FitAmplitude[idx],
            lambdaRange, overlapThreshold);
    tplCatalog.m_logsampling = currentSampling;
    if (spcmodelPtr == nullptr)
      THROWG(ErrorCode::INTERNAL_ERROR, "Couldnt compute spectrum model");
    extremaResult->m_savedModelSpectrumResults[i] = std::move(spcmodelPtr);
  }

  return extremaResult;
}

void CTplCombinationSolve::StoreExtremaResults(
    std::shared_ptr<COperatorResultStore> resultStore,
    std::shared_ptr<const TplCombinationExtremaResult> &extremaResult) const {
  resultStore->StoreScopedGlobalResult("extrema_results", extremaResult);
  Log.LogInfo("TplCombination, saving extrema results");

  return;
}
