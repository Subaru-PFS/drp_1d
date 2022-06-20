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
#include "RedshiftLibrary/method/tplcombinationsolve.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/operator/tplcombinationresult.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/statistics/zprior.h"

#include <cfloat>

using namespace NSEpic;
using namespace std;

CTplcombinationSolve::CTplcombinationSolve(TScopeStack &scope,
                                           string objectType)
    : CObjectSolve("TplcombinationSolve", scope, objectType) {}

std::shared_ptr<CSolveResult>
CTplcombinationSolve::compute(std::shared_ptr<const CInputContext> inputContext,
                              std::shared_ptr<COperatorResultStore> resultStore,
                              TScopeStack &scope)

{

  const CSpectrum &spc = *(inputContext->GetSpectrum());
  const CTemplateCatalog &tplCatalog = *(inputContext->GetTemplateCatalog());

  bool storeResult = false;
  m_redshiftSeparation = inputContext->GetParameterStore()->Get<Float64>(
      "extremaredshiftseparation");
  m_opt_maxCandidate =
      inputContext->GetParameterStore()->GetScoped<int>("extremacount");
  m_opt_pdfcombination =
      inputContext->GetParameterStore()->GetScoped<std::string>(
          "pdfcombination");
  std::string opt_interp =
      inputContext->GetParameterStore()->GetScoped<std::string>(
          "interpolation");
  bool opt_dustFit =
      inputContext->GetParameterStore()->GetScoped<bool>("dustfit");
  Float64 overlapThreshold =
      inputContext->GetParameterStore()->GetScoped<Float64>("overlapThreshold");
  bool opt_extinction =
      inputContext->GetParameterStore()->GetScoped<bool>("extinction");

  std::string opt_spcComponent =
      inputContext->GetParameterStore()->GetScoped<std::string>(
          "spectrum.component");

  std::vector<CMask> maskList;

  std::string scopeStr = "tplcombination";

  EType _type = nType_raw;
  if (opt_spcComponent == "raw") {
    _type = nType_raw;
  } else if (opt_spcComponent == "nocontinuum") {
    _type = nType_noContinuum;
    scopeStr = "tplcombination_nocontinuum";
  } else if (opt_spcComponent == "continuum") {
    _type = nType_continuumOnly;
    scopeStr = "tplcombination_continuum";
  } else if (opt_spcComponent == "all") {
    _type = nType_all;
  } else {
    THROWG(INTERNAL_ERROR, "Unknown spectrum component");
  }

  // for now interp must be 'lin'. pfg not availbale for now...
  if (opt_interp != "lin") {
    THROWG(INTERNAL_ERROR, "interpolation parameter must be 'lin'");
  }

  Log.LogInfo("Method parameters:");
  Log.LogInfo("    -interpolation: %s", opt_interp.c_str());
  Log.LogInfo("    -overlapThreshold: %.3f", overlapThreshold);
  Log.LogInfo("    -component: %s", opt_spcComponent.c_str());
  Log.LogInfo("    -IGM extinction: %s", opt_extinction ? "true" : "false");
  Log.LogInfo("    -ISM dust-fit: %s", opt_dustFit ? "true" : "false");
  // Log.LogInfo( "    -pdfcombination: %s", m_opt_pdfcombination.c_str());
  Log.LogInfo("");

  Solve(resultStore, spc, tplCatalog, m_lambdaRange, m_redshifts,
        overlapThreshold, maskList, _type, opt_interp, opt_extinction,
        opt_dustFit);

  COperatorPdfz pdfz(m_opt_pdfcombination, m_redshiftSeparation, 0.0,
                     m_opt_maxCandidate);

  std::shared_ptr<PdfCandidatesZResult> candidateResult =
      pdfz.Compute(BuildChisquareArray(resultStore, scopeStr));

  // save in resultstore pdf results
  resultStore->StoreScopedGlobalResult(
      "pdf", pdfz.m_postmargZResult); // need to store this pdf with this exact
                                      // same name so that zqual can load it.
                                      // see zqual.cpp/ExtractFeaturesPDF

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
                          clampedLbdaRange, overlapThreshold, opt_interp);
  // store extrema results
  StoreExtremaResults(resultStore, extremaResult);

  std::shared_ptr<CTplCombinationSolveResult> solveResult =
      std::make_shared<CTplCombinationSolveResult>(
          resultStore->GetCurrentScopeName(),
          extremaResult->m_ranked_candidates[0].second, m_opt_pdfcombination,
          pdfz.m_postmargZResult->valMargEvidenceLog);

  return solveResult;
}

bool CTplcombinationSolve::Solve(
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
  if (m_categoryList.size() > 1) {
    THROWG(INTERNAL_ERROR, "Multiple categories are passed for "
                           "tplcombinationsolve. Only one is required");
  }

  const TTemplateConstRefList &tplList =
      tplCatalog.GetTemplateList(m_categoryList);
  if (tplList.empty()) {
    THROWG(BAD_TEMPLATECATALOG, Formatter()
                                    << "Empty template catalog for category "
                                    << m_categoryList[0]);
  }

  // check all templates have same spectralAxis
  const CSpectrumSpectralAxis &refSpcAxis = tplList[0]->GetSpectralAxis();
  Int32 axisSize = refSpcAxis.GetSamplesCount();
  for (Int32 ktpl = 1; ktpl < tplList.size(); ktpl++) {
    const CSpectrumSpectralAxis &currentSpcAxis =
        tplList[ktpl]->GetSpectralAxis();
    if (axisSize != tplList[ktpl]->GetSampleCount()) {
      THROWG(INTERNAL_ERROR, "templates do not have same size");
    }
    for (Int32 i = 0; i < axisSize; i++) {
      if (std::abs(refSpcAxis[i] - currentSpcAxis[i]) > 1E-8) {
        THROWG(INTERNAL_ERROR, "templates do not have same spectralAxis");
      }
    }
  }

  // case: nType_all
  if (spctype == nType_all) {
    _ntype = 3;
  }

  const CSpectrum::EType save_spcType = spc.GetType();

  std::vector<CSpectrum::EType> save_tplTypes;
  for (std::shared_ptr<const NSEpic::CTemplate> tpl : tplList) {
    save_tplTypes.push_back(tpl->GetType());
  }

  for (Int32 i = 0; i < _ntype; i++) {
    if (spctype == nType_all) {
      _spctype = _spctypetab[i];
    } else {
      _spctype = static_cast<CSpectrum::EType>(spctype);
    }

    spc.SetType(_spctype);
    for (std::shared_ptr<const NSEpic::CTemplate> tpl : tplList)
      tpl->SetType(_spctype);

    if (_spctype == CSpectrum::nType_continuumOnly) {
      // use continuum only
      scopeStr = "tplcombination_continuum";

    } else if (_spctype == CSpectrum::nType_raw) {
      // use full spectrum
      scopeStr = "tplcombination";

    } else if (_spctype == CSpectrum::nType_noContinuum) {
      // use spectrum without continuum
      scopeStr = "tplcombination_nocontinuum";
      enable_dustFitting = 0;
    } else {
      THROWG(INTERNAL_ERROR, "Unknown spectrum component");
    }

    // Compute merit function
    auto result = std::dynamic_pointer_cast<CTplCombinationResult>(
        m_tplcombinationOperator.Compute(
            spc, tplList, lambdaRange, redshifts, overlapThreshold, maskList,
            opt_interp, enable_extinction, enable_dustFitting));

    if (!result) {
      // Log.LogError( "Failed to compute chi square value");
      return false;
    } else {
      // Store results
      Log.LogDetail("tplcombinationsolve: Save tplcombination results");
      resultStore->StoreScopedGlobalResult(scopeStr.c_str(), result);
    }
  }

  // restore component types
  spc.SetType(save_spcType);
  auto it = save_tplTypes.begin();
  for (std::shared_ptr<const NSEpic::CTemplate> tpl : tplList) {
    tpl->SetType(*it);
    it++;
  }

  return true;
}

ChisquareArray CTplcombinationSolve::BuildChisquareArray(
    std::shared_ptr<COperatorResultStore> store,
    const std::string &scopeStr) const {
  ChisquareArray chisquarearray;

  Log.LogDetail("tplcombinationsolve: build chisquare array");
  Log.LogDetail(Formatter()
                << "    tplcombinationsolve: using results in scope: "
                << store->GetScopedName(scopeStr));

  auto results = store->GetScopedGlobalResult(scopeStr.c_str());
  if (results.expired()) {
    THROWG(INTERNAL_ERROR, "Unable to retrieve tplcombination results");
  }
  std::shared_ptr<const CTplCombinationResult> result =
      std::dynamic_pointer_cast<const CTplCombinationResult>(results.lock());

  chisquarearray.cstLog = -1;

  Int32 retPdfz = -1;
  {
    Int32 nISM = -1;
    Int32 nIGM = -1;
    if (result->ChiSquareIntermediate.size() > 0) {
      nISM = result->ChiSquareIntermediate[0].size();
      if (result->ChiSquareIntermediate[0].size() > 0) {
        nIGM = result->ChiSquareIntermediate[0][0].size();
      }
    }
    if (chisquarearray.cstLog == -1) {
      chisquarearray.cstLog = result->CstLog;
      Log.LogInfo("tplcombinationsolve: using cstLog = %f",
                  chisquarearray.cstLog);
    } else if (chisquarearray.cstLog != result->CstLog) {
      THROWG(INTERNAL_ERROR,
             Formatter() << "cstLog values do not match in results: val1="
                         << chisquarearray.cstLog
                         << " != val2=" << result->CstLog);
    }
    if (chisquarearray.redshifts.size() == 0) {
      chisquarearray.redshifts = result->Redshifts;
    }

    // check chi2 results status for this template
    {
      bool foundBadStatus = 0;
      for (Int32 kz = 0; kz < result->Redshifts.size(); kz++) {
        if (result->Status[kz] != COperator::nStatus_OK) {
          foundBadStatus = 1;
          break;
        }
      }
      if (foundBadStatus) {
        THROWG(INTERNAL_ERROR, " Found bad status result");
      }
    }

    CZPrior zpriorhelper;
    for (Int32 kism = 0; kism < nISM; kism++) {
      for (Int32 kigm = 0; kigm < nIGM; kigm++) {
        chisquarearray.zpriors.push_back(
            zpriorhelper.GetConstantLogZPrior(result->Redshifts.size()));

        // correct chi2 for ampl. marg. if necessary: todo add switch, currently
        // deactivated
        chisquarearray.chisquares.emplace_back(
            result->ChiSquareIntermediate.size(), DBL_MAX);
        TFloat64List &logLikelihoodCorrected = chisquarearray.chisquares.back();
        for (Int32 kz = 0; kz < result->Redshifts.size(); kz++) {
          logLikelihoodCorrected[kz] =
              result->ChiSquareIntermediate
                  [kz][kism]
                  [kigm]; // + resultXXX->ScaleMargCorrectionTplshapes[][]?;
        }
        Log.LogDetail("    tplcombinationsolve: Pdfz combine - prepared merit "
                      "#%d for ism=%d, igm=%d",
                      chisquarearray.chisquares.size() - 1, kism, kigm);
      }
    }
  }

  return chisquarearray;
}
std::shared_ptr<const TplCombinationExtremaResult>
CTplcombinationSolve::buildExtremaResults(
    std::shared_ptr<const COperatorResultStore> store,
    const std::string &scopeStr, const TCandidateZbyRank &ranked_zCandidates,
    const CSpectrum &spc, const CTemplateCatalog &tplCatalog,
    const TFloat64Range &lambdaRange, Float64 overlapThreshold,
    std::string opt_interp) {

  Log.LogDetail(
      "CTplCombinationSolve::buildExtremaResults: building chisquare array");
  Log.LogDetail(Formatter()
                << "    tplCombinationSolve: using results in scope: "
                << store->GetScopedName(scopeStr));
  // in contrast to linemodel and TF, there is no perTemplateResults for
  // tplCombination
  auto results = store->GetScopedGlobalResult(scopeStr.c_str());
  if (results.expired()) {
    THROWG(INTERNAL_ERROR, "Unable to retrieve tplcombination results");
  }
  auto TplFitResult =
      std::dynamic_pointer_cast<const CTplCombinationResult>(results.lock());
  const TFloat64List &redshifts = TplFitResult->Redshifts;

  bool foundRedshiftAtLeastOnce = false;

  if (TplFitResult->ChiSquare.size() != redshifts.size()) {
    THROWG(INTERNAL_ERROR, "Size do not match among templatefitting results");
  }

  bool foundBadStatus = false;

  if (foundBadStatus) {
    THROWG(INTERNAL_ERROR, "Bad status result");
  }

  // prepare the list of components/templates
  if (m_categoryList.size() > 1) {
    THROWG(INTERNAL_ERROR, "Multiple categories are passed for "
                           "tplcombinationsolve. Only one is required");
  }

  const TTemplateConstRefList &tplList =
      tplCatalog.GetTemplateList(m_categoryList);

  std::shared_ptr<TplCombinationExtremaResult> extremaResult =
      make_shared<TplCombinationExtremaResult>(ranked_zCandidates);
  Int32 extremumCount = ranked_zCandidates.size();
  for (Int32 i = 0; i < extremumCount; i++) {
    Float64 z = ranked_zCandidates[i].second->Redshift;
    auto itZ = std::find(redshifts.begin(), redshifts.end(), z);
    const Int32 idx = std::distance(redshifts.begin(), itZ);

    // Fill extrema Result
    extremaResult->m_ranked_candidates[i].second->FittedTplMerit =
        TplFitResult->ChiSquare[idx];
    extremaResult->m_ranked_candidates[i].second->FittedTplMeritPhot =
        TplFitResult->ChiSquarePhot[idx];
    extremaResult->m_ranked_candidates[i].second->FittedTplMeiksinIdx =
        TplFitResult->FitMeiksinIdx[idx];
    extremaResult->m_ranked_candidates[i].second->FittedTplEbmvCoeff =
        TplFitResult->FitEbmvCoeff[idx];
    extremaResult->m_ranked_candidates[i].second->FittedTplAmplitudeList =
        TplFitResult->FitAmplitude[idx];
    extremaResult->m_ranked_candidates[i].second->FittedTplAmplitudeErrorList =
        TplFitResult->FitAmplitudeError[idx];
    extremaResult->m_ranked_candidates[i].second->FittedTplAmplitudeSigmaList =
        TplFitResult->FitAmplitudeSigma[idx];
    extremaResult->m_ranked_candidates[i].second->FittedTplCovMatrix =
        TplFitResult->FitCOV[idx];
    extremaResult->m_ranked_candidates[i].second->FittedTplLogPrior = NAN;
    extremaResult->m_ranked_candidates[i].second->FittedTplSNR =
        TplFitResult->SNR[idx];
    // make sure tpl is non-rebinned
    bool currentSampling = tplCatalog.m_logsampling;
    tplCatalog.m_logsampling = false;
    std::shared_ptr<CModelSpectrumResult> spcmodelPtr =
        m_tplcombinationOperator.ComputeSpectrumModel(
            spc, tplList, z, TplFitResult->FitEbmvCoeff[idx],
            TplFitResult->FitMeiksinIdx[idx], TplFitResult->FitAmplitude[idx],
            opt_interp, lambdaRange, overlapThreshold);
    tplCatalog.m_logsampling = currentSampling;
    if (spcmodelPtr == nullptr)
      THROWG(INTERNAL_ERROR, "Couldnt compute spectrum model");
    extremaResult->m_savedModelSpectrumResults[i] = std::move(spcmodelPtr);
  }

  return extremaResult;
}

void CTplcombinationSolve::StoreExtremaResults(
    std::shared_ptr<COperatorResultStore> resultStore,
    std::shared_ptr<const TplCombinationExtremaResult> &extremaResult) const {
  resultStore->StoreScopedGlobalResult("extrema_results", extremaResult);
  Log.LogInfo("TplCombination, saving extrema results");

  return;
}
