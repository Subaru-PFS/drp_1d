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
#include "RedshiftLibrary/method/templatefittingsolve.h"
#include "RedshiftLibrary/method/templatefittingsolveresult.h"

#include "RedshiftLibrary/operator/modelspectrumresult.h"

#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/operator/modelphotvalueresult.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/operator/templatefitting.h"
#include "RedshiftLibrary/operator/templatefittinglog.h"
#include "RedshiftLibrary/operator/templatefittingwithphot.h"
#include "RedshiftLibrary/photometry/photometricdata.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/statistics/zprior.h"

using namespace NSEpic;
using namespace std;

CTemplateFittingSolve::CTemplateFittingSolve()
    : CTwoPassSolve("templateFittingSolve") {}

const std::unordered_map<CTemplateFittingSolve::EType, CSpectrum::EType>
    CTemplateFittingSolve::fittingTypeToSpectrumType = {
        {EType::raw, CSpectrum::EType::raw},
        {EType::noContinuum, CSpectrum::EType::noContinuum},
        {EType::continuumOnly, CSpectrum::EType::continuumOnly},
};

void CTemplateFittingSolve::PopulateParameters(
    const std::shared_ptr<const CParameterStore> &parameterStore) {
  m_redshiftSeparation =
      parameterStore->Get<Float64>("extremaRedshiftSeparation");
  m_overlapThreshold = parameterStore->GetScoped<Float64>("overlapThreshold");
  m_spcComponent = parameterStore->GetScoped<std::string>("spectrum.component");
  m_interpolation = parameterStore->GetScoped<std::string>("interpolation");
  m_extinction = parameterStore->GetScoped<bool>("igmFit");
  m_dustFit = parameterStore->GetScoped<bool>("ismFit");

  m_fftProcessing = parameterStore->GetScoped<bool>("fftProcessing");

  if (parameterStore->HasScoped<bool>("enablePhotometry")) {
    m_usePhotometry = parameterStore->GetScoped<bool>("enablePhotometry");
    if (m_usePhotometry)
      m_photometryWeight =
          parameterStore->GetScoped<Float64>("photometry.weight");
  }
  m_opt_maxCandidate = parameterStore->GetScoped<int>("extremaCount");
  m_opt_pdfcombination =
      parameterStore->GetScoped<std::string>("pdfCombination");
  m_opt_skipsecondpass = parameterStore->GetScoped<bool>("skipSecondPass");
  if (useTwoPass()) {
    m_secondPassContinuumFit = str2ContinuumFit.at(
        parameterStore->GetScoped<std::string>("secondPass.continuumFit"));

    m_secondPass_halfwindowsize =
        parameterStore->GetScoped<Float64>("secondPass.halfWindowSize");
  }
}

void CTemplateFittingSolve::UpdateParamsAndCatalogForFft(
    const std::shared_ptr<const CInputContext> &inputContext,
    const CTemplateCatalog &tplCatalog) {
  if (m_fftProcessing && m_usePhotometry)
    THROWG(ErrorCode::FFT_WITH_PHOTOMETRY_NOTIMPLEMENTED,
           "fftProcessing not "
           "implemented with photometry enabled");

  if (m_fftProcessing) {
    m_continuumFittingOperator =
        std::make_shared<COperatorTemplateFittingLog>(m_redshifts);
    tplCatalog.m_logsampling = true;
  } else {
    if (m_usePhotometry) {
      m_continuumFittingOperator =
          std::make_shared<COperatorTemplateFittingPhot>(
              inputContext->GetPhotBandCatalog(), m_photometryWeight,
              m_redshifts);
    } else {
      m_continuumFittingOperator =
          std::make_shared<COperatorTemplateFitting>(m_redshifts);
    }
    tplCatalog.m_logsampling = false;
  }
}

void CTemplateFittingSolve::LogParameters() {
  Log.LogInfo("Method parameters:");
  Log.LogInfo(Formatter() << "    -m_overlapThreshold: " << m_overlapThreshold);
  Log.LogInfo(Formatter() << "    -component: " << m_spcComponent);
  Log.LogInfo(Formatter() << "    -interp: " << m_interpolation);
  Log.LogInfo(Formatter() << "    -IGM extinction: "
                          << (m_extinction ? "true" : "false"));
  Log.LogInfo(Formatter() << "    -ISM dust-fit: "
                          << (m_dustFit ? "true" : "false"));
  Log.LogInfo(Formatter() << "    -pdfcombination: " << m_opt_pdfcombination);
  Log.LogInfo("");
}

void CTemplateFittingSolve::CheckTemplateCatalog(
    const CTemplateCatalog &tplCatalog) {
  if (tplCatalog.GetTemplateCount(m_category) == 0) {
    THROWG(ErrorCode::BAD_TEMPLATECATALOG,
           Formatter() << "Empty template catalog for category "
                       << m_category[0]);
  }
}

std::string
CTemplateFittingSolve::MakeScopeStrFromTypeAndPass(EType type,
                                                   bool isFirstPass) const {
  std::unordered_map<EType, std::string> spectrumTypeToStr = {
      {EType::raw, "templatefitting"},
      {EType::continuumOnly, "templatefitting_continuum"},
      {EType::noContinuum, "templatefitting_nocontinuum"},
      {EType::all, "templatefitting"},
  };

  std::string scopeStr = spectrumTypeToStr.at(type);
  if (useTwoPass() && isFirstPass)
    scopeStr += "_firstpass";

  return scopeStr;
}

CTemplateFittingSolve::EType
CTemplateFittingSolve::MakeTypeFromSpectrumComponentParam(
    std::string component) {
  std::unordered_map<std::string, EType> stringToEnumType{
      {"raw", EType::raw},
      {"nocontinuum", EType::noContinuum},
      {"continuum", EType::continuumOnly},
      {"all", EType::all},
  };
  return stringToEnumType.at(component);
};

std::shared_ptr<CSolveResult> CTemplateFittingSolve::compute() {

  auto const &inputContext = Context.GetInputContext();
  auto const &resultStore = Context.GetResultStore();
  const CTemplateCatalog &tplCatalog = *(inputContext->GetTemplateCatalog());
  const std::shared_ptr<const CParameterStore> parameterStore =
      inputContext->GetParameterStore();

  PopulateParameters(parameterStore);
  UpdateParamsAndCatalogForFft(inputContext, tplCatalog);

  LogParameters();
  CheckTemplateCatalog(tplCatalog);

  Log.LogInfo(Formatter() << "Iterating over " << m_category.size()
                          << " tplCategories");
  Log.LogInfo(Formatter() << "Trying " << m_category.c_str() << " ("
                          << tplCatalog.GetTemplateCount(m_category)
                          << " templates)");

  // Calculates first pass for each template
  EType type = MakeTypeFromSpectrumComponentParam(m_spcComponent);
  std::vector<CMask> maskList;
  bool hasResult = false;
  for (auto tpl : tplCatalog.GetTemplateList(m_category)) {
    Solve(resultStore, tpl, m_overlapThreshold, maskList, type, m_interpolation,
          m_extinction, m_dustFit);
    hasResult = true;
  }
  // Question : a priori maskList ne sert à rien et pourrait être supprimé. À
  // confirmer
  if (!hasResult)
    THROWG(ErrorCode::INTERNAL_ERROR, "no result for any template");

  std::shared_ptr<const ExtremaResult> extremaResult;
  TZGridListParams zgridParams = {};
  std::shared_ptr<CTemplateFittingSolveResult> TemplateFittingSolveResult;

  // TODO split in two methods depending if use two pass or not
  // TODO clean
  if (useTwoPass()) {
    COperatorPdfz pdfz(m_opt_pdfcombination,
                       2 * m_secondPass_halfwindowsize, // peak separation
                       -1, m_opt_maxCandidate, m_zLogSampling, "FPE", true, 0);
    std::string scopeStr = MakeScopeStrFromTypeAndPass(type);
    std::shared_ptr<PdfCandidatesZResult> candidateResult =
        pdfz.Compute(BuildChisquareArray(resultStore, scopeStr));

    std::shared_ptr<COperatorTemplateFitting> templateFittingOperator =
        std::dynamic_pointer_cast<COperatorTemplateFitting>(
            m_continuumFittingOperator);
    templateFittingOperator->SetFirstPassCandidates(
        candidateResult->m_ranked_candidates);

    extremaResult = templateFittingOperator->BuildFirstPassExtremaResults(
        resultStore->GetScopedPerTemplateResult(scopeStr));

    std::string firstpassExtremaResultsStr = scopeStr;
    firstpassExtremaResultsStr.append("_extrema");
    resultStore->StoreScopedGlobalResult(firstpassExtremaResultsStr.c_str(),
                                         extremaResult);
    // Store first pass pdf results
    resultStore->StoreScopedGlobalResult("firstpass_pdf",
                                         pdfz.m_postmargZResult);
    resultStore->StoreScopedGlobalResult("firstpass_pdf_params",
                                         pdfz.m_postmargZResult);
    templateFittingOperator->SetRedshifts(m_redshifts);
    templateFittingOperator->init(m_secondPass_halfwindowsize, m_zLogSampling,
                                  m_redshifts, m_redshiftStep);
    templateFittingOperator->buildExtendedRedshifts();
    // TODO move this
    size_t nResults = m_results.size();
    // TODO size could be better initalized
    // TODO revoir la gestion dans le cas de plusieurs resultIdxs (i.e. cas all)
    // a priori fait nimp pour le moment
    for (size_t resultIdx = 0; resultIdx < nResults; ++resultIdx) {
      templateFittingOperator->updateRedshiftGridAndResults(
          m_results[resultIdx]);

      // TODO see if some of thi can be removed
      m_redshifts = m_results[resultIdx]->Redshifts;
      templateFittingOperator->m_redshifts = m_results[resultIdx]->Redshifts;

      // TODO check that the order is the same for
      // fpExtremaResult->m_ranked_candidates and
      // templateFittingOperator->m_firstpass_extremaResult
      // TODO is ugly - homogenize
      for (Int32 candidateIdx = 0;
           candidateIdx < extremaResult->m_ranked_candidates.size();
           ++candidateIdx) {
        // TODO here get zIdxsToCompute
        std::vector<Int32> zIdxsToCompute =
            templateFittingOperator->getzIdxsToCompute(
                m_redshifts,
                templateFittingOperator->m_extendedRedshifts[candidateIdx]);
        auto candidate = extremaResult->m_ranked_candidates[candidateIdx];
        std::shared_ptr<const CTemplate> tpl = tplCatalog.GetTemplateByName(
            {m_category}, candidate.second->fittedContinuum.name);
        Int32 igmIdx = candidate.second->fittedContinuum.meiksinIdx;
        Int32 ismIdx = Context.getFluxCorrectionCalzetti()->GetEbmvIndex(
            candidate.second->fittedContinuum.ebmvCoef);
        std::vector<CMask> maskList;
        // TODO check that the second pass redshift grid is used
        Solve(resultStore, tpl, m_overlapThreshold, maskList, type,
              m_interpolation, m_extinction, m_dustFit, ismIdx, igmIdx,
              candidate.first, zIdxsToCompute);
      }
    }
    zgridParams = templateFittingOperator->getSPZGridParams();
    scopeStr = MakeScopeStrFromTypeAndPass(type, false);
    COperatorPdfz pdfz2(m_opt_pdfcombination, m_redshiftSeparation, 0.0,
                        m_opt_maxCandidate, m_zLogSampling);
    std::shared_ptr<const ExtremaResult> fpExtremaResult =
        std::dynamic_pointer_cast<const ExtremaResult>(
            m_continuumFittingOperator->getFirstPassExtremaResults());
    auto chi2array = BuildChisquareArray(resultStore, scopeStr, fpExtremaResult,
                                         zgridParams);
    candidateResult = pdfz2.Compute(chi2array);

    extremaResult = buildExtremaResults(
        resultStore, scopeStr, candidateResult->m_ranked_candidates, tplCatalog,
        m_overlapThreshold, fpExtremaResult);
    resultStore->StoreScopedGlobalResult("pdf", pdfz2.m_postmargZResult);
    resultStore->StoreScopedGlobalResult("pdf_params", pdfz2.m_postmargZResult);
    resultStore->StoreScopedGlobalResult("candidatesresult", candidateResult);

    // for each extrema, get best model by reading from datastore and selecting
    // best fit
    /////////////////////////////////////////////////////////////////////////////////////
    if (m_fftProcessing) {
      tplCatalog.m_logsampling = false;
    }
    StoreExtremaResults(resultStore, extremaResult);
    TemplateFittingSolveResult = std::make_shared<CTemplateFittingSolveResult>(
        extremaResult->m_ranked_candidates[0].second, m_opt_pdfcombination,
        pdfz.m_postmargZResult->valMargEvidenceLog);

  } else {
    COperatorPdfz pdfz(m_opt_pdfcombination, m_redshiftSeparation, 0.0,
                       m_opt_maxCandidate, m_zLogSampling);
    std::string scopeStr = MakeScopeStrFromTypeAndPass(type);
    std::shared_ptr<PdfCandidatesZResult> candidateResult =
        pdfz.Compute(BuildChisquareArray(resultStore, scopeStr));

    extremaResult = buildExtremaResults(resultStore, scopeStr,
                                        candidateResult->m_ranked_candidates,
                                        tplCatalog, m_overlapThreshold);
    resultStore->StoreScopedGlobalResult("pdf", pdfz.m_postmargZResult);
    resultStore->StoreScopedGlobalResult("pdf_params", pdfz.m_postmargZResult);
    resultStore->StoreScopedGlobalResult("candidatesresult", candidateResult);

    // for each extrema, get best model by reading from datastore and selecting
    // best fit
    /////////////////////////////////////////////////////////////////////////////////////
    if (m_fftProcessing) {
      tplCatalog.m_logsampling = false;
    }
    StoreExtremaResults(resultStore, extremaResult);
    TemplateFittingSolveResult = std::make_shared<CTemplateFittingSolveResult>(
        extremaResult->m_ranked_candidates[0].second, m_opt_pdfcombination,
        pdfz.m_postmargZResult->valMargEvidenceLog);
  }
  // TODO add   m_opt_candidatesLogprobaCutThreshold =
  // parameterStore->GetScoped<Float64>("lineModel.extremaCutProbaThreshold");
  // option

  return TemplateFittingSolveResult;
}

void CTemplateFittingSolve::Solve(
    std::shared_ptr<COperatorResultStore> resultStore,
    const std::shared_ptr<const CTemplate> &tpl, Float64 m_overlapThreshold,
    std::vector<CMask> maskList, EType fittingType, std::string m_interpolation,
    bool m_extinction, bool m_dustFitting, Int32 FitEbmvIdx,
    Int32 FitMeiksinIdx, std::string parentId,
    std::vector<Int32> zIdxsToCompute) {

  // For saving initial spectra fitting types and template type
  std::vector<CSpectrum::EType> save_spcTypes;
  CSpectrum::EType save_tplType;
  for (auto spc : Context.getSpectra())
    save_spcTypes.push_back(spc->GetType());
  save_tplType = tpl->GetType();

  // If fitting type is all, loop on all spectrum fitting types
  // otherwise, just use the corresponding one
  std::vector<CSpectrum::EType> spectrumTypes;
  if (fittingType == EType::all) {
    spectrumTypes = {CSpectrum::EType::raw, CSpectrum::EType::noContinuum,
                     CSpectrum::EType::continuumOnly};
  } else
    spectrumTypes = {fittingTypeToSpectrumType.at(fittingType)};
  const bool isFirstPass = FitEbmvIdx == undefIdx;

  if (isFirstPass)
    m_results = std::vector<std::shared_ptr<CTemplateFittingResult>>(
        spectrumTypes.size(), std::make_shared<CTemplateFittingResult>(
                                  CTemplateFittingResult(m_redshifts.size())));
  for (size_t spectrumTypeIdx = 0; spectrumTypeIdx < spectrumTypes.size();
       ++spectrumTypeIdx) {
    // Adapts spectra and template fitting type if necessary
    CSpectrum::EType spectrumType = spectrumTypes[spectrumTypeIdx];
    if (fittingType == EType::all) {
      for (auto spc : Context.getSpectra())
        spc->SetType(spectrumType);
      tpl->SetType(spectrumType);
    }

    if (fittingType == EType::noContinuum)
      m_dustFitting = false;
    tpl->setRebinInterpMethod(m_interpolation);

    // TODO
    // ici on pourrait plutôt initialiser m_results à la bonne taille et update
    // le result situé à l'index visé
    // TODO find a more elegant way to test if is COPeratorTemplateFitting
    std::shared_ptr<CTemplateFittingResult> templateFittingResult;
    if (!m_fftProcessing) {
      templateFittingResult =
          std::dynamic_pointer_cast<COperatorTemplateFitting>(
              m_continuumFittingOperator)
              ->Compute(tpl, m_overlapThreshold, m_interpolation, m_extinction,
                        m_dustFitting, 0, CPriorHelper::TPriorZEList(),
                        FitEbmvIdx, FitMeiksinIdx, m_results[spectrumTypeIdx],
                        isFirstPass, zIdxsToCompute);
    } else {
      // Là est le point des FitEbmvIdx et FitMeiksinIdx
      // TODO impélmenter cette part fftprocesing ?
      templateFittingResult = std::dynamic_pointer_cast<CTemplateFittingResult>(
          m_continuumFittingOperator->Compute(
              tpl, m_overlapThreshold, m_interpolation, m_extinction,
              m_dustFitting, 0, CPriorHelper::TPriorZEList(), FitEbmvIdx,
              FitMeiksinIdx));
    }
    // TODO see how to deal with this error
    if (!templateFittingResult)
      THROWG(ErrorCode::INTERNAL_ERROR,
             "no results returned by templateFittingOperator");

    // Store results
    // TODO find parent ID if possible ?
    std::string scopeStr =
        MakeScopeStrFromTypeAndPass(fittingType, isFirstPass);
    if (parentId != "")
      scopeStr += "_" + parentId;
    // TODO add spectrumtype
    resultStore->StoreScopedPerTemplateResult(tpl, scopeStr.c_str(),
                                              templateFittingResult);
  }

  // Reset spectra fitting types and template type
  int i = 0;
  for (auto spc : Context.getSpectra())
    spc->SetType(save_spcTypes[i++]);
  tpl->SetType(save_tplType);
  // Question add timer ?
}

ChisquareArray CTemplateFittingSolve::BuildChisquareArray(
    std::shared_ptr<const COperatorResultStore> store,
    const std::string &scopeStr, std::shared_ptr<const ExtremaResult> fpResults,
    TZGridListParams zgridParams) const {
  ChisquareArray chisquarearray;

  Log.LogDetail("templatefittingsolver: building chisquare array");
  Log.LogDetail(Formatter()
                << "    templatefittingsolver: using results in scope: "
                << store->GetScopedName(scopeStr));

  // TODO pass parentZCand ?
  bool isSecondPass = bool(fpResults);
  TOperatorResultMap templatesResultsMap;
  if (isSecondPass) {
    chisquarearray.parentCandidates =
        m_continuumFittingOperator->getFirstPassCandidatesZByRank();
    ;
    TStringList fpResultsIds = fpResults->GetIDs();
    for (std::string id : fpResultsIds) {
      const TOperatorResultMap &tmpMap =
          store->GetScopedPerTemplateResult(scopeStr + "_" + id);
      templatesResultsMap.insert(tmpMap.begin(), tmpMap.end());
    }
  } else {
    templatesResultsMap = store->GetScopedPerTemplateResult(scopeStr);
  }
  // TODO check what is contained is the zgrid results
  // if ok set chi2 redshifts grid

  chisquarearray.cstLog = -1;
  // TODO adapt le zstep etc
  chisquarearray.zstep = m_coarseRedshiftStep;
  chisquarearray.zgridParams = zgridParams;

  for (TOperatorResultMap::const_iterator templateResultMap =
           templatesResultsMap.begin();
       templateResultMap != templatesResultsMap.end(); ++templateResultMap) {
    auto templateResult =
        std::dynamic_pointer_cast<const CTemplateFittingResult>(
            (*templateResultMap).second);
    Int32 nISM = -1;
    Int32 nIGM = -1;
    if (isSecondPass) {
      nISM = 1;
      nIGM = 1;
    } else {
      if (templateResult->ChiSquareIntermediate.size() > 0) {
        nISM = templateResult->ChiSquareIntermediate[0].size();
        if (templateResult->ChiSquareIntermediate[0].size() > 0) {
          nIGM = templateResult->ChiSquareIntermediate[0][0].size();
        }
      }
    }
    if (chisquarearray.cstLog == -1) {
      chisquarearray.cstLog = templateResult->CstLog;
      Log.LogInfo(Formatter() << "templatefittingsolver: using cstLog = "
                              << chisquarearray.cstLog);
    } else if (chisquarearray.cstLog != templateResult->CstLog) {
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "cstLog values do not correspond: val1="
                         << chisquarearray.cstLog
                         << " != val2=" << templateResult->CstLog);
    }
    // TODO move this and replace by a check ?
    if (chisquarearray.redshifts.size() == 0) {
      chisquarearray.redshifts = templateResult->Redshifts;
    }

    CZPrior zpriorhelper;
    for (Int32 kism = 0; kism < nISM; kism++) {
      for (Int32 kigm = 0; kigm < nIGM; kigm++) {
        chisquarearray.zpriors.push_back(zpriorhelper.GetConstantLogZPrior(
            templateResult->Redshifts.size()));

        // correct chi2 for ampl. marg. if necessary: todo add switch, currently
        // deactivated
        chisquarearray.chisquares.emplace_back(
            templateResult->ChiSquareIntermediate.size(), DBL_MAX);
        TFloat64List &logLikelihoodCorrected = chisquarearray.chisquares.back();
        for (Int32 kz = 0; kz < templateResult->Redshifts.size(); kz++) {
          logLikelihoodCorrected[kz] =
              templateResult->ChiSquareIntermediate
                  [kz][kism]
                  [kigm]; // + resultXXX->ScaleMargCorrectionTplratios[][]?;
        }
        Log.LogDetail(
            Formatter()
            << "    templatefittingsolver: Pdfz combine - prepared merit #"
            << chisquarearray.chisquares.size() - 1
            << " for model : " << ((*templateResultMap).first));
      }
    }
  }

  return chisquarearray;
}

// TODO move this in operator ?
std::shared_ptr<const ExtremaResult> CTemplateFittingSolve::buildExtremaResults(
    std::shared_ptr<const COperatorResultStore> store,
    const std::string &scopeStr, const TCandidateZbyRank &ranked_zCandidates,
    const CTemplateCatalog &tplCatalog, Float64 m_overlapThreshold,
    std::shared_ptr<const ExtremaResult> fpResults) {

  Log.LogDetail(
      "CTemplateFittingSolve::buildExtremaResults: building chisquare array");
  Log.LogDetail(
      Formatter()
      << "    templatefittingsolve:r using tplFitResultsMap in scope: "
      << store->GetScopedName(scopeStr));

  TOperatorResultMap tplFitResultsMap =
      store->GetScopedPerTemplateResult(scopeStr);
  bool isSecondPass = bool(fpResults);

  // TODO this is a duplicate with buildchi2 array. To put in common
  if (isSecondPass) {
    TStringList fpResultsIds = fpResults->GetIDs();
    for (std::string id : fpResultsIds) {
      const TOperatorResultMap &tmpMap =
          store->GetScopedPerTemplateResult(scopeStr + "_" + id);
      tplFitResultsMap.insert(tmpMap.begin(), tmpMap.end());
    }
  } else {
    tplFitResultsMap = store->GetScopedPerTemplateResult(scopeStr);
  }
  // end of duplicate

  Int32 extremumCount = ranked_zCandidates.size();

  auto firstResult = std::dynamic_pointer_cast<const CTemplateFittingResult>(
      (*tplFitResultsMap.begin()).second);
  const TFloat64List &redshifts = firstResult->Redshifts;

  // check all tplFitResultsMap  status
  for (auto &result : tplFitResultsMap) {
    auto tplFitResult =
        std::dynamic_pointer_cast<const CTemplateFittingResult>(result.second);
    if (tplFitResult->ChiSquare.size() != redshifts.size()) {
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "Size do not match among templatefitting "
                            "tplFitResultsMap, for tpl="
                         << result.first);
    }
    for (Int32 kz = 0; kz < tplFitResult->Redshifts.size(); kz++) {
      if (tplFitResult->Redshifts[kz] != redshifts[kz]) {
        THROWG(ErrorCode::INTERNAL_ERROR,
               Formatter() << "redshift vector is not the same for tpl="
                           << result.first);
      }
    }
  }

  std::shared_ptr<ExtremaResult> extremaResult =
      make_shared<ExtremaResult>(ranked_zCandidates);

  for (Int32 iExtremum = 0; iExtremum < extremumCount; iExtremum++) {
    Float64 z = ranked_zCandidates[iExtremum].second->Redshift;

    // find the corresponding Z
    auto itZ = std::find(redshifts.begin(), redshifts.end(), z);
    const Int32 zIndex = std::distance(redshifts.begin(), itZ);

    // find the min chisquare at corresponding redshift
    Float64 ChiSquare = DBL_MAX;
    std::string name = "";
    for (auto &r : tplFitResultsMap) {
      auto TplFitResult =
          std::dynamic_pointer_cast<const CTemplateFittingResult>(r.second);

      if (TplFitResult->ChiSquare[zIndex] < ChiSquare) {
        ChiSquare = TplFitResult->ChiSquare[zIndex];
        name = r.first;
      };
    }

    // Fill extrema Result
    auto TplFitResult = std::dynamic_pointer_cast<const CTemplateFittingResult>(
        tplFitResultsMap[name]);
    extremaResult->m_ranked_candidates[iExtremum]
        .second->fittedContinuum.merit = ChiSquare;
    extremaResult->m_ranked_candidates[iExtremum]
        .second->fittedContinuum.tplMeritPhot =
        TplFitResult->ChiSquarePhot[zIndex];
    extremaResult->m_ranked_candidates[iExtremum].second->fittedContinuum.name =
        name;
    extremaResult->m_ranked_candidates[iExtremum]
        .second->fittedContinuum.reducedChi2 =
        TplFitResult->ReducedChiSquare[zIndex];
    extremaResult->m_ranked_candidates[iExtremum]
        .second->fittedContinuum.meiksinIdx =
        TplFitResult->FitMeiksinIdx[zIndex];
    extremaResult->m_ranked_candidates[iExtremum]
        .second->fittedContinuum.ebmvCoef = TplFitResult->FitEbmvCoeff[zIndex];
    extremaResult->m_ranked_candidates[iExtremum]
        .second->fittedContinuum.tplAmplitude =
        TplFitResult->FitAmplitude[zIndex];
    extremaResult->m_ranked_candidates[iExtremum]
        .second->fittedContinuum.tplAmplitudeError =
        TplFitResult->FitAmplitudeError[zIndex];
    extremaResult->m_ranked_candidates[iExtremum]
        .second->fittedContinuum.tplDtM = TplFitResult->FitDtM[zIndex];
    extremaResult->m_ranked_candidates[iExtremum]
        .second->fittedContinuum.tplMtM = TplFitResult->FitMtM[zIndex];
    extremaResult->m_ranked_candidates[iExtremum].second->fittedContinuum.SNR =
        TplFitResult->SNR[zIndex];
    extremaResult->m_ranked_candidates[iExtremum]
        .second->fittedContinuum.tplLogPrior = TplFitResult->LogPrior[zIndex];

    // make sure tpl is non-rebinned
    bool currentSampling = tplCatalog.m_logsampling;
    tplCatalog.m_logsampling = false;
    std::shared_ptr<const CTemplate> tpl =
        tplCatalog.GetTemplateByName({m_category}, name);

    std::shared_ptr<CModelSpectrumResult> spcmodelPtr =
        std::make_shared<CModelSpectrumResult>();
    for (int spcIndex = 0; spcIndex < Context.getSpectra().size(); spcIndex++) {
      const std::string &obsId = Context.getSpectra()[spcIndex]->getObsID();

      TPhotVal values = m_continuumFittingOperator->ComputeSpectrumModel(
          tpl, z, TplFitResult->FitEbmvCoeff[zIndex],
          TplFitResult->FitMeiksinIdx[zIndex],
          TplFitResult->FitAmplitude[zIndex], m_overlapThreshold, spcIndex,
          spcmodelPtr);

      if (spcmodelPtr == nullptr)
        THROWG(ErrorCode::INTERNAL_ERROR, "Could not "
                                          "compute spectrum model");
      tplCatalog.m_logsampling = currentSampling;

      extremaResult->m_modelPhotValues[iExtremum] =
          std::make_shared<const CModelPhotValueResult>(values);
    }
    extremaResult->m_savedModelSpectrumResults[iExtremum] = spcmodelPtr;
  }

  return extremaResult;
}

void CTemplateFittingSolve::StoreExtremaResults(
    std::shared_ptr<COperatorResultStore> resultStore,
    std::shared_ptr<const ExtremaResult> &extremaResult) const {
  resultStore->StoreScopedGlobalResult("extrema_results", extremaResult);
  Log.LogInfo("CTemplateFittingSolve::StoreExtremaResults: Templatefitting, "
              "saving extrema results");

  return;
}

// TODO put un common with linemodelsolve ?
void CTemplateFittingSolve::initTwoPassZStepFactor() {
  // NB: To be used only if second pass is enabled
  m_twoPassZStepFactor =
      Context.GetInputContext()->GetParameterStore()->GetScoped<Int32>(
          "firstPass.largeGridStepRatio");
};
