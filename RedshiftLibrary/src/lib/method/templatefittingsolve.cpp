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

#include "RedshiftLibrary/operator/modelspectrumresult.h"

#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/operator/modelphotvalueresult.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/operator/templatefitting.h"
#include "RedshiftLibrary/operator/templatefittinglog.h"
#include "RedshiftLibrary/operator/templatefittingresult.h"
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
  std::string spcComponent =
      parameterStore->GetScoped<std::string>("spectrum.component");
  m_spectrumType = getFitTypeFromParam(spcComponent);
  m_interpolation = parameterStore->GetScoped<std::string>("interpolation");
  m_extinction = parameterStore->GetScoped<bool>("igmFit");
  m_dustFit = parameterStore->GetScoped<bool>("ismFit");

  if (m_spcComponent == "noContinuum")
    {
	if (m_dustFit)
	  THROWG(ErrorCode::BAD_PARAMETER_VALUE, "noContinuum option incompatible with ismFit");
	if (m_extinction)
	  THROWG(ErrorCode::BAD_PARAMETER_VALUE, "noContinuum option incompatible with ismFit");
    }
		     
  m_fftProcessing = parameterStore->GetScoped<bool>("fftProcessing");

  if (parameterStore->HasScoped<bool>("enablePhotometry")) {
    m_usePhotometry = parameterStore->GetScoped<bool>("enablePhotometry");
    if (m_usePhotometry)
      m_photometryWeight =
          parameterStore->GetScoped<Float64>("photometry.weight");
  }
  m_opt_extremacount = parameterStore->GetScoped<int>("extremaCount");
  m_opt_pdfcombination =
      parameterStore->GetScoped<std::string>("pdfCombination");
  m_opt_singlePass = parameterStore->GetScoped<bool>("singlePass");
  if (!isSinglePass()) {
    m_opt_maxCandidate =
        parameterStore->GetScoped<int>("firstPass.extremaCount");
    m_secondPassContinuumFit = str2ContinuumFit.at(
        parameterStore->GetScoped<std::string>("secondPass.continuumFit"));

    m_secondPass_halfwindowsize =
        parameterStore->GetScoped<Float64>("secondPass.halfWindowSize");
  }
}

void CTemplateFittingSolve::InitFittingOperator() {
  const CTemplateCatalog &tplCatalog = *(Context.GetTemplateCatalog());
  if (m_fftProcessing && m_usePhotometry)
    THROWG(ErrorCode::FFT_WITH_PHOTOMETRY_NOTIMPLEMENTED,
           "fftProcessing not "
           "implemented with photometry enabled");

  if (m_fftProcessing) {
    m_templateFittingOperator =
        std::make_shared<COperatorTemplateFittingLog>(m_redshifts);
    tplCatalog.m_logsampling = true;
  } else {
    if (m_usePhotometry) {
      m_templateFittingOperator =
          std::make_shared<COperatorTemplateFittingPhot>(
              Context.GetPhotBandCatalog(), m_redshifts, m_photometryWeight);
    } else {
      m_templateFittingOperator =
          std::make_shared<COperatorTemplateFitting>(m_redshifts);
    }
    tplCatalog.m_logsampling = false;
  }
}

void CTemplateFittingSolve::LogParameters() {
  Log.LogInfo("Method parameters:");
  Log.LogInfo(Formatter() << "    -m_overlapThreshold: " << m_overlapThreshold);
  Log.LogInfo(Formatter() << "    -component: "
                          << static_cast<Int32>(m_spectrumType));
  Log.LogInfo(Formatter() << "    -interp: " << m_interpolation);
  Log.LogInfo(Formatter() << "    -IGM extinction: "
                          << (m_extinction ? "true" : "false"));
  Log.LogInfo(Formatter() << "    -ISM dust-fit: "
                          << (m_dustFit ? "true" : "false"));
  Log.LogInfo(Formatter() << "    -pdfcombination: " << m_opt_pdfcombination);
  Log.LogInfo("");
}

void CTemplateFittingSolve::CheckTemplateCatalog() {
  const auto &tplCatalog = *(Context.GetTemplateCatalog());
  if (tplCatalog.GetTemplateCount(m_category) == 0) {
    THROWG(ErrorCode::BAD_TEMPLATECATALOG,
           Formatter() << "Empty template catalog for category "
                       << m_category[0]);
  }
}

std::string CTemplateFittingSolve::getScopeStr() const {
  std::unordered_map<EType, std::string> spectrumTypeToStr = {
      {EType::raw, "templatefitting"},
      {EType::continuumOnly, "templatefitting_continuum"},
      {EType::noContinuum, "templatefitting_nocontinuum"},
  };

  std::string scopeStr = spectrumTypeToStr.at(m_spectrumType);
  if (twoPassIsActive() && m_isFirstPass)
    scopeStr += "_firstpass";

  return scopeStr;
}

CTemplateFittingSolve::EType
CTemplateFittingSolve::getFitTypeFromParam(const std::string &component) {
  std::unordered_map<std::string, EType> stringToEnumType{
      {"raw", EType::raw},
      {"noContinuum", EType::noContinuum},
      {"continuum", EType::continuumOnly},
  };
  return stringToEnumType.at(component);
};

std::shared_ptr<CSolveResult> CTemplateFittingSolve::compute() {
  const std::shared_ptr<const CParameterStore> parameterStore =
      Context.GetParameterStore();

  PopulateParameters(parameterStore);
  InitFittingOperator();
  LogParameters();
  CheckTemplateCatalog();
  const auto tplCatalog = *(Context.GetTemplateCatalog());
  Log.LogInfo(Formatter() << "Iterating over " << m_category.size()
                          << " tplCategories");
  Log.LogInfo(Formatter() << "Trying " << m_category.c_str() << " ("
                          << tplCatalog.GetTemplateCount(m_category)
                          << " templates)");

  std::shared_ptr<CTemplateFittingSolveResult> TemplateFittingSolveResult;
  // NB: this m_castedTemplateFittingOperator should disapper when 2 pass will
  // be implemented with FFT too
  if (!m_fftProcessing)
    m_castedTemplateFittingOperator =
        std::dynamic_pointer_cast<COperatorTemplateFitting>(
            m_templateFittingOperator);
  if (!isSinglePass()) {
    TemplateFittingSolveResult = computeTwoPass();
  } else {
    TemplateFittingSolveResult = computeSinglePass();
  }

  if (m_fftProcessing)
    tplCatalog.m_logsampling = false;

  return TemplateFittingSolveResult;
}

std::shared_ptr<CTemplateFittingSolveResult>
CTemplateFittingSolve::computeSinglePass() {
  computeFirstPass();
  COperatorPdfz pdfz(m_opt_pdfcombination, m_redshiftSeparation,
                     m_opt_candidatesLogprobaCutThreshold, m_opt_extremacount,
                     m_zLogSampling);
  std::string fpScopeStr = getScopeStr();
  std::shared_ptr<PdfCandidatesZResult> candidateResult =
      pdfz.Compute(BuildChisquareArray(fpScopeStr));
  std::shared_ptr<const ExtremaResult> extremaResult = buildExtremaResults(
      fpScopeStr, candidateResult->m_ranked_candidates, m_overlapThreshold);
  storeSinglePassResults(pdfz, extremaResult);
  auto templateFittingSolveResult =
      std::make_shared<CTemplateFittingSolveResult>(
          extremaResult->getRankedCandidateCPtr(0), m_opt_pdfcombination,
          pdfz.m_postmargZResult->valMargEvidenceLog);
  return templateFittingSolveResult;
}

void CTemplateFittingSolve::storeSinglePassResults(
    const COperatorPdfz &pdfz,
    std::shared_ptr<const ExtremaResult> extremaResult) {
  auto const &resultStore = Context.GetResultStore();
  resultStore->StoreScopedGlobalResult("pdf", pdfz.m_postmargZResult);
  resultStore->StoreScopedGlobalResult("pdf_params", pdfz.m_postmargZResult);
  resultStore->StoreScopedGlobalResult("extrema_results", extremaResult);
  Log.LogInfo("CTemplateFittingSolve::StoreExtremaResults: Templatefitting, "
              "saving extrema results");
}

std::shared_ptr<CTemplateFittingSolveResult>
CTemplateFittingSolve::computeTwoPass() {
  // First pass
  computeFirstPass();

  COperatorPdfz pdfz(m_opt_pdfcombination,
                     2 * m_secondPass_halfwindowsize, // peak separation
                     m_opt_candidatesLogprobaCutThreshold, m_opt_maxCandidate,
                     m_zLogSampling, "FPE", true, 0);

  auto results = computeFirstPassResults(pdfz);
  auto extremaResult = results.second;

  storeFirstPassResults(pdfz, extremaResult);

  // Second pass
  m_isFirstPass = false;
  computeSecondPass(extremaResult);

  TZGridListParams zgridParams =
      m_castedTemplateFittingOperator->getSPZGridParams();
  COperatorPdfz pdfz2(m_opt_pdfcombination, 0.0,
                      m_opt_candidatesLogprobaCutThreshold, m_opt_extremacount,
                      m_zLogSampling, "SPE", false, 1);
  results = computeSecondPassResults(pdfz2, zgridParams);
  extremaResult = results.second;
  storeSecondPassResults(pdfz2, extremaResult);

  auto templateFittingSolveResult =
      std::make_shared<CTemplateFittingSolveResult>(
          extremaResult->getRankedCandidateCPtr(0), m_opt_pdfcombination,
          pdfz.m_postmargZResult->valMargEvidenceLog);
  return templateFittingSolveResult;
}

void CTemplateFittingSolve::computeFirstPass() {
  bool hasResult = false;
  auto const &resultStore = Context.GetResultStore();
  const CTemplateCatalog &tplCatalog = *(Context.GetTemplateCatalog());

  for (auto tpl : tplCatalog.GetTemplateList(m_category)) {
    Solve(resultStore, tpl, m_overlapThreshold, m_interpolation, m_extinction,
          m_dustFit);
    hasResult = true;
  }
  if (!hasResult)
    THROWG(ErrorCode::INTERNAL_ERROR, "no result for any template");
}

std::pair<std::shared_ptr<PdfCandidatesZResult>,
          std::shared_ptr<const ExtremaResult>>
CTemplateFittingSolve::computeFirstPassResults(COperatorPdfz &pdfz) {
  std::string scopeStr = getScopeStr();
  auto chi2array = BuildChisquareArray(scopeStr);
  std::shared_ptr<PdfCandidatesZResult> candidateResult =
      pdfz.Compute(chi2array);

  m_castedTemplateFittingOperator->SetFirstPassCandidates(
      candidateResult->m_ranked_candidates);

  auto const &resultStore = Context.GetResultStore();
  std::shared_ptr<const ExtremaResult> extremaResult =
      m_castedTemplateFittingOperator->BuildFirstPassExtremaResults(
          resultStore->GetScopedPerTemplateResult(scopeStr));
  return std::pair<std::shared_ptr<PdfCandidatesZResult>,
                   std::shared_ptr<const ExtremaResult>>(candidateResult,
                                                         extremaResult);
}

std::pair<std::shared_ptr<PdfCandidatesZResult>,
          std::shared_ptr<const ExtremaResult>>
CTemplateFittingSolve::computeSecondPassResults(
    COperatorPdfz &pdfz, const TZGridListParams &zgridParams) {
  std::string scopeStr = getScopeStr();
  std::shared_ptr<const ExtremaResult> fpExtremaResult =
      m_castedTemplateFittingOperator->getFirstPassExtremaResults();
  std::shared_ptr<PdfCandidatesZResult> candidateResult =
      pdfz.Compute(BuildChisquareArray(scopeStr, fpExtremaResult, zgridParams));
  m_castedTemplateFittingOperator->SetFirstPassCandidates(
      candidateResult->m_ranked_candidates);

  std::shared_ptr<const ExtremaResult> extremaResult =
      buildExtremaResults(scopeStr, candidateResult->m_ranked_candidates,
                          m_overlapThreshold, fpExtremaResult);
  return std::pair<std::shared_ptr<PdfCandidatesZResult>,
                   std::shared_ptr<const ExtremaResult>>(candidateResult,
                                                         extremaResult);
}

void CTemplateFittingSolve::storeFirstPassResults(
    const COperatorPdfz &pdfz,
    std::shared_ptr<const ExtremaResult> extremaResult) {
  std::string firstpassExtremaResultsStr = getScopeStr();
  auto const &resultStore = Context.GetResultStore();
  firstpassExtremaResultsStr.append("_extrema");
  resultStore->StoreScopedGlobalResult(firstpassExtremaResultsStr.c_str(),
                                       extremaResult);
  resultStore->StoreScopedGlobalResult("firstpass_pdf", pdfz.m_postmargZResult);
  resultStore->StoreScopedGlobalResult("firstpass_pdf_params",
                                       pdfz.m_postmargZResult);
}

void CTemplateFittingSolve::storeSecondPassResults(
    const COperatorPdfz &pdfz,
    std::shared_ptr<const ExtremaResult> extremaResult

) {
  auto const &resultStore = Context.GetResultStore();
  resultStore->StoreScopedGlobalResult("pdf", pdfz.m_postmargZResult);
  resultStore->StoreScopedGlobalResult("pdf_params", pdfz.m_postmargZResult);
  resultStore->StoreScopedGlobalResult("extrema_results", extremaResult);
  Log.LogInfo("CTemplateFittingSolve::StoreExtremaResults: Templatefitting, "
              "saving extrema results");
}

void CTemplateFittingSolve::computeSecondPass(
    std::shared_ptr<const ExtremaResult> extremaResult) {
  const CTemplateCatalog &tplCatalog = *(Context.GetTemplateCatalog());
  auto const &resultStore = Context.GetResultStore();

  m_castedTemplateFittingOperator->SetRedshifts(m_redshifts);
  m_castedTemplateFittingOperator->COperatorTwoPass::setClassVariables(
      m_secondPass_halfwindowsize, m_zLogSampling, m_redshifts, m_redshiftStep);
  m_castedTemplateFittingOperator->buildExtendedRedshifts();
  m_castedTemplateFittingOperator->updateRedshiftGridAndResults(m_result);

  m_redshifts = m_result->Redshifts;
  m_castedTemplateFittingOperator->SetRedshifts(m_result->Redshifts);
  for (Int32 candidateIdx = 0; candidateIdx < extremaResult->size();
       ++candidateIdx) {
    std::vector<Int32> zIdxsToCompute =
        m_castedTemplateFittingOperator->getzIdxsToCompute(
            m_redshifts, m_castedTemplateFittingOperator
                             ->getExtendedRedshifts()[candidateIdx]);
    auto candidate = extremaResult->getRankedCandidateCPtr(candidateIdx);
    const std::string &candidateName = extremaResult->ID(candidateIdx);
    std::shared_ptr<const CTemplate> tpl = tplCatalog.GetTemplateByName(
        {m_category}, candidate->fittedContinuum.name);
    Int32 igmIdx = candidate->fittedContinuum.meiksinIdx;
    Int32 ismIdx = Context.getFluxCorrectionCalzetti()->GetEbmvIndex(
        candidate->fittedContinuum.ebmvCoef);
    Solve(resultStore, tpl, m_overlapThreshold, m_interpolation, m_extinction,
          m_dustFit, ismIdx, igmIdx, candidateName, zIdxsToCompute);
  }
}

void CTemplateFittingSolve::Solve(
    std::shared_ptr<COperatorResultStore> resultStore,
    const std::shared_ptr<const CTemplate> &tpl, Float64 m_overlapThreshold,
    std::string m_interpolation, bool m_extinction, bool m_dustFitting,
    Int32 FitEbmvIdx, Int32 FitMeiksinIdx, std::string parentId,
    std::vector<Int32> zIdxsToCompute) {

  // For saving initial spectra fitting types and template type
  std::vector<CSpectrum::EType> save_spcTypes;
  CSpectrum::EType save_tplType;
  for (auto spc : Context.getSpectra())
    save_spcTypes.push_back(spc->GetType());
  save_tplType = tpl->GetType();

  // If fitting type is all, loop on all spectrum fitting types
  // otherwise, just use the corresponding one
  CSpectrum::EType spectrumType = fittingTypeToSpectrumType.at(m_spectrumType);

  if (m_isFirstPass)
    m_result = std::make_shared<CTemplateFittingResult>(
        CTemplateFittingResult(m_redshifts.size()));

  for (auto spc : Context.getSpectra())
    spc->SetType(spectrumType);
  tpl->SetType(spectrumType);

  if (m_spectrumType == EType::noContinuum)
    m_dustFitting = false;
  tpl->setRebinInterpMethod(m_interpolation);
  std::shared_ptr<CTemplateFittingResult> templateFittingResult;
  if (!m_fftProcessing) {
    templateFittingResult = m_castedTemplateFittingOperator->Compute(
        tpl, m_overlapThreshold, m_interpolation, m_extinction, m_dustFitting,
        0, CPriorHelper::TPriorZEList(), FitEbmvIdx, FitMeiksinIdx, m_result,
        m_isFirstPass, zIdxsToCompute);
  } else {
    templateFittingResult = m_templateFittingOperator->Compute(
        tpl, m_overlapThreshold, m_interpolation, m_extinction, m_dustFitting,
        0, CPriorHelper::TPriorZEList(), FitEbmvIdx, FitMeiksinIdx);
  }
  if (!templateFittingResult)
    THROWG(ErrorCode::INTERNAL_ERROR,
           "no results returned by templateFittingOperator");

  // Store results
  std::string scopeStr = getScopeStr();
  if (parentId != "")
    scopeStr += "_" + parentId;
  resultStore->StoreScopedPerTemplateResult(tpl, scopeStr.c_str(),
                                            templateFittingResult);

  // Reset spectra fitting types and template type
  int i = 0;
  for (auto spc : Context.getSpectra())
    spc->SetType(save_spcTypes[i++]);
  tpl->SetType(save_tplType);
}

ChisquareArray CTemplateFittingSolve::BuildChisquareArray(
    const std::string &scopeStr, std::shared_ptr<const ExtremaResult> fpResults,
    TZGridListParams zgridParams) const {
  ChisquareArray chisquarearray;

  Log.LogDetail("templatefittingsolver: building chisquare array");

  TOperatorResultMap templatesResultsMap =
      createPerTemplateResultMap(scopeStr, fpResults);
  bool isSecondPass = bool(fpResults);
  if (isSecondPass)
    chisquarearray.parentCandidates =
        m_templateFittingOperator->getFirstPassCandidatesZByRank();

  // Should set cstLog ? What does it correspond to ?
  chisquarearray.cstLog = -1;
  // Question : why coarse step for second pass too ? This is valid for second
  // pass only ?
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

    // NB for the moment this if else is useless as cstLog set to -1 at the
    // beginning of this method
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
    chisquarearray.redshifts = templateResult->Redshifts;

    CZPrior zpriorhelper;
    for (Int32 kism = 0; kism < nISM; kism++) {
      for (Int32 kigm = 0; kigm < nIGM; kigm++) {
        chisquarearray.zpriors.push_back(zpriorhelper.GetConstantLogZPrior(
            templateResult->Redshifts.size()));

        // Correct chi2 for ampl. marg
        chisquarearray.chisquares.emplace_back(
            templateResult->ChiSquareIntermediate.size(), DBL_MAX);
        TFloat64List &logLikelihoodCorrected = chisquarearray.chisquares.back();
        for (Int32 kz = 0; kz < templateResult->Redshifts.size(); kz++) {
          logLikelihoodCorrected[kz] =
              templateResult->ChiSquareIntermediate[kz][kism][kigm];
        }
      }
    }
  }

  return chisquarearray;
}

std::shared_ptr<const ExtremaResult> CTemplateFittingSolve::buildExtremaResults(
    const std::string &scopeStr, const TCandidateZbyRank &ranked_zCandidates,
    Float64 m_overlapThreshold,
    std::shared_ptr<const ExtremaResult> fpResults) {

  Log.LogDetail(
      "CTemplateFittingSolve::buildExtremaResults: building chisquare array");

  TOperatorResultMap tplFitResultsMap =
      createPerTemplateResultMap(scopeStr, fpResults);

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
    auto candidate = extremaResult->getRankedCandidatePtr(iExtremum);
    candidate->fittedContinuum.merit = ChiSquare;
    candidate->fittedContinuum.tplMeritPhot =
        TplFitResult->ChiSquarePhot[zIndex];
    candidate->fittedContinuum.name = name;
    candidate->fittedContinuum.reducedChi2 =
        TplFitResult->ReducedChiSquare[zIndex];
    candidate->fittedContinuum.meiksinIdx = TplFitResult->FitMeiksinIdx[zIndex];
    candidate->fittedContinuum.ebmvCoef = TplFitResult->FitEbmvCoeff[zIndex];
    candidate->fittedContinuum.tplAmplitude =
        TplFitResult->FitAmplitude[zIndex];
    candidate->fittedContinuum.tplAmplitudeError =
        TplFitResult->FitAmplitudeError[zIndex];
    candidate->fittedContinuum.tplDtM = TplFitResult->FitDtM[zIndex];
    candidate->fittedContinuum.tplMtM = TplFitResult->FitMtM[zIndex];
    candidate->fittedContinuum.SNR = TplFitResult->SNR[zIndex];
    candidate->fittedContinuum.tplLogPrior = TplFitResult->LogPrior[zIndex];

    // make sure tpl is non-rebinned
    const CTemplateCatalog &tplCatalog = *(Context.GetTemplateCatalog());
    bool currentSampling = tplCatalog.m_logsampling;
    tplCatalog.m_logsampling = false;
    std::shared_ptr<const CTemplate> tpl =
        tplCatalog.GetTemplateByName({m_category}, name);

    std::shared_ptr<CModelSpectrumResult> spcmodelPtr =
        std::make_shared<CModelSpectrumResult>();
    for (int spcIndex = 0; spcIndex < Context.getSpectra().size(); spcIndex++) {
      const std::string &obsId = Context.getSpectra()[spcIndex]->getObsID();

      TPhotVal values = m_templateFittingOperator->ComputeSpectrumModel(
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

void CTemplateFittingSolve::initSkipSecondPass() {
  m_opt_skipsecondpass = false;
  m_opt_singlePass =
      Context.GetInputContext()->GetParameterStore()->GetScoped<bool>(
          "singlePass");
};

void CTemplateFittingSolve::initTwoPassZStepFactor() {
  // NB: To be used only if second pass is enabled
  m_twoPassZStepFactor =
      Context.GetInputContext()->GetParameterStore()->GetScoped<Int32>(
          "firstPass.largeGridStepRatio");
};

TOperatorResultMap CTemplateFittingSolve::createPerTemplateResultMap(
    const std::string &scopeStr,
    std::shared_ptr<const ExtremaResult> fpResults) const {
  auto const &resultStore = Context.GetResultStore();
  TOperatorResultMap tplFitResultsMap =
      resultStore->GetScopedPerTemplateResult(scopeStr);
  bool isSecondPass = bool(fpResults);
  if (isSecondPass) {
    TStringList fpResultsIds = fpResults->GetIDs();
    for (std::string id : fpResultsIds) {
      const TOperatorResultMap &tmpMap =
          resultStore->GetScopedPerTemplateResult(scopeStr + "_" + id);
      tplFitResultsMap.insert(tmpMap.begin(), tmpMap.end());
    }
  } else {
    tplFitResultsMap = resultStore->GetScopedPerTemplateResult(scopeStr);
  }
  return tplFitResultsMap;
}
