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

#include "RedshiftLibrary/common/size.h"
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

  if (m_spectrumType == EType::noContinuum) {
    if (m_dustFit)
      THROWG(ErrorCode::BAD_PARAMETER_VALUE,
             "noContinuum option incompatible with ismFit");
    if (m_extinction)
      THROWG(ErrorCode::BAD_PARAMETER_VALUE,
             "noContinuum option incompatible with ismFit");
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

std::string CTemplateFittingSolve::getResultName() const {
  std::unordered_map<EType, std::string> spectrumTypeToStr = {
      {EType::raw, "templatefitting"},
      {EType::continuumOnly, "templatefitting_continuum"},
      {EType::noContinuum, "templatefitting_nocontinuum"},
  };

  std::string resultName = spectrumTypeToStr.at(m_spectrumType);
  if (twoPassIsActive() && m_isFirstPass)
    resultName += "_firstpass";

  return resultName;
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

  auto extremaResult = computeResults(pdfz);
  storeResults(pdfz, extremaResult);
  auto templateFittingSolveResult =
      std::make_shared<CTemplateFittingSolveResult>(
          extremaResult->getRankedCandidateCPtr(0), m_opt_pdfcombination,
          pdfz.m_postmargZResult->valMargEvidenceLog);
  return templateFittingSolveResult;
}

std::shared_ptr<CTemplateFittingSolveResult>
CTemplateFittingSolve::computeTwoPass() {
  // First pass
  computeFirstPass();

  COperatorPdfz pdfz(m_opt_pdfcombination,
                     2 * m_secondPass_halfwindowsize, // peak separation
                     m_opt_candidatesLogprobaCutThreshold, m_opt_maxCandidate,
                     m_zLogSampling, "FPE", true, 0);

  auto extremaResult = computeResults(pdfz);

  storeFirstPassResults(pdfz, extremaResult);

  // Second pass
  computeSecondPass(extremaResult);

  TZGridListParams zgridParams = m_templateFittingOperator->getSPZGridParams();
  COperatorPdfz pdfz2(m_opt_pdfcombination, 0.0,
                      m_opt_candidatesLogprobaCutThreshold, m_opt_extremacount,
                      m_zLogSampling, "SPE", false, 1);
  extremaResult = computeResults(pdfz2, zgridParams);
  storeResults(pdfz2, extremaResult);

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
    auto const tplFitResult = Solve(resultStore, tpl);
    hasResult = true;
    // Store results
    std::string resultName = getResultName();
    resultStore->StoreScopedPerTemplateResult(tpl, resultName, tplFitResult);
  }
  if (!hasResult)
    THROWG(ErrorCode::INTERNAL_ERROR, "no result for any template");
}

std::shared_ptr<const ExtremaResult>
CTemplateFittingSolve::computeResults(COperatorPdfz &pdfz,
                                      const TZGridListParams &zgridParams) {
  std::string resultName = getResultName();
  auto const &candidateResult =
      pdfz.Compute(BuildChisquareArray(resultName, zgridParams));

  auto const &extremaResult =
      buildExtremaResults(resultName, candidateResult->m_ranked_candidates);

  if (m_isFirstPass && twoPassIsActive())
    m_templateFittingOperator->SetFirstPassExtremaResults(extremaResult);

  return extremaResult;
}

void CTemplateFittingSolve::storeFirstPassResults(
    const COperatorPdfz &pdfz,
    std::shared_ptr<const ExtremaResult> const &extremaResult) {
  std::string const firstpassExtremaResultsStr =
      getResultName().append("_extrema");
  auto const &resultStore = Context.GetResultStore();
  resultStore->StoreScopedGlobalResult(firstpassExtremaResultsStr,
                                       extremaResult);
  resultStore->StoreScopedGlobalResult("firstpass_pdf", pdfz.m_postmargZResult);
  resultStore->StoreScopedGlobalResult("firstpass_pdf_params",
                                       pdfz.m_postmargZResult);
}

void CTemplateFittingSolve::storeResults(
    const COperatorPdfz &pdfz,
    std::shared_ptr<const ExtremaResult> const &extremaResult

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

  m_templateFittingOperator->setTwoPassParameters(
      m_secondPass_halfwindowsize, m_zLogSampling, m_redshiftStep,
      m_twoPassZStepFactor);
  m_templateFittingOperator->buildExtendedRedshifts();

  auto const insertionIndices = m_templateFittingOperator->updateRedshiftGrid();
  auto templatesResultsMap = getPerTemplateResultMapCopy(getResultName());
  for (auto const &[tplName, tplFitResult] : templatesResultsMap) {
    m_templateFittingOperator->updateVectors(tplFitResult, insertionIndices);
  }

  m_isFirstPass = false;

  for (Int32 candidateIdx = 0; candidateIdx < extremaResult->size();
       ++candidateIdx) {
    auto candidate = extremaResult->getRankedCandidateCPtr(candidateIdx);
    const std::string &candidateName = extremaResult->ID(candidateIdx);
    std::shared_ptr<const CTemplate> tpl = tplCatalog.GetTemplateByName(
        {m_category}, candidate->fittedContinuum.name);
    Int32 igmIdx = candidate->fittedContinuum.meiksinIdx;
    Int32 ismIdx = Context.getFluxCorrectionCalzetti()->GetEbmvIndex(
        candidate->fittedContinuum.ebmvCoef);
    auto const tplFitResult =
        templatesResultsMap[candidate->fittedContinuum.name];
    Solve(resultStore, tpl, ismIdx, igmIdx, candidateName, candidateIdx,
          tplFitResult);
  }

  // save all template results
  for (auto const [tplName, tplFitResult] : templatesResultsMap) {
    std::shared_ptr<const CTemplate> tpl =
        tplCatalog.GetTemplateByName({m_category}, tplName);
    resultStore->StoreScopedPerTemplateResult(tpl, getResultName(),
                                              tplFitResult);
  }
}

std::shared_ptr<CTemplateFittingResult> CTemplateFittingSolve::Solve(
    std::shared_ptr<COperatorResultStore> resultStore,
    const std::shared_ptr<const CTemplate> &tpl, Int32 FitEbmvIdx,
    Int32 FitMeiksinIdx, std::string parentId, Int32 candidateIdx,
    std::shared_ptr<CTemplateFittingResult> const &result) {

  // For saving initial spectra fitting types and template type
  std::vector<CSpectrum::EType> save_spcTypes;
  CSpectrum::EType save_tplType;
  for (auto spc : Context.getSpectra())
    save_spcTypes.push_back(spc->GetType());
  save_tplType = tpl->GetType();

  // If fitting type is all, loop on all spectrum fitting types
  // otherwise, just use the corresponding one
  CSpectrum::EType spectrumType = fittingTypeToSpectrumType.at(m_spectrumType);

  for (auto spc : Context.getSpectra())
    spc->SetType(spectrumType);
  tpl->SetType(spectrumType);

  if (m_spectrumType == EType::noContinuum)
    m_dustFit = false;
  tpl->setRebinInterpMethod(m_interpolation);

  TInt32Range zIdxRangeToCompute =
      candidateIdx == undefIdx
          ? TInt32Range(undefIdx, undefIdx)
          : m_templateFittingOperator->getzIdxRangeToCompute(candidateIdx);

  std::shared_ptr<CTemplateFittingResult> templateFittingResult =
      m_templateFittingOperator->Compute(
          tpl, m_overlapThreshold, m_interpolation, m_extinction, m_dustFit, 0,
          CPriorHelper::TPriorZEList(), FitEbmvIdx, FitMeiksinIdx,
          zIdxRangeToCompute, result);

  if (!templateFittingResult)
    THROWG(ErrorCode::INTERNAL_ERROR,
           "no results returned by templateFittingOperator");

  // Reset spectra fitting types and template type
  int i = 0;
  for (auto spc : Context.getSpectra())
    spc->SetType(save_spcTypes[i++]);
  tpl->SetType(save_tplType);

  return templateFittingResult;
}

ChisquareArray
CTemplateFittingSolve::BuildChisquareArray(const std::string &resultName,
                                           TZGridListParams zgridParams) const {
  ChisquareArray chisquarearray;

  Log.LogDetail("templatefittingsolver: building chisquare array");

  bool isSecondPass = !zgridParams.empty();
  if (isSecondPass)
    chisquarearray.parentCandidates =
        m_templateFittingOperator->getFirstPassCandidatesZByRank();

  // Question : why coarse step for second pass too ? This is valid for second
  // pass only ?
  chisquarearray.zstep = m_coarseRedshiftStep;
  chisquarearray.zgridParams = zgridParams;

  auto templatesResultsMap = getPerTemplateResultMap(resultName);
  chisquarearray.cstLog = templatesResultsMap.begin()->second->CstLog;
  Log.LogInfo(Formatter() << "templatefittingsolver: using cstLog = "
                          << chisquarearray.cstLog);
  chisquarearray.redshifts = templatesResultsMap.begin()->second->Redshifts;
  auto const zprior =
      CZPrior().GetConstantLogZPrior(chisquarearray.redshifts.size());
  for (auto const &[_, templateResult] : templatesResultsMap) {
    auto const &[nISM, nIGM] = templateResult->getIsmIgmSizes();

    if (chisquarearray.cstLog != templateResult->CstLog) {
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "cstLog values do not correspond: val1="
                         << chisquarearray.cstLog
                         << " != val2=" << templateResult->CstLog);
    }

    for (Int32 kism = 0; kism < nISM; kism++) {
      for (Int32 kigm = 0; kigm < nIGM; kigm++) {
        chisquarearray.zpriors.push_back(zprior);

        // Correct chi2 for ampl. marg
        chisquarearray.chisquares.emplace_back(
            templateResult->ChiSquareIntermediate.size(), DBL_MAX);
        TFloat64List &chisquare = chisquarearray.chisquares.back();
        for (Int32 kz = 0; kz < ssize(templateResult->Redshifts); ++kz) {
          chisquare[kz] = templateResult->ChiSquareIntermediate[kz][kism][kigm];
        }
      }
    }
  }

  return chisquarearray;
}

std::shared_ptr<ExtremaResult> CTemplateFittingSolve::buildExtremaResults(
    const std::string &resultName,
    const TCandidateZbyRank &ranked_zCandidates) {

  Log.LogDetail("CTemplateFittingSolve::buildExtremaResults");

  auto const tplFitResultsMap = getPerTemplateResultMap(resultName);

  Int32 extremumCount = ranked_zCandidates.size();

  auto firstResult = (*tplFitResultsMap.begin()).second;
  const TFloat64List &redshifts = firstResult->Redshifts;

  // check all tplFitResultsMap  status
  for (auto const &[tplName, tplFitResult] : tplFitResultsMap) {
    if (tplFitResult->ChiSquare.size() != redshifts.size()) {
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "Size do not match among templatefitting "
                            "tplFitResultsMap, for tpl="
                         << tplName);
    }
    for (Int32 kz = 0; kz < ssize(tplFitResult->Redshifts); kz++) {
      if (tplFitResult->Redshifts[kz] != redshifts[kz]) {
        THROWG(ErrorCode::INTERNAL_ERROR,
               Formatter() << "redshift vector is not the same for tpl="
                           << tplName);
      }
    }
  }

  std::shared_ptr<ExtremaResult> extremaResult =
      make_shared<ExtremaResult>(ranked_zCandidates);

  for (Int32 iExtremum = 0; iExtremum < extremumCount; iExtremum++) {
    auto candidate = extremaResult->getRankedCandidatePtr(iExtremum);
    Float64 z = candidate->Redshift;

    // find the corresponding Z
    auto const zIndex = CIndexing<Float64>::getIndex(redshifts, z);

    // find the min chisquare at corresponding redshift
    using TPairTplFitResult =
        std::pair<std::string, std::shared_ptr<const CTemplateFittingResult>>;
    auto const &[bestName, bestResult] = *std::min_element(
        tplFitResultsMap.cbegin(), tplFitResultsMap.cend(),
        [zIndex](TPairTplFitResult const &l, TPairTplFitResult const &r) {
          return l.second->ChiSquare[zIndex] < r.second->ChiSquare[zIndex];
        });

    // Fill extrema Result
    // only usefull attributes for 1st pass in two-pass mode, to build chisquare
    // Array
    candidate->fittedContinuum.name = bestName;
    candidate->fittedContinuum.meiksinIdx = bestResult->FitMeiksinIdx[zIndex];
    candidate->fittedContinuum.ebmvCoef = bestResult->FitEbmvCoeff[zIndex];
    if (m_isFirstPass && twoPassIsActive())
      continue;

    // Fill all the remaing attributes
    candidate->fittedContinuum.merit = bestResult->ChiSquare[zIndex];
    candidate->fittedContinuum.tplMeritPhot = bestResult->ChiSquarePhot[zIndex];
    candidate->fittedContinuum.reducedChi2 =
        bestResult->ReducedChiSquare[zIndex];
    candidate->fittedContinuum.pValue = bestResult->pValue[zIndex];
    candidate->fittedContinuum.tplAmplitude = bestResult->FitAmplitude[zIndex];
    candidate->fittedContinuum.tplAmplitude = bestResult->FitAmplitude[zIndex];
    candidate->fittedContinuum.tplAmplitudeError =
        bestResult->FitAmplitudeError[zIndex];
    candidate->fittedContinuum.tplDtM = bestResult->FitDtM[zIndex];
    candidate->fittedContinuum.tplMtM = bestResult->FitMtM[zIndex];
    candidate->fittedContinuum.SNR = bestResult->SNR[zIndex];
    candidate->fittedContinuum.tplLogPrior = bestResult->LogPrior[zIndex];

    // make sure tpl is non-rebinned
    const CTemplateCatalog &tplCatalog = *(Context.GetTemplateCatalog());
    bool currentSampling = tplCatalog.m_logsampling;
    tplCatalog.m_logsampling = false;
    std::shared_ptr<const CTemplate> tpl =
        tplCatalog.GetTemplateByName({m_category}, bestName);

    std::shared_ptr<CModelSpectrumResult> spcmodelPtr =
        std::make_shared<CModelSpectrumResult>();
    for (int spcIndex = 0; spcIndex < ssize(Context.getSpectra()); spcIndex++) {
      const std::string &obsId = Context.getSpectra()[spcIndex]->getObsID();

      TPhotVal values = m_templateFittingOperator->ComputeSpectrumModel(
          tpl, z, bestResult->FitEbmvCoeff[zIndex],
          bestResult->FitMeiksinIdx[zIndex], bestResult->FitAmplitude[zIndex],
          m_overlapThreshold, spcIndex, spcmodelPtr);

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

TConstTemplateFittingResultMap CTemplateFittingSolve::getPerTemplateResultMap(
    const std::string &resultName) const {
  auto const &resultStore = Context.GetResultStore();
  auto const resultsMap = resultStore->GetScopedPerTemplateResult(resultName);

  TConstTemplateFittingResultMap tplFitResultMap;
  for (auto const &[name, result] : resultsMap)
    tplFitResultMap[name] =
        std::dynamic_pointer_cast<const CTemplateFittingResult>(result);

  return tplFitResultMap;
}

TTemplateFittingResultMap CTemplateFittingSolve::getPerTemplateResultMapCopy(
    const std::string &resultName) const {
  auto const &resultStore = Context.GetResultStore();
  auto const resultsMap = resultStore->GetScopedPerTemplateResult(resultName);

  TTemplateFittingResultMap tplFitResultMap;
  for (auto const &[name, result] : resultsMap)
    tplFitResultMap[name] = std::make_shared<CTemplateFittingResult>(
        *std::dynamic_pointer_cast<const CTemplateFittingResult>(result));

  return tplFitResultMap;
}
