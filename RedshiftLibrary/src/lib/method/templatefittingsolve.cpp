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
    : CObjectSolve("templateFittingSolve") {}

std::shared_ptr<CSolveResult> CTemplateFittingSolve::compute() {

  auto const &inputContext = Context.GetInputContext();
  auto const &resultStore = Context.GetResultStore();
  const CTemplateCatalog &tplCatalog = *(inputContext->GetTemplateCatalog());

  m_redshiftSeparation = inputContext->GetParameterStore()->Get<Float64>(
      "extremaRedshiftSeparation"); // todo: deci

  bool storeResult = false;
  Float64 overlapThreshold =
      inputContext->GetParameterStore()->GetScoped<Float64>("overlapThreshold");
  std::string opt_spcComponent =
      inputContext->GetParameterStore()->GetScoped<std::string>(
          "spectrum.component");
  std::string opt_interp =
      inputContext->GetParameterStore()->GetScoped<std::string>(
          "interpolation");
  const bool opt_extinction =
      inputContext->GetParameterStore()->GetScoped<bool>("igmFit");
  bool opt_dustFit =
      inputContext->GetParameterStore()->GetScoped<bool>("ismFit");

  bool fft_processing =
      inputContext->GetParameterStore()->GetScoped<bool>("fftProcessing");

  bool use_photometry = false;
  Float64 photometry_weight = NAN;
  if (inputContext->GetParameterStore()->HasScoped<bool>("enablePhotometry")) {
    use_photometry =
        inputContext->GetParameterStore()->GetScoped<bool>("enablePhotometry");
    if (use_photometry)
      photometry_weight = inputContext->GetParameterStore()->GetScoped<Float64>(
          "photometry.weight");
  }
  if (fft_processing && use_photometry)
    THROWG(ErrorCode::FFT_WITH_PHOTOMETRY_NOTIMPLEMENTED,
           "fftProcessing not "
           "implemented with photometry enabled");

  if (fft_processing) {
    m_continuumFittingOperator =
        std::make_shared<COperatorTemplateFittingLog>(m_redshifts);
    tplCatalog.m_logsampling = true;
  } else {
    if (use_photometry) {
      m_continuumFittingOperator =
          std::make_shared<COperatorTemplateFittingPhot>(
              inputContext->GetPhotBandCatalog(), photometry_weight,
              m_redshifts);
    } else {
      m_continuumFittingOperator =
          std::make_shared<COperatorTemplateFitting>(m_redshifts);
    }
    tplCatalog.m_logsampling = false;
  }

  // prepare the unused masks
  std::vector<CMask> maskList;
  // define the redshift search grid
  //         Log.LogInfo("Stellar fitting redshift range = [%.5f, %.5f],
  //         step=%.6f", starRedshiftRange.GetBegin(),
  //         starRedshiftRange.GetEnd(), starRedshiftStep);

  std::string scopeStr = "templatefitting";

  EType _type = nType_raw;
  if (opt_spcComponent == "raw") {
    _type = nType_raw;
  } else if (opt_spcComponent == "nocontinuum") {
    _type = nType_noContinuum;
    scopeStr = "templatefitting_nocontinuum";
  } else if (opt_spcComponent == "continuum") {
    _type = nType_continuumOnly;
    scopeStr = "templatefitting_continuum";
  } else if (opt_spcComponent == "all") {
    _type = nType_all;
  } else {
    THROWG(ErrorCode::INTERNAL_ERROR, "Unknown spectrum component");
  }

  m_opt_maxCandidate =
      inputContext->GetParameterStore()->GetScoped<int>("extremaCount");
  m_opt_pdfcombination =
      inputContext->GetParameterStore()->GetScoped<std::string>(
          "pdfCombination");

  Log.LogInfo("Method parameters:");
  Log.LogInfo(Formatter() << "    -overlapThreshold: " << overlapThreshold);
  Log.LogInfo(Formatter() << "    -component: " << opt_spcComponent);
  Log.LogInfo(Formatter() << "    -interp: " << opt_interp);
  Log.LogInfo(Formatter() << "    -IGM extinction: "
                          << (opt_extinction ? "true" : "false"));
  Log.LogInfo(Formatter() << "    -ISM dust-fit: "
                          << (opt_dustFit ? "true" : "false"));
  Log.LogInfo(Formatter() << "    -pdfcombination: " << m_opt_pdfcombination);
  Log.LogInfo("");

  if (tplCatalog.GetTemplateCount(m_category) == 0) {
    THROWG(ErrorCode::BAD_TEMPLATECATALOG,
           Formatter() << "Empty template catalog for category "
                       << m_category[0]);
  }
  Log.LogInfo(Formatter() << "Iterating over " << m_category.size()
                          << " tplCategories");

  Log.LogInfo(Formatter() << "Trying " << m_category.c_str() << " ("
                          << tplCatalog.GetTemplateCount(m_category)
                          << " templates)");
  for (Int32 j = 0; j < tplCatalog.GetTemplateCount(m_category); j++) {
    std::shared_ptr<const CTemplate> tpl =
        tplCatalog.GetTemplate(m_category, j);

    Solve(resultStore, tpl, overlapThreshold, maskList, _type, opt_interp,
          opt_extinction, opt_dustFit);

    storeResult = true;
  }

  if (!storeResult)
    THROWG(ErrorCode::INTERNAL_ERROR, "no result for any template");

  COperatorPdfz pdfz(m_opt_pdfcombination, m_redshiftSeparation, 0.0,
                     m_opt_maxCandidate, m_redshiftSampling == "log");

  std::shared_ptr<PdfCandidatesZResult> candidateResult =
      pdfz.Compute(BuildChisquareArray(resultStore, scopeStr));

  resultStore->StoreScopedGlobalResult("pdf", pdfz.m_postmargZResult);
  resultStore->StoreScopedGlobalResult("pdf_params", pdfz.m_postmargZResult);

  // save in resultstore candidates results
  resultStore->StoreScopedGlobalResult("candidatesresult", candidateResult);

  // for each extrema, get best model by reading from datastore and selecting
  // best fit
  /////////////////////////////////////////////////////////////////////////////////////
  if (fft_processing) {
    tplCatalog.m_logsampling = false;
  }
  std::shared_ptr<const ExtremaResult> extremaResult = buildExtremaResults(
      resultStore, scopeStr, candidateResult->m_ranked_candidates, tplCatalog,
      overlapThreshold);

  // store extrema results
  StoreExtremaResults(resultStore, extremaResult);

  std::shared_ptr<CTemplateFittingSolveResult> TemplateFittingSolveResult =
      std::make_shared<CTemplateFittingSolveResult>(
          extremaResult->m_ranked_candidates[0].second, m_opt_pdfcombination,
          pdfz.m_postmargZResult->valMargEvidenceLog);

  return TemplateFittingSolveResult;
}

void CTemplateFittingSolve::Solve(
    std::shared_ptr<COperatorResultStore> resultStore,
    const std::shared_ptr<const CTemplate> &tpl, Float64 overlapThreshold,
    std::vector<CMask> maskList, EType spctype, std::string opt_interp,
    bool opt_extinction, bool opt_dustFitting) {

  std::string scopeStr = "templatefitting";
  Int32 _ntype = 1;
  CSpectrum::EType _spctype = CSpectrum::nType_raw;
  CSpectrum::EType _spctypetab[3] = {CSpectrum::nType_raw,
                                     CSpectrum::nType_noContinuum,
                                     CSpectrum::nType_continuumOnly};

  // case: nType_all
  if (spctype == nType_all) {
    _ntype = 3;
  }

  std::vector<CSpectrum::EType> save_spcTypes;
  for (auto spc : Context.getSpectra())
    save_spcTypes.push_back(spc->GetType());
  const CSpectrum::EType save_tplType = tpl->GetType();

  for (Int32 i = 0; i < _ntype; i++) {
    if (spctype == nType_all) {
      _spctype = _spctypetab[i];
    } else {
      _spctype = static_cast<CSpectrum::EType>(spctype);
    }
    for (auto spc : Context.getSpectra())
      spc->SetType(_spctype);
    tpl->SetType(_spctype);
    tpl->setRebinInterpMethod(opt_interp);

    if (_spctype == CSpectrum::nType_continuumOnly) {
      // use continuum only
      scopeStr = "templatefitting_continuum";

    } else if (_spctype == CSpectrum::nType_raw) {
      // use full spectrum
      scopeStr = "templatefitting";

    } else if (_spctype == CSpectrum::nType_noContinuum) {
      // use spectrum without continuum
      scopeStr = "templatefitting_nocontinuum";
      //
      opt_dustFitting = false;
    } else {
      // unknown type
      THROWG(ErrorCode::INTERNAL_ERROR, "Unknown spectrum component");
    }

    // Compute merit function
    auto templateFittingResult =
        std::dynamic_pointer_cast<CTemplateFittingResult>(
            m_continuumFittingOperator->Compute(tpl, overlapThreshold,
                                                opt_interp, opt_extinction,
                                                opt_dustFitting));

    if (!templateFittingResult)
      THROWG(ErrorCode::INTERNAL_ERROR,
             "no results returned by templateFittingOperator");

    // Store results
    resultStore->StoreScopedPerTemplateResult(tpl, scopeStr.c_str(),
                                              templateFittingResult);
  }

  int i = 0;
  for (auto spc : Context.getSpectra())
    spc->SetType(save_spcTypes[i++]);
  tpl->SetType(save_tplType);
}

ChisquareArray CTemplateFittingSolve::BuildChisquareArray(
    std::shared_ptr<const COperatorResultStore> store,
    const std::string &scopeStr) const {
  ChisquareArray chisquarearray;

  Log.LogDetail("templatefittingsolver: building chisquare array");
  Log.LogDetail(Formatter()
                << "    templatefittingsolver: using results in scope: "
                << store->GetScopedName(scopeStr));

  TOperatorResultMap meritResults =
      store->GetScopedPerTemplateResult(scopeStr.c_str());

  chisquarearray.cstLog = -1;
  chisquarearray.zstep = m_redshiftStep;
  for (TOperatorResultMap::const_iterator it = meritResults.begin();
       it != meritResults.end(); ++it) {
    auto meritResult =
        std::dynamic_pointer_cast<const CTemplateFittingResult>((*it).second);
    Int32 nISM = -1;
    Int32 nIGM = -1;
    if (meritResult->ChiSquareIntermediate.size() > 0) {
      nISM = meritResult->ChiSquareIntermediate[0].size();
      if (meritResult->ChiSquareIntermediate[0].size() > 0) {
        nIGM = meritResult->ChiSquareIntermediate[0][0].size();
      }
    }
    if (chisquarearray.cstLog == -1) {
      chisquarearray.cstLog = meritResult->CstLog;
      Log.LogInfo(Formatter() << "templatefittingsolver: using cstLog = "
                              << chisquarearray.cstLog);
    } else if (chisquarearray.cstLog != meritResult->CstLog) {
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "cstLog values do not correspond: val1="
                         << chisquarearray.cstLog
                         << " != val2=" << meritResult->CstLog);
    }
    if (chisquarearray.redshifts.size() == 0) {
      chisquarearray.redshifts = meritResult->Redshifts;
    }

    CZPrior zpriorhelper;
    for (Int32 kism = 0; kism < nISM; kism++) {
      for (Int32 kigm = 0; kigm < nIGM; kigm++) {
        chisquarearray.zpriors.push_back(
            zpriorhelper.GetConstantLogZPrior(meritResult->Redshifts.size()));

        // correct chi2 for ampl. marg. if necessary: todo add switch, currently
        // deactivated
        chisquarearray.chisquares.emplace_back(
            meritResult->ChiSquareIntermediate.size(), DBL_MAX);
        TFloat64List &logLikelihoodCorrected = chisquarearray.chisquares.back();
        for (Int32 kz = 0; kz < meritResult->Redshifts.size(); kz++) {
          logLikelihoodCorrected[kz] =
              meritResult->ChiSquareIntermediate
                  [kz][kism]
                  [kigm]; // + resultXXX->ScaleMargCorrectionTplratios[][]?;
        }
        Log.LogDetail(
            Formatter()
            << "    templatefittingsolver: Pdfz combine - prepared merit #"
            << chisquarearray.chisquares.size() - 1
            << " for model : " << ((*it).first));
      }
    }
  }

  return chisquarearray;
}

std::shared_ptr<const ExtremaResult> CTemplateFittingSolve::buildExtremaResults(
    std::shared_ptr<const COperatorResultStore> store,
    const std::string &scopeStr, const TCandidateZbyRank &ranked_zCandidates,
    const CTemplateCatalog &tplCatalog, Float64 overlapThreshold) {

  Log.LogDetail(
      "CTemplateFittingSolve::buildExtremaResults: building chisquare array");
  Log.LogDetail(Formatter()
                << "    templatefittingsolve:r using results in scope: "
                << store->GetScopedName(scopeStr));

  TOperatorResultMap results = store->GetScopedPerTemplateResult(scopeStr);

  Int32 extremumCount = ranked_zCandidates.size();

  auto firstResult = std::dynamic_pointer_cast<const CTemplateFittingResult>(
      (*results.begin()).second);
  const TFloat64List &redshifts = firstResult->Redshifts;

  // check all results  status
  for (auto &r : results) {
    auto TplFitResult =
        std::dynamic_pointer_cast<const CTemplateFittingResult>(r.second);
    if (TplFitResult->ChiSquare.size() != redshifts.size()) {
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter()
                 << "Size do not match among templatefitting results, for tpl="
                 << r.first);
    }
    for (Int32 kz = 0; kz < TplFitResult->Redshifts.size(); kz++) {
      if (TplFitResult->Redshifts[kz] != redshifts[kz]) {
        THROWG(ErrorCode::INTERNAL_ERROR,
               Formatter() << "redshift vector is not the same for tpl="
                           << r.first);
      }
    }
  }

  std::shared_ptr<ExtremaResult> extremaResult =
      make_shared<ExtremaResult>(ranked_zCandidates);

  for (Int32 iExtremum = 0; iExtremum < extremumCount; iExtremum++) {
    // std::string Id = ranked_zCandidates[iExtremum].first;
    Float64 z = ranked_zCandidates[iExtremum].second->Redshift;

    // find the corresponding Z
    auto itZ = std::find(redshifts.begin(), redshifts.end(), z);
    const Int32 zIndex = std::distance(redshifts.begin(), itZ);

    // find the min chisquare at corresponding redshift
    Float64 ChiSquare = DBL_MAX;
    std::string tplName = "";
    for (auto &r : results) {
      auto TplFitResult =
          std::dynamic_pointer_cast<const CTemplateFittingResult>(r.second);

      if (TplFitResult->ChiSquare[zIndex] < ChiSquare) {
        ChiSquare = TplFitResult->ChiSquare[zIndex];
        tplName = r.first;
      };
    }

    // Fill extrema Result
    auto TplFitResult = std::dynamic_pointer_cast<const CTemplateFittingResult>(
        results[tplName]);
    extremaResult->m_ranked_candidates[iExtremum].second->fittedTpl.tplMerit =
        ChiSquare;
    extremaResult->m_ranked_candidates[iExtremum]
        .second->fittedTpl.tplMeritPhot = TplFitResult->ChiSquarePhot[zIndex];
    extremaResult->m_ranked_candidates[iExtremum].second->fittedTpl.tplName =
        tplName;
    extremaResult->m_ranked_candidates[iExtremum]
        .second->fittedTpl.tplMeiksinIdx = TplFitResult->FitMeiksinIdx[zIndex];
    extremaResult->m_ranked_candidates[iExtremum]
        .second->fittedTpl.tplEbmvCoeff = TplFitResult->FitEbmvCoeff[zIndex];
    extremaResult->m_ranked_candidates[iExtremum]
        .second->fittedTpl.tplAmplitude = TplFitResult->FitAmplitude[zIndex];
    extremaResult->m_ranked_candidates[iExtremum]
        .second->fittedTpl.tplAmplitudeError =
        TplFitResult->FitAmplitudeError[zIndex];
    extremaResult->m_ranked_candidates[iExtremum].second->fittedTpl.tplDtM =
        TplFitResult->FitDtM[zIndex];
    extremaResult->m_ranked_candidates[iExtremum].second->fittedTpl.tplMtM =
        TplFitResult->FitMtM[zIndex];
    extremaResult->m_ranked_candidates[iExtremum].second->fittedTpl.tplSNR =
        TplFitResult->SNR[zIndex];
    extremaResult->m_ranked_candidates[iExtremum]
        .second->fittedTpl.tplLogPrior = TplFitResult->LogPrior[zIndex];

    // make sure tpl is non-rebinned
    bool currentSampling = tplCatalog.m_logsampling;
    tplCatalog.m_logsampling = false;
    std::shared_ptr<const CTemplate> tpl =
        tplCatalog.GetTemplateByName({m_category}, tplName);

    std::shared_ptr<CModelSpectrumResult> spcmodelPtr =
        std::make_shared<CModelSpectrumResult>();
    for (int spcIndex = 0; spcIndex < Context.getSpectra().size(); spcIndex++) {
      const std::string &obsId = Context.getSpectra()[spcIndex]->getObsID();

      TPhotVal values = m_continuumFittingOperator->ComputeSpectrumModel(
          tpl, z, TplFitResult->FitEbmvCoeff[zIndex],
          TplFitResult->FitMeiksinIdx[zIndex],
          TplFitResult->FitAmplitude[zIndex], overlapThreshold, spcIndex,
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
