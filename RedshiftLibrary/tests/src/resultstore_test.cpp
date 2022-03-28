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
#include "RedshiftLibrary/continuum/indexesresult.h"
#include "RedshiftLibrary/linemodel/continuummodelsolution.h"
#include "RedshiftLibrary/linemodel/linemodelextremaresult.h"
#include "RedshiftLibrary/linemodel/linemodelsolution.h"
#include "RedshiftLibrary/method/classificationresult.h"
#include "RedshiftLibrary/method/reliabilityresult.h"
#include "RedshiftLibrary/operator/extremaresult.h"
#include "RedshiftLibrary/operator/flagResult.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/operator/operator.h"
#include "RedshiftLibrary/operator/pdfMargZLogResult.h"
#include "RedshiftLibrary/operator/spectraFluxResult.h"
#include "RedshiftLibrary/operator/tplCombinationExtremaResult.h"
#include "RedshiftLibrary/processflow/resultstore.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/statistics/pdfcandidatesz.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(ResultStore)

typedef std::map<std::string, std::shared_ptr<const COperatorResult>>
    TResultsMap;
typedef std::map<std::string, TResultsMap> TPerTemplateResultsMap;

//---------------------------------------------------------------
// RESULTS MODEL
//---------------------------------------------------------------

// Create scopeStack for tests
TScopeStack getScopeStack() {
  TScopeStack scopeStack;
  scopeStack.push_back("object");
  scopeStack.push_back("method");
  return scopeStack;
};

// Create Flag result
std::shared_ptr<const CFlagLogResult> getFlagResult() {
  std::shared_ptr<const CFlagLogResult> result_in =
      std::make_shared<const CFlagLogResult>(Flag.getBitMask(),
                                             Flag.getListMessages());
  return result_in;
}

// create candidates
TCandidateZbyRank getZCandidates() {
  TCandidateZ zcandidates;
  TCandidateZbyRank zcandByRank;
  zcandByRank.push_back(std::make_pair("EXT0", zcandidates));
  return zcandByRank;
}

// create LineModel Solution
std::shared_ptr<CLineModelSolution> getLineModelSolution() {
  std::shared_ptr<CLineModelSolution> result_in =
      std::make_shared<CLineModelSolution>();
  return result_in;
}

// create Model Spectrum Result>
std::shared_ptr<CModelSpectrumResult> getModelSpectrumResult() {
  CSpectrum spc;
  std::shared_ptr<CModelSpectrumResult> result_in =
      std::make_shared<CModelSpectrumResult>(spc);
  return result_in;
}

// create Spectra Flux result
std::shared_ptr<CSpectraFluxResult> getSpectraFluxresult() {
  std::shared_ptr<CSpectraFluxResult> result_in =
      (std::shared_ptr<CSpectraFluxResult>)new CSpectraFluxResult();
  return result_in;
}

// Create LineModel Extrema Result
std::shared_ptr<LineModelExtremaResult> getLineModelExtremaResult() {

  std::shared_ptr<LineModelExtremaResult> result_in =
      make_shared<LineModelExtremaResult>(getZCandidates());

  result_in->m_savedModelContinuumSpectrumResults[0] = getSpectraFluxresult();

  std::vector<std::shared_ptr<CModelSpectrumResult>> modelSpectrimResultList;
  std::shared_ptr<CModelSpectrumResult> modelSpectrimResult =
      getModelSpectrumResult();
  modelSpectrimResultList.push_back(modelSpectrimResult);
  result_in->m_savedModelSpectrumResults[0] = modelSpectrimResultList[0];

  std::vector<std::shared_ptr<CLineModelSolution>> lineModelresultList;
  std::shared_ptr<CLineModelSolution> lineModelresult = getLineModelSolution();
  lineModelresultList.push_back(lineModelresult);
  result_in->m_savedModelFittingResults[0] = lineModelresultList[0];

  return result_in;
}

// Create Tpl Combination Extrema Result
std::shared_ptr<TplCombinationExtremaResult> getTplCombinationExtremaResult() {
  std::shared_ptr<TplCombinationExtremaResult> result_in =
      make_shared<TplCombinationExtremaResult>(getZCandidates());
  return result_in;
}

// Create Extrema Result
std::shared_ptr<ExtremaResult> getExtremaResult() {
  std::shared_ptr<ExtremaResult> result_in =
      make_shared<ExtremaResult>(getZCandidates());
  return result_in;
}

//---------------------------------------------------------------
// TESTS
//---------------------------------------------------------------

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(StoreResult_test) {
  TScopeStack scopeStack_1 = getScopeStack();

  std::shared_ptr<ExtremaResult> result_in = getExtremaResult();

  COperatorResultStore store_1(scopeStack_1);

  // test store in context
  store_1.StoreResult(store_1.m_GlobalResults, store_1.GetCurrentScopeName(),
                      "extremaResult", result_in);
  BOOST_CHECK(store_1.GetScopedName("extremaResult") ==
              "object.method.extremaResult");

  std::shared_ptr<const TExtremaResult> result_out =
      store_1.GetExtremaResult("object", "method", "extremaResult", 0);
  BOOST_CHECK(result_out->getType() == "TExtremaResult");

  // test store outside context
  Flag.warning(Flag.CRANGE_NO_INTERSECTION, "Test code 4");

  TScopeStack scopeStack_2;
  COperatorResultStore store_2(scopeStack_2);

  store_2.StoreResult(store_2.m_GlobalResults, "", "warningFlag",
                      std::make_shared<const CFlagLogResult>(
                          Flag.getBitMask(), Flag.getListMessages()));
  BOOST_CHECK(store_2.GetScopedName("warningFlag") == "warningFlag");

  std::shared_ptr<const COperatorResult> result_out_2 =
      store_2.GetFlagResult("warningFlag", "", "warningFlag");
  BOOST_CHECK(result_out_2->getType() == "CFlagLogResult");
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(StoreGlobalResult_test) {
  TScopeStack scopeStack_1 = getScopeStack();

  std::shared_ptr<ExtremaResult> result_in = getExtremaResult();

  COperatorResultStore store_1(scopeStack_1);

  // test store in context
  store_1.StoreGlobalResult(store_1.GetCurrentScopeName(), "extremaResult",
                            result_in);
  BOOST_CHECK(store_1.GetScopedName("extremaResult") ==
              "object.method.extremaResult");

  std::shared_ptr<const TExtremaResult> result_out =
      store_1.GetExtremaResult("object", "method", "extremaResult", 0);
  BOOST_CHECK(result_out->getType() == "TExtremaResult");

  // test store outside context
  Flag.warning(Flag.CRANGE_NO_INTERSECTION, "Test code 4");

  TScopeStack scopeStack_2;
  COperatorResultStore store_2(scopeStack_2);

  store_2.StoreGlobalResult("", "warningFlag",
                            std::make_shared<const CFlagLogResult>(
                                Flag.getBitMask(), Flag.getListMessages()));
  BOOST_CHECK(store_2.GetScopedName("warningFlag") == "warningFlag");

  std::shared_ptr<const COperatorResult> result_out_2 =
      store_2.GetFlagResult("warningFlag", "", "warningFlag");
  BOOST_CHECK(result_out_2->getType() == "CFlagLogResult");
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(StoreFlagMethods_test) {
  TScopeStack scopeStack_1 = getScopeStack();

  Flag.warning(Flag.CRANGE_NO_INTERSECTION, "Test code 4");
  std::shared_ptr<const CFlagLogResult> result_in = getFlagResult();

  COperatorResultStore store_1(scopeStack_1);

  // test store flag in context
  store_1.StoreScopedGlobalResult("warningFlag", result_in);
  BOOST_CHECK(store_1.GetScopedName("warningFlag") ==
              "object.method.warningFlag");

  std::shared_ptr<const COperatorResult> result_out =
      store_1.GetFlagResult("object", "method", "warningFlag");
  BOOST_CHECK(result_out->getType() == "CFlagLogResult");

  // test store flag outside context
  TScopeStack scopeStack_2;
  COperatorResultStore store_2(scopeStack_2);

  store_2.StoreFlagResult("warningFlag", Flag.getBitMask());
  BOOST_CHECK(store_2.GetScopedName("warningFlag") == "warningFlag");

  result_out = store_2.GetFlagResult("warningFlag", "", "warningFlag");
  BOOST_CHECK(result_out->getType() == "CFlagLogResult");
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(StoreTemplateMethods_test) {
  TScopeStack scopeStack = getScopeStack();

  Flag.warning(Flag.CRANGE_NO_INTERSECTION, "Test code 4");

  std::shared_ptr<const CFlagLogResult> result_in = getFlagResult();

  std::shared_ptr<const CTemplate> tpl =
      std::make_shared<const CTemplate>("template_1", "category");

  // StoreScopedPerTemplateResult
  COperatorResultStore store(scopeStack);
  store.StoreScopedPerTemplateResult(tpl, "warningFlag", result_in);

  std::weak_ptr<const COperatorResult> result_out =
      store.GetPerTemplateResult(tpl, "warningFlag");
  BOOST_CHECK(result_out.lock()->getType() == "CFlagLogResult");
  BOOST_CHECK_THROW(store.GetPerTemplateResult(tpl, "warningFlag_"),
                    GlobalException);

  TOperatorResultMap result_out_2 = store.GetPerTemplateResult("warningFlag");
  TOperatorResultMap::const_iterator it;
  for (it = result_out_2.begin(); it != result_out_2.end(); ++it) {
    std::string tplName = (*it).first;
    BOOST_CHECK(result_out_2[tplName]->getType() == "CFlagLogResult");
  }

  // empty map
  result_out_2 = store.GetPerTemplateResult("warningFlag__");
  BOOST_CHECK(result_out_2.size() == 0);

  // test store and get outside context
  TScopeStack scopeStack_2;
  COperatorResultStore store_2(scopeStack_2);
  store_2.StoreScopedPerTemplateResult(tpl, "warningFlag", result_in);

  result_out = store_2.GetPerTemplateResult(tpl, "warningFlag");
  BOOST_CHECK(result_out.lock()->getType() == "CFlagLogResult");

  // StorePerTemplateResult
  store.StorePerTemplateResult(tpl, store.GetCurrentScopeName(),
                               "warningFlag_2", result_in);

  result_out = store.GetPerTemplateResult(tpl, "warningFlag_2");
  BOOST_CHECK(result_out.lock()->getType() == "CFlagLogResult");

  result_out_2 = store.GetPerTemplateResult("warningFlag_2");
  for (it = result_out_2.begin(); it != result_out_2.end(); ++it) {
    std::string tplName = (*it).first;
    BOOST_CHECK(result_out_2[tplName]->getType() == "CFlagLogResult");
  }
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(GetMethods_test) {
  TScopeStack scopeStack = getScopeStack();

  Flag.warning(Flag.CRANGE_NO_INTERSECTION, "Test code 4");

  std::shared_ptr<const CFlagLogResult> result_in = getFlagResult();

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("warningFlag", result_in);
  std::weak_ptr<const COperatorResult> result_out =
      store.GetScopedGlobalResult("warningFlag");
  BOOST_CHECK(result_out.lock()->getType() == "CFlagLogResult");
  BOOST_CHECK_THROW(result_out = store.GetGlobalResult("warningFlag__"),
                    GlobalException);

  result_out = store.GetGlobalResult("object", "method", "warningFlag");
  BOOST_CHECK(result_out.lock()->getType() == "CFlagLogResult");

  store.StoreScopedGlobalResult("warningFlag_2", result_in);
  result_out = store.GetGlobalResult("object", "method", "warningFlag_", 2);
  BOOST_CHECK(result_out.lock()->getType() == "CFlagLogResult");
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(GetGlobalResultType_test) {
  TScopeStack scopeStack = getScopeStack();

  std::shared_ptr<const CFlagLogResult> result_in = getFlagResult();

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("warningFlag", result_in);
  std::weak_ptr<const COperatorResult> result_out =
      store.GetScopedGlobalResult("warningFlag");
  BOOST_CHECK(store.GetGlobalResultType("object", "method", "warningFlag") ==
              "CFlagLogResult");
  BOOST_CHECK_THROW(
      store.GetGlobalResultType("object", "method", "warningFlag__"),
      GlobalException);
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(GetSolveResult_test) {
  TScopeStack scopeStack = getScopeStack();

  std::shared_ptr<const CFlagLogResult> result_in = getFlagResult();

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("solveResult", result_in);
  std::weak_ptr<const COperatorResult> result_out =
      store.GetSolveResult("object", "method");
  BOOST_CHECK(result_out.lock()->getType() == "CFlagLogResult");
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(GetScope_test) {
  TScopeStack scopeStack = getScopeStack();

  std::shared_ptr<const CFlagLogResult> result_in = getFlagResult();

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("warningFlag", result_in);

  std::string scope = store.GetScope(*result_in.get());
  BOOST_CHECK(scope == "object.method.");
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(GetClassificationResult_test) {
  TScopeStack scopeStack = getScopeStack();

  std::shared_ptr<CClassificationResult> result_in =
      std::make_shared<CClassificationResult>();
  result_in->SetTypeLabel("label");

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("classification", result_in);

  std::shared_ptr<const CClassificationResult> result_out =
      store.GetClassificationResult("object", "method", "classification");
  BOOST_CHECK(result_out->getType() == "CClassificationResult");
  BOOST_CHECK(result_out->m_TypeLabel == "label");
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(GetReliabilityResult_test) {
  TScopeStack scopeStack = getScopeStack();

  std::shared_ptr<const CReliabilityResult> result_in =
      std::make_shared<const CReliabilityResult>();

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("reliability", result_in);

  std::shared_ptr<const CReliabilityResult> result_out =
      store.GetReliabilityResult("object", "method", "reliability");
  BOOST_CHECK(result_out->getType() == "CReliabilityResult");
  BOOST_CHECK(result_out->m_ReliabilityLabel == "C6");
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(GetPdfMargZLogResult_test) {
  TScopeStack scopeStack = getScopeStack();

  TFloat64List redshifts = {1, 2, 3, 4};
  std::shared_ptr<const CPdfMargZLogResult> result_in =
      std::make_shared<const CPdfMargZLogResult>(redshifts);

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("pdfMarg", result_in);

  std::shared_ptr<const CPdfMargZLogResult> result_out =
      store.GetPdfMargZLogResult("object", "method", "pdfMarg");
  BOOST_CHECK(result_out->getType() == "CPdfMargZLogResult");
  BOOST_CHECK(result_out->getIndex(2) == 1);
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(GetLineModelResult_test) {
  TScopeStack scopeStack = getScopeStack();

  std::shared_ptr<LineModelExtremaResult> result_in =
      getLineModelExtremaResult();

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("lineModel", result_in);

  std::shared_ptr<const TLineModelResult> result_out =
      store.GetLineModelResult("object", "method", "lineModel", 0);
  BOOST_CHECK(result_out->getType() == "TLineModelResult");
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(GetTplCombinationResult_test) {
  TScopeStack scopeStack = getScopeStack();

  std::shared_ptr<TplCombinationExtremaResult> result_in =
      getTplCombinationExtremaResult();

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("tplCombination", result_in);

  std::shared_ptr<const TTplCombinationResult> result_out =
      store.GetTplCombinationResult("object", "method", "tplCombination", 0);
  BOOST_CHECK(result_out->getType() == "TTplCombinationResult");
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(GetExtremaResult_test) {
  TScopeStack scopeStack = getScopeStack();

  std::shared_ptr<ExtremaResult> result_in = getExtremaResult();

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("extremaResult", result_in);

  std::shared_ptr<const TExtremaResult> result_out =
      store.GetExtremaResult("object", "method", "extremaResult", 0);
  BOOST_CHECK(result_out->getType() == "TExtremaResult");
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(GetLineModelSolution_test) {
  TScopeStack scopeStack = getScopeStack();

  std::shared_ptr<CLineModelSolution> result_in = getLineModelSolution();

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("lineModel", result_in);

  std::shared_ptr<const CLineModelSolution> result_out =
      store.GetLineModelSolution("object", "method", "lineModel");
  BOOST_CHECK(result_out->getType() == "CLineModelSolution");
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(GetSpectraFluxResult_test) {
  TScopeStack scopeStack = getScopeStack();

  std::shared_ptr<LineModelExtremaResult> result_in =
      getLineModelExtremaResult();

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("spectraFlux", result_in);

  std::shared_ptr<const CSpectraFluxResult> result_out =
      store.GetSpectraFluxResult("object", "method", "spectraFlux", 0);
  BOOST_CHECK(result_out->getType() == "CSpectraFluxResult");
}

//---------------------------------------------------------------v
BOOST_AUTO_TEST_CASE(GetModelSpectrumResult_test) {
  TScopeStack scopeStack = getScopeStack();

  std::shared_ptr<LineModelExtremaResult> result_in =
      getLineModelExtremaResult();

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("modelSpectrum", result_in);

  std::shared_ptr<const CModelSpectrumResult> result_out =
      store.GetModelSpectrumResult("object", "method", "modelSpectrum", 0);
  BOOST_CHECK(result_out->getType() == "CModelSpectrumResult");
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(GetLineModelSolution_test_2) {
  TScopeStack scopeStack = getScopeStack();

  std::shared_ptr<LineModelExtremaResult> result_in =
      getLineModelExtremaResult();

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("lineModel", result_in);

  std::shared_ptr<const CLineModelSolution> result_out =
      store.GetLineModelSolution("object", "method", "lineModel", 0);

  BOOST_CHECK(result_out->getType() == "CLineModelSolution");
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(GetModelSpectrumResult_test_2) {
  TScopeStack scopeStack = getScopeStack();

  CSpectrum spc;
  std::shared_ptr<const CModelSpectrumResult> result_in =
      getModelSpectrumResult();

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("modelSpectrum", result_in);

  std::shared_ptr<const CModelSpectrumResult> result_out =
      store.GetModelSpectrumResult("object", "method", "modelSpectrum");
  BOOST_CHECK(result_out->getType() == "CModelSpectrumResult");
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(GetCandidateResultType_test) {
  TScopeStack scopeStack = getScopeStack();

  std::shared_ptr<const TplCombinationExtremaResult> result_in =
      getTplCombinationExtremaResult();

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("tplCombination", result_in);

  std::string result_out = store.GetCandidateResultType(
      "object", "method", "tplCombination", "model_parameters");
  BOOST_CHECK(result_out == "TTplCombinationResult");
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(HasCandidateDataset_test) {
  TScopeStack scopeStack = getScopeStack();

  std::shared_ptr<const TplCombinationExtremaResult> result_in =
      getTplCombinationExtremaResult();

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("tplCombination", result_in);

  bool result_out = store.HasCandidateDataset(
      "object", "method", "tplCombination", "model_parameters");
  BOOST_CHECK(result_out == true);

  result_out = store.HasCandidateDataset("object", "method", "tplCombination_2",
                                         "model");
  BOOST_CHECK(result_out == false);
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(HasDataset_test) {
  TScopeStack scopeStack = getScopeStack();

  TFloat64List redshifts = {1, 2, 3, 4};
  std::shared_ptr<const CPdfMargZLogResult> result_in =
      std::make_shared<const CPdfMargZLogResult>(redshifts);

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("pdfMarg", result_in);

  bool result_out = store.HasDataset("object", "method", "pdfMarg");
  BOOST_CHECK(result_out == true);

  result_out = store.HasDataset("object", "method", "pdfMarg_2");
  BOOST_CHECK(result_out == false);
}

//---------------------------------------------------------------
BOOST_AUTO_TEST_CASE(getNbRedshiftCandidates_test) {
  TScopeStack scopeStack = getScopeStack();

  TFloat64List redshifts = {1, 2, 3, 4};
  std::shared_ptr<CPdfMargZLogResult> result_in =
      std::make_shared<CPdfMargZLogResult>(redshifts);

  COperatorResultStore store(scopeStack);
  store.StoreScopedGlobalResult("extrema_results", result_in);

  int result_out = store.getNbRedshiftCandidates("object", "method");
  BOOST_CHECK(result_out == 0);

  std::shared_ptr<PdfCandidatesZResult> result_in_2 =
      make_shared<PdfCandidatesZResult>();
  result_in_2->m_ranked_candidates = getZCandidates();

  store.StoreScopedGlobalResult("extrema_results", result_in_2);

  result_out = store.getNbRedshiftCandidates("object", "method");
  BOOST_CHECK(result_out == 1);

  std::shared_ptr<ExtremaResult> result_in_3 = getExtremaResult();

  store.StoreScopedGlobalResult("extrema_results", result_in_3);
  result_out = store.getNbRedshiftCandidates("object", "method");
  BOOST_CHECK(result_out == 1);

  std::shared_ptr<LineModelExtremaResult> result_in_4 =
      getLineModelExtremaResult();

  store.StoreScopedGlobalResult("extrema_results", result_in_4);
  result_out = store.getNbRedshiftCandidates("object", "method");
  BOOST_CHECK(result_out == 1);

  std::shared_ptr<TplCombinationExtremaResult> result_in_5 =
      getTplCombinationExtremaResult();

  store.StoreScopedGlobalResult("extrema_results", result_in_5);
  result_out = store.getNbRedshiftCandidates("object", "method");
  BOOST_CHECK(result_out == 1);
}

//---------------------------------------------------------------
// CreateResultStorage not used
BOOST_AUTO_TEST_CASE(CreateResultStorage_test) {
  TScopeStack scopeStack = getScopeStack();

  std::fstream stream;
  const boost::filesystem::path path("/test/test.txt");
  const boost::filesystem::path baseDir(".");

  COperatorResultStore store(scopeStack);
  Int32 result_out = store.CreateResultStorage(stream, path, baseDir);

  BOOST_CHECK(result_out == 0 || result_out == 1);
}

//---------------------------------------------------------------
// getProcessedObjectTypes ???

BOOST_AUTO_TEST_SUITE_END()