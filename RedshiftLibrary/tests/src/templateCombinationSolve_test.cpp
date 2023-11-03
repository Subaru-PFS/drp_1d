
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
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/method/templatefittingsolveresult.h"
#include "RedshiftLibrary/method/tplcombinationsolve.h"
#include "RedshiftLibrary/operator/extremaresult.h"
#include "RedshiftLibrary/processflow/context.h"
#include "tests/src/tool/inputContextLight.h"
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

const std::string jsonString = "{\"lambdaRange\" : [ 7200, 7800 ],"
                               "\"smoothWidth\" : 0.0,"
                               "\"templateCatalog\" : {"
                               "\"continuumRemoval\" : {"
                               "\"method\" : \"irregularSamplingMedian\","
                               "\"medianKernelWidth\" : 75,"
                               "\"medianEvenReflection\" : true}},"
                               "\"continuumRemoval\" : {"
                               "\"method\" : \"irregularSamplingMedian\","
                               "\"medianKernelWidth\" : 400,"
                               "\"medianEvenReflection\" : true,"
                               "\"decompScales\" : 9},"
                               "\"extremaRedshiftSeparation\" : 0.01,"
                               "\"spectrumModels\" : [\"qso\"],"
                               "\"autoCorrectInput\" : false,"
                               "\"qso\" : {"
                               "\"redshiftRange\" : [ 5.13, 5.16 ],"
                               "\"redshiftStep\" : 0.001,"
                               "\"redshiftSampling\" : \"log\","
                               "\"method\" : \"tplCombinationSolver\","
                               "\"linemeas_method\" : null,"
                               "\"tplCombinationSolver\" : {"
                               "\"overlapThreshold\" : 1,"
                               "\"interpolation\" : \"lin\","
                               "\"pdfCombination\" : \"marg\","
                               "\"extremaCount\" : 5,"
                               "\"ismFit\" : true,"
                               "\"igmFit\" : true,";

const std::string jsonStringSpcComponentRaw =
    "\"spectrum\" : {\"component\" : \"raw\"}}}}";

const std::string jsonStringSpcComponentNoContinuum =
    "\"spectrum\" : {\"component\" : \"noContinuum\"}}}}";

const std::string jsonStringSpcComponentContinuum =
    "\"spectrum\" : {\"component\" : \"continuum\"}}}}";

const std::string jsonStringSpcComponentAll =
    "\"spectrum\" : {\"component\" : \"all\"}}}}";

class fixture_TplCombinationTest {
public:
  fixture_Context ctx;
  TScopeStack scopeStack;
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin =
      fixture_MeiskinCorrection().igmCorrectionMeiksin;
  std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti =
      fixture_CalzettiCorrection().ismCorrectionCalzetti;
  std::shared_ptr<CLSF> LSF =
      fixture_LSFGaussianConstantResolution(scopeStack).LSF;
  std::shared_ptr<CSpectrum> spc = fixture_SharedSpectrumQso().spc;
  std::shared_ptr<CTemplateCatalog> catalog =
      fixture_sharedTemplateCatalog().catalog;

  void fillCatalog() {
    catalog->Add(fixture_SharedQsoTemplate().tpl_c1);
    catalog->Add(fixture_SharedQsoTemplate().tpl_c2);
    catalog->Add(fixture_SharedQsoTemplate().tpl_c3);
    catalog->Add(fixture_SharedQsoTemplate().tpl_c4);
    catalog->Add(fixture_SharedQsoTemplate().tpl_mean);
  }
};

class fixture_TplCombinationTestRaw : public fixture_TplCombinationTest {
public:
  fixture_TplCombinationTestRaw() {
    fillCatalog();
    ctx.loadParameterStore(jsonString + jsonStringSpcComponentRaw);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.addSpectrum(spc, LSF);
    ctx.initContext();
  }
};

class fixture_TplCombinationTestNoContinuum
    : public fixture_TplCombinationTest {
public:
  fixture_TplCombinationTestNoContinuum() {
    fillCatalog();
    ctx.loadParameterStore(jsonString + jsonStringSpcComponentNoContinuum);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.addSpectrum(spc, LSF);
    ctx.initContext();
  }
};

class fixture_TplCombinationTestContinuum : public fixture_TplCombinationTest {
public:
  fixture_TplCombinationTestContinuum() {
    fillCatalog();
    ctx.loadParameterStore(jsonString + jsonStringSpcComponentContinuum);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.addSpectrum(spc, LSF);
    ctx.initContext();
  }
};

class fixture_TplCombinationTestAll : public fixture_TplCombinationTest {
public:
  fixture_TplCombinationTestAll() {
    fillCatalog();
    ctx.loadParameterStore(jsonString + jsonStringSpcComponentAll);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.addSpectrum(spc, LSF);
    ctx.initContext();
  }
};

BOOST_AUTO_TEST_SUITE(tplCombinationSolve_test)

BOOST_FIXTURE_TEST_CASE(computeRaw_test, fixture_TplCombinationTestRaw) {
  CTplcombinationSolve tplcombinationSolve(Context.m_ScopeStack, "qso");
  BOOST_CHECK_NO_THROW(tplcombinationSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("qso", "tplCombinationSolver");
  BOOST_CHECK(result_out.lock()->getType() == "CTplCombinationSolveResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "qso", "tplCombinationSolver", "pdf");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "qso", "tplCombinationSolver", "pdf_params");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  std::string resType = Context.GetResultStore()->GetCandidateResultType(
      "qso", "tplCombinationSolver", "extrema_results", "model_parameters");
  BOOST_CHECK(resType == "TTplCombinationResult");

  std::shared_ptr<const TExtremaResult> res =
      Context.GetResultStore()->GetExtremaResult("qso", "tplCombinationSolver",
                                                 "extrema_results",
                                                 "model_parameters", 0);

  Float64 z = res->Redshift;
  BOOST_CHECK_CLOSE(z, 5.1299999999999999, 1e-6);

  ctx.reset();
}

BOOST_FIXTURE_TEST_CASE(computeNoContinuum_test,
                        fixture_TplCombinationTestNoContinuum) {
  CTplcombinationSolve tplcombinationSolve(Context.m_ScopeStack, "qso");
  BOOST_CHECK_NO_THROW(tplcombinationSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("qso", "tplCombinationSolver");
  BOOST_CHECK(result_out.lock()->getType() == "CTplCombinationSolveResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "qso", "tplCombinationSolver", "pdf");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "qso", "tplCombinationSolver", "pdf_params");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  std::string resType = Context.GetResultStore()->GetCandidateResultType(
      "qso", "tplCombinationSolver", "extrema_results", "model_parameters");
  BOOST_CHECK(resType == "TTplCombinationResult");

  std::shared_ptr<const TExtremaResult> res =
      Context.GetResultStore()->GetExtremaResult("qso", "tplCombinationSolver",
                                                 "extrema_results",
                                                 "model_parameters", 0);

  Float64 z = res->Redshift;
  BOOST_CHECK_CLOSE(z, 5.1545691054521061, 1e-6);

  ctx.reset();
}

BOOST_FIXTURE_TEST_CASE(computeContinuum_test,
                        fixture_TplCombinationTestContinuum) {
  CTplcombinationSolve tplcombinationSolve(Context.m_ScopeStack, "qso");
  BOOST_CHECK_NO_THROW(tplcombinationSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("qso", "tplCombinationSolver");
  BOOST_CHECK(result_out.lock()->getType() == "CTplCombinationSolveResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "qso", "tplCombinationSolver", "pdf");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "qso", "tplCombinationSolver", "pdf_params");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  std::string resType = Context.GetResultStore()->GetCandidateResultType(
      "qso", "tplCombinationSolver", "extrema_results", "model_parameters");
  BOOST_CHECK(resType == "TTplCombinationResult");

  std::shared_ptr<const TExtremaResult> res =
      Context.GetResultStore()->GetExtremaResult("qso", "tplCombinationSolver",
                                                 "extrema_results",
                                                 "model_parameters", 0);

  Float64 z = res->Redshift;
  BOOST_CHECK_CLOSE(z, 5.1299999999999999, 1e-6);

  ctx.reset();
}

BOOST_FIXTURE_TEST_CASE(computeAll_test, fixture_TplCombinationTestAll) {
  CTplcombinationSolve tplcombinationSolve(Context.m_ScopeStack, "qso");
  BOOST_CHECK_NO_THROW(tplcombinationSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("qso", "tplCombinationSolver");
  BOOST_CHECK(result_out.lock()->getType() == "CTplCombinationSolveResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "qso", "tplCombinationSolver", "pdf");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "qso", "tplCombinationSolver", "pdf_params");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  std::string resType = Context.GetResultStore()->GetCandidateResultType(
      "qso", "tplCombinationSolver", "extrema_results", "model_parameters");
  BOOST_CHECK(resType == "TTplCombinationResult");

  std::shared_ptr<const TExtremaResult> res =
      Context.GetResultStore()->GetExtremaResult("qso", "tplCombinationSolver",
                                                 "extrema_results",
                                                 "model_parameters", 0);

  Float64 z = res->Redshift;
  BOOST_CHECK_CLOSE(z, 5.1299999999999999, 1e-6);

  ctx.reset();
}

BOOST_AUTO_TEST_SUITE_END()
