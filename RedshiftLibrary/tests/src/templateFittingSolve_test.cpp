
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
#include "RedshiftLibrary/method/templatefittingsolve.h"
#include "RedshiftLibrary/method/templatefittingsolveresult.h"
#include "RedshiftLibrary/operator/extremaresult.h"
#include "RedshiftLibrary/processflow/context.h"
#include "tests/src/tool/inputContextLight.h"
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

const std::string jsonString =
    "{\"lambdarange\" : [ 4680, 4712 ],"
    "\"smoothWidth\" : 0.0,"
    "\"templateCatalog\" : {"
    "\"continuumRemoval\" : {"
    "\"method\" : \"zero\","
    "\"medianKernelWidth\" : 75,"
    "\"medianEvenReflection\" : true}},"
    "\"ebmv\" : {\"start\" : 0, \"step\" : 0.1, \"count\" : 10},"
    "\"continuumRemoval\" : {"
    "\"method\" : \"IrregularSamplingMedian\","
    "\"medianKernelWidth\" : 400,"
    "\"medianEvenReflection\" : true,"
    "\"decompScales\" : 9},"
    "\"LSF\" : {\"LSFType\" : \"GaussianConstantResolution\", \"resolution\" : "
    "4300},"
    "\"extremaredshiftseparation\" : 0.01,"
    "\"objects\" : [\"galaxy\"],"
    "\"calibrationDir\" : \"\","
    "\"autocorrectinput\" : false,"
    "\"airvacuum_method\" : \"default\",";

const std::string jsonStringFFT = {
    "\"galaxy\" : {"
    "\"redshiftrange\" : [ 2.84, 2.88 ],"
    "\"redshiftstep\" : 0.0001,"
    "\"redshiftsampling\" : \"log\","
    "\"method\" : \"TemplateFittingSolve\","
    "\"linemeas_method\" : null,"
    "\"TemplateFittingSolve\" : {"
    "\"extremacount\" : 5,"
    "\"overlapThreshold\" : 1,"
    "\"spectrum\" : {\"component\" : \"raw\"},"
    "\"fftprocessing\" : true,"
    "\"interpolation\" : \"precomputedfinegrid\","
    "\"extinction\" : true,"
    "\"dustfit\" : true,"
    "\"pdfcombination\" : \"marg\","
    "\"enablephotometry\" : false}}}"};

const std::string jsonStringNoFFT = {
    "\"galaxy\" : {"
    "\"redshiftrange\" : [ 2.84, 2.88 ],"
    "\"redshiftstep\" : 0.0001,"
    "\"redshiftsampling\" : \"log\","
    "\"method\" : \"TemplateFittingSolve\","
    "\"linemeas_method\" : null,"
    "\"TemplateFittingSolve\" : {"
    "\"extremacount\" : 5,"
    "\"overlapThreshold\" : 1,"
    "\"spectrum\" : {\"component\" : \"raw\"},"
    "\"fftprocessing\" : false,"
    "\"interpolation\" : \"precomputedfinegrid\","
    "\"extinction\" : true,"
    "\"dustfit\" : true,"
    "\"pdfcombination\" : \"marg\","
    "\"enablephotometry\" : true,"
    "\"photometry\": {\"weight\" : 1.0}}}}"};

class fixture_TemplateFittingSolveTestNoFFT {
public:
  fixture_Context ctx;
  fixture_TemplateFittingSolveTestNoFFT() {
    fillCatalog();
    ctx.loadParameterStore(jsonString + jsonStringNoFFT);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.setPhotoBandCatalog(photoBandCatalog);
    spc->SetPhotData(photoData);
    ctx.addSpectrum(spc, LSF);
    ctx.initContext();
  }

  TScopeStack scopeStack;
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin =
      fixture_MeiskinCorrection().igmCorrectionMeiksin;
  std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti =
      fixture_CalzettiCorrection().ismCorrectionCalzetti;
  std::shared_ptr<CLSF> LSF =
      fixture_LSFGaussianConstantResolution(scopeStack).LSF;
  std::shared_ptr<CSpectrum> spc = fixture_SharedSpectrum().spc;
  std::shared_ptr<CTemplateCatalog> catalog =
      fixture_sharedTemplateCatalog().catalog;
  std::shared_ptr<CPhotBandCatalog> photoBandCatalog =
      fixture_PhotoBandCatalog().photoBandCatalog;
  std::shared_ptr<CPhotometricData> photoData = fixture_PhotoData().photoData;

  void fillCatalog() {
    catalog->Add(fixture_SharedGalaxyTemplate().tpl);
    catalog->m_logsampling = 1;
    catalog->Add(fixture_SharedGalaxyTemplate().tpl);
  }
};

class fixture_TemplateFittingSolveTestFFT {
public:
  fixture_Context ctx;
  fixture_TemplateFittingSolveTestFFT() {
    fillCatalog();
    ctx.loadParameterStore(jsonString + jsonStringFFT);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.addSpectrum(spc, LSF);
    ctx.initContext();
  }

  TScopeStack scopeStack;
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin =
      fixture_MeiskinCorrection().igmCorrectionMeiksin;
  std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti =
      fixture_CalzettiCorrection().ismCorrectionCalzetti;
  std::shared_ptr<CLSF> LSF =
      fixture_LSFGaussianConstantResolution(scopeStack).LSF;
  std::shared_ptr<CSpectrum> spc = fixture_SharedSpectrumExtended().spc;
  std::shared_ptr<CTemplateCatalog> catalog =
      fixture_sharedTemplateCatalog().catalog;

  void fillCatalog() {
    catalog->Add(fixture_SharedGalaxyTemplate().tpl);
    catalog->m_logsampling = 1;
    catalog->Add(fixture_SharedGalaxyTemplate().tpl);
  }
};

BOOST_AUTO_TEST_SUITE(templateFittingSolve_test)

BOOST_FIXTURE_TEST_CASE(computeNoFFT_test,
                        fixture_TemplateFittingSolveTestNoFFT) {
  CTemplateFittingSolve templateFittingSolve(Context.m_ScopeStack, "galaxy");
  BOOST_CHECK_NO_THROW(templateFittingSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("galaxy",
                                               "TemplateFittingSolve");
  BOOST_CHECK(result_out.lock()->getType() == "CTemplateFittingSolveResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "TemplateFittingSolve", "pdf");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "TemplateFittingSolve", "pdf_params");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  std::string resType = Context.GetResultStore()->GetCandidateResultType(
      "galaxy", "TemplateFittingSolve", "extrema_results", "model_parameters");
  BOOST_CHECK(resType == "TExtremaResult");

  std::shared_ptr<const TExtremaResult> res =
      Context.GetResultStore()->GetExtremaResult(
          "galaxy", "TemplateFittingSolve", "extrema_results",
          "model_parameters", 0);
  Float64 z = res->Redshift;
  BOOST_CHECK_CLOSE(z, 2.8604060282076706, 1e-6);

  ctx.reset();
}

BOOST_FIXTURE_TEST_CASE(computeFFT_test, fixture_TemplateFittingSolveTestFFT) {
  CTemplateFittingSolve templateFittingSolve(Context.m_ScopeStack, "galaxy");
  BOOST_CHECK_NO_THROW(templateFittingSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("galaxy",
                                               "TemplateFittingSolve");
  BOOST_CHECK(result_out.lock()->getType() == "CTemplateFittingSolveResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "TemplateFittingSolve", "pdf");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "TemplateFittingSolve", "pdf_params");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  std::string resType = Context.GetResultStore()->GetCandidateResultType(
      "galaxy", "TemplateFittingSolve", "extrema_results", "model_parameters");
  BOOST_CHECK(resType == "TExtremaResult");

  std::shared_ptr<const TExtremaResult> res =
      Context.GetResultStore()->GetExtremaResult(
          "galaxy", "TemplateFittingSolve", "extrema_results",
          "model_parameters", 0);

  Float64 z = res->Redshift;
  BOOST_CHECK_CLOSE(z, 2.8777553864462204, 1e-6);

  ctx.reset();
}

BOOST_AUTO_TEST_SUITE_END()
