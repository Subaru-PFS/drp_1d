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
#include <boost/test/unit_test.hpp>

#include "RedshiftLibrary/method/classificationresult.h"
#include "RedshiftLibrary/method/classificationsolve.h"
#include "RedshiftLibrary/operator/operator.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/processflow/resultstore.h"
#include "tests/src/tool/inputContextLight.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(classificationSolve_test)

Float64 precision = 1e-12;

// JsonFile in string format
std::string jsonString =
    "{\"lambdaRange\" : [ 4630, 4815 ],"
    "\"smoothWidth\" : 0.0,"
    "\"templateCatalog\" : {"
    "\"continuumRemoval\" : {"
    "\"method\" : \"zero\","
    "\"medianKernelWidth\" : 75,"
    "\"medianEvenReflection\" : true}},"
    "\"continuumRemoval\" : {"
    "\"method\" : \"irregularSamplingMedian\","
    "\"medianKernelWidth\" : 400,"
    "\"medianEvenReflection\" : true,"
    "\"decompScales\" : 9},"
    "\"autoCorrectInput\" : false,"
    "\"lsf\" : {\"lsfType\" : \"gaussianConstantResolution\", \"resolution\" : "
    "4300},"
    "\"spectrumModels\" : [\"galaxy\", \"star\", \"qso\"],"
    "\"galaxy\" : {\"redshiftSolver\": { \"method\" : \"templateFittingSolve\" "
    "}}, "
    "\"star\" : {\"redshiftSolver\": { \"method\" : \"templateFittingSolve\" "
    "}}, "
    "\"qso\" : {\"redshiftSolver\": { \"method\" : \"templateFittingSolve\"}}}";
std::string type = "test";

class fixture_classificationSolveTest {
public:
  fixture_Context ctx;
  fixture_classificationSolveTest() {
    fillCatalog();
    ctx.reset();
    ctx.loadParameterStore(jsonString);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.addSpectrum(spc, LSF);
    ctx.initContext();
  }

  std::shared_ptr<CScopeStack> scopeStack = std::make_shared<CScopeStack>();
  std::shared_ptr<CPhotBandCatalog> photoBandCatalog =
      fixture_PhotoBandCatalog().photoBandCatalog;
  std::shared_ptr<CPhotometricData> photoData = fixture_PhotoData().photoData;

  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin =
      fixture_MeiskinCorrection().igmCorrectionMeiksin;
  std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti =
      fixture_CalzettiCorrection().ismCorrectionCalzetti;
  std::shared_ptr<CLSF> LSF =
      fixture_LSFGaussianConstantResolution(scopeStack).LSF;
  std::shared_ptr<CSpectrum> spc = fixture_SharedSpectrum().spc;
  std::shared_ptr<CTemplateCatalog> catalog =
      fixture_sharedTemplateCatalog().catalog;

  void fillCatalog() {
    catalog->Add(fixture_SharedGalaxyTemplate().tpl);
    catalog->m_logsampling = 1;
    catalog->Add(fixture_SharedGalaxyTemplate().tpl);
  }
};

BOOST_FIXTURE_TEST_CASE(compute_test, fixture_classificationSolveTest) {

  // CAutoScope scope_test(Context.m_ScopeStack, "param");
  auto const &resultStore = Context.GetResultStore();

  // create candidateZ
  std::shared_ptr<TCandidateZ> candidateZ = std::make_shared<TCandidateZ>();

  // create Pdf solve result for galaxy
  std::shared_ptr<CPdfSolveResult> result_in =
      std::make_shared<CPdfSolveResult>(type, candidateZ, "marg", 1.);

  // create OperatorResult
  resultStore->StoreGlobalResult("galaxy.redshiftSolver.templateFittingSolve",
                                 "solveResult", result_in);

  // create Pdf solve result for star
  result_in = std::make_shared<CPdfSolveResult>(type, candidateZ, "marg", 1.5);
  resultStore->StoreGlobalResult("star.redshiftSolver.templateFittingSolve",
                                 "solveResult", result_in);

  // qso has no results : compute function catch this and does not take count of
  // qso objects in classification
  // create Classification solve
  std::string spectrumModel = "galaxy";
  CAutoScope spectrumModel_autoscope(Context.m_ScopeStack, spectrumModel,
                                     ScopeType::SPECTRUMMODEL);
  CAutoScope stage_autoscope(Context.m_ScopeStack, "classification",
                             ScopeType::STAGE);
  CClassificationSolve cSolve;

  std::shared_ptr<CSolveResult> classifResult;
  BOOST_REQUIRE_NO_THROW(classifResult = cSolve.compute());
  resultStore->StoreGlobalResult("object.stage.method", "classification",
                                 classifResult);

  auto result_out = resultStore->GetClassificationResult(
      "object", "stage", "method", "classification");

  // check results
  Float64 logEvidence_ref_1 = 1.;
  Float64 logEvidence_ref_2 = 1.5;

  Float64 proba_ref_1 = exp(logEvidence_ref_1 - logEvidence_ref_2);
  Float64 proba_ref_2 = exp(logEvidence_ref_2 - logEvidence_ref_2);

  Float64 sum_ref = proba_ref_1 + proba_ref_2;

  BOOST_CHECK_CLOSE(proba_ref_1 / sum_ref, result_out->m_proba.at("galaxy"),
                    precision);
  BOOST_CHECK_CLOSE(proba_ref_2 / sum_ref, result_out->m_proba.at("star"),
                    precision);

  BOOST_CHECK_CLOSE(logEvidence_ref_1, result_out->m_evidences.at("galaxy"),
                    precision);
  BOOST_CHECK_CLOSE(logEvidence_ref_2, result_out->m_evidences.at("star"),
                    precision);

  BOOST_CHECK(result_out->m_TypeLabel == "star");

  // Test with one probabilitie undefined
  result_in = std::make_shared<CPdfSolveResult>(type, candidateZ, "marg", 1.);
  resultStore->reset();
  resultStore->StoreGlobalResult("galaxy.redshiftSolver.templateFittingSolve",
                                 "solveResult", result_in);

  result_in =
      std::make_shared<CPdfSolveResult>(type, candidateZ, "marg", -INFINITY);
  resultStore->StoreGlobalResult("star.redshiftSolver.templateFittingSolve",
                                 "solveResult", result_in);

  classifResult = cSolve.compute();
  resultStore->StoreGlobalResult("object.stage.method", "classification_2",
                                 classifResult);
  result_out = resultStore->GetClassificationResult("object", "stage", "method",
                                                    "classification_2");

  logEvidence_ref_1 = 1.;
  logEvidence_ref_2 = -INFINITY;

  BOOST_CHECK_CLOSE(1., result_out->m_proba.at("galaxy"), precision);
  BOOST_CHECK_CLOSE(0., result_out->m_proba.at("star"), precision);

  BOOST_CHECK_CLOSE(logEvidence_ref_1, result_out->m_evidences.at("galaxy"),
                    precision);
  BOOST_CHECK(logEvidence_ref_2 == result_out->m_evidences.at("star"));

  BOOST_CHECK(result_out->m_TypeLabel == "galaxy");

  // Test with all probabilitie undefined
  result_in =
      std::make_shared<CPdfSolveResult>(type, candidateZ, "marg", -INFINITY);
  resultStore->reset();
  resultStore->StoreGlobalResult("galaxy.redshiftSolver.templateFittingSolve",
                                 "solveResult", result_in);

  result_in =
      std::make_shared<CPdfSolveResult>(type, candidateZ, "marg", -INFINITY);
  resultStore->StoreGlobalResult("star.redshiftSolver.templateFittingSolve",
                                 "solveResult", result_in);

  BOOST_CHECK_THROW(cSolve.compute(), AmzException);
}

BOOST_AUTO_TEST_SUITE_END()