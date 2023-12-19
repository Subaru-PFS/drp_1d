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
#include "RedshiftLibrary/log/consolehandler.h"
#include "RedshiftLibrary/log/filehandler.h"
#include "RedshiftLibrary/method/linemodelsolve.h"
#include "RedshiftLibrary/processflow/context.h"
#include "tests/src/tool/inputContextLight.h"
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

const std::string lambdaString = "{\"lambdaRange\" : [ 4680, 4712 ],";
const std::string multiLambdaString =
    "{\"lambdaRange\" :{ \"A\":[ 4680, 4697 ],\"B\":[4695,4712]}, ";

const std::string jsonString =
    "\"smoothWidth\" : 0.0,"
    "\"templateCatalog\" : {"
    "\"continuumRemoval\" : {"
    "\"method\" : \"zero\","
    "\"medianKernelWidth\" : 75,"
    "\"medianEvenReflection\" : true}},"
    "\"ebmv\" : {\"start\" : 0, \"step\" : 0.1, \"count\" : 10},"
    "\"continuumRemoval\" : {"
    "\"method\" : \"irregularSamplingMedian\","
    "\"medianKernelWidth\" : 400,"
    "\"medianEvenReflection\" : true,"
    "\"decompScales\" : 9},"
    "\"lsf\" : {\"lsfType\" : \"gaussianConstantResolution\", "
    "\"resolution\" : "
    "4300},"
    "\"extremaRedshiftSeparation\" : 0.01,"
    "\"spectrumModels\" : [\"galaxy\"],"
    "\"autoCorrectInput\" : false,"
    "\"airVacuumMethod\" : \"default\","
    "\"galaxy\" : {"
    "\"redshiftRange\" : [ 0.24, 0.3 ],"
    "\"redshiftStep\" : 0.0001,"
    "\"redshiftSampling\" : \"log\","
    "\"stages\" : [\"redshiftSolver\"],"
    "\"redshiftSolver\" : {"
    "\"method\" : \"lineModelSolve\","
    "\"lineModelSolve\" : {"
    "\"lineModel\" : {"
    "\"continuumReestimation\" : \"no\","
    "\"velocityFit\" : true,"
    "\"emVelocityFitMin\" : 100,"
    "\"emVelocityFitMax\" : 700, "
    "\"emVelocityFitStep\" : 20,"
    "\"absVelocityFitMin\" : 150,"
    "\"absVelocityFitMax\" : 500, "
    "\"absVelocityFitStep\" : 50,"
    "\"ampOffsetFit\": \"false\","
    "\"lbdaOffsetFit\": \"false\","
    "\"extremaCount\" : 5,"
    "\"extremaCountB\" : 3,"
    "\"nSigmaSupport\" : 8,"
    "\"hAlphaPrior\" : 0.5,"
    "\"nOfZPriorStrength\" : 1.0,"
    "\"extremaCutProbaThreshold\" : -1,"
    "\"skipSecondPass\" : false,"
    "\"secondPassLcFittingMethod\" : -1,"
    "\"useLogLambdaSampling\": false,"
    "\"lyaForceFit\": false,"
    "\"lyaForceDisableFit\": false,"
    "\"strongLinesPrior\" : 1.0,"
    "\"fittingMethod\": \"svd\","
    "\"lineWidthType\": \"combined\","
    "\"velocityEmission\" : 100,"
    "\"velocityAbsorption\": 100,"
    "\"lineTypeFilter\" : \"no\","
    "\"lineForceFilter\" : \"no\","
    "\"lyaFit\": {"
    "\"asymFitMin\" : 0,"
    "\"asymFitMax\" : 4, \"asymFitStep\" : 1, "
    "\"widthFitMin\" : 1,"
    "\"widthFitMax\" : 4, \"widthFitStep\" : 1, "
    "\"deltaFitMin\" : 0,"
    "\"deltaFitMax\" : 0, \"deltaStepMax\" : 1}, "
    "\"tplRatio\": { \"priors\": {"
    "\"betaA\" : 1,    \"betaTE\" : 1, \"betaZ\" : 1, "
    "\"catalogDirPath\" : \"\"}}, "
    "\"firstPass\": { \"fittingMethod\" : \"individual\", "
    "\"tplRatioIsmFit\" : true,"
    "\"largeGridStepRatio\" : 6, "
    "\"multipleContinuumFitDisable\": true},"
    "\"secondPass\" : {\"halfWindowSize\" : 0.001, "
    "\"continuumFit\" : \"reFitFirstPass\"},";

const std::string jsonStringS =
    "\"smoothWidth\" : 0.0,"
    "\"templateCatalog\" : {"
    "\"continuumRemoval\" : {"
    "\"method\" : \"zero\","
    "\"medianKernelWidth\" : 75,"
    "\"medianEvenReflection\" : true}},"
    "\"ebmv\" : {\"start\" : 0, \"step\" : 0.1, \"count\" : 10},"
    "\"continuumRemoval\" : {"
    "\"method\" : \"irregularSamplingMedian\","
    "\"medianKernelWidth\" : 400,"
    "\"medianEvenReflection\" : true,"
    "\"decompScales\" : 9},"
    "\"lsf\" : {\"lsfType\" : \"gaussianConstantResolution\", "
    "\"resolution\" : "
    "4300},"
    "\"extremaRedshiftSeparation\" : 0.01,"
    "\"spectrumModels\" : [\"galaxy\"],"
    "\"autoCorrectInput\" : false,"
    "\"airVacuumMethod\" : \"default\","
    "\"galaxy\" : {"
    "\"redshiftRange\" : [ 0.24, 0.3 ],"
    "\"redshiftStep\" : 0.0001,"
    "\"redshiftSampling\" : \"log\","
    "\"method\" : \"lineModelSolver\","
    "\"lineModelSolver\" : {"
    "\"lineModel\" : {"
    "\"continuumReestimation\" : \"no\","
    "\"velocityFit\" : true,"
    "\"emVelocityFitMin\" : 100,"
    "\"emVelocityFitMax\" : 200, "
    "\"emVelocityFitStep\" : 50,"
    "\"absVelocityFitMin\" : 150,"
    "\"absVelocityFitMax\" : 200, "
    "\"absVelocityFitStep\" : 50,"
    "\"ampOffsetFit\": \"false\","
    "\"lbdaOffsetFit\": \"false\","
    "\"extremaCount\" : 5,"
    "\"extremaCountB\" : 3,"
    "\"nSigmaSupport\" : 8,"
    "\"hAlphaPrior\" : 0.5,"
    "\"nOfZPriorStrength\" : 1.0,"
    "\"extremaCutProbaThreshold\" : -1,"

    "\"secondPassLcFittingMethod\" : -1,"
    "\"useLogLambdaSampling\": false,"
    "\"lyaForceFit\": false,"
    "\"lyaForceDisableFit\": false,"
    "\"strongLinesPrior\" : 1.0,"
    "\"fittingMethod\": \"individual\","
    "\"lineWidthType\": \"combined\","
    "\"velocityEmission\" : 100,"
    "\"velocityAbsorption\": 100,"
    "\"lineTypeFilter\" : \"no\","
    "\"lineForceFilter\" : \"no\","
    "\"lyaFit\": {"
    "\"asymFitMin\" : 0,"
    "\"asymFitMax\" : 4, \"asymFitStep\" : 1, "
    "\"widthFitMin\" : 1,"
    "\"widthFitMax\" : 4, \"widthFitStep\" : 1, "
    "\"deltaFitMin\" : 0,"
    "\"deltaFitMax\" : 0, \"deltaStepMax\" : 1}, "
    "\"tplRatio\": { \"priors\": {"
    "\"betaA\" : 1,    \"betaTE\" : 1, \"betaZ\" : 1, "
    "\"catalogDirPath\" : \"\"}}, "
    "\"firstPass\": { \"fittingMethod\" : \"individual\", "
    "\"tplRatio_ismfit\" : true,"
    "\"largeGridStepRatio\" : 10, "
    "\"multipleContinuumFitDisable\": true},"
    "\"secondPass\" : {\"halfWindowSize\" : 0.001, "
    "\"continuumFit\" : \"reFitFirstPass\"},";

const std::string jsonStringTplFitRules =
    "\"skipSecondPass\" : false,"
    "\"continuumComponent\" : \"tplFit\","
    "\"pdfCombination\" : \"marg\","
    "\"tplRatioIsmFit\" : true,"
    "\"rules\" : \"all\","
    "\"improveBalmerFit\" : true,"
    "\"lineRatioType\": \"rules\","
    "\"enablePhotometry\" : true, "
    "\"photometry\" : {\"weight\" : 1},"
    "\"continuumFit\" : { \"ignoreLineSupport\": false,"
    "\"negativeThreshold\": -5.0,"
    "\"count\" : 1,"
    "\"nullThreshold\": 3,"
    "\"ismFit\" : true,"
    "\"igmFit\" : true,"
    "\"fftProcessing\": false, "
    "\"priors\": { \"betaA\" : 1, \"betaTE\" : 1, \"betaZ\" : 1,"
    "\"catalogDirPath\" : \"\"}}}}}}}";

const std::string jsonStringTplFitTplRatio =
    "\"skipSecondPass\" : false,"
    "\"continuumComponent\" : \"tplFit\","
    "\"pdfCombination\" : \"marg\","
    "\"tplRatioIsmFit\" : true,"
    "\"rules\" : \"all\","
    "\"improveBalmerFit\" : true,"
    "\"lineRatioType\": \"tplRatio\","
    "\"enablePhotometry\" : true, "
    "\"photometry\" : {\"weight\" : 1},"
    "\"continuumFit\" : { \"ignoreLineSupport\": false,"
    "\"negativeThreshold\": -5.0,"
    "\"count\" : 1,"
    "\"nullThreshold\": 3,"
    "\"ismFit\" : true,"
    "\"igmFit\" : true,"
    "\"fftProcessing\": false, "
    "\"priors\": { \"betaA\" : 1, \"betaTE\" : 1, \"betaZ\" : 1,"
    "\"catalogDirPath\" : \"\"}}}}}}}";

const std::string jsonStringNoContinuumTplRatio =
    "\"skipSecondPass\" : true,"
    "\"continuumComponent\" : \"nocontinuum\","
    "\"pdfCombination\" : \"marg\","
    "\"tplRatioIsmFit\" : true,"
    "\"rules\" : \"all\","
    "\"improveBalmerFit\" : true,"
    "\"lineRatioType\": \"tplratio\","
    "\"enablePhotometry\" : false, "
    "\"photometry\" : {\"weight\" : 1},"
    "\"continuumFit\" : { \"ignoreLineSupport\": false,"
    "\"negativeThreshold\": -5.0,"
    "\"count\" : 1,"
    "\"nullThreshold\": 3,"
    "\"ismFit\" : true,"
    "\"igmFit\" : true,"
    "\"fftProcessing\": false, "
    "\"priors\": { \"betaA\" : 1, \"betaTE\" : 1, \"betaZ\" : 1,"
    "\"catalogDirPath\" : \"\"}}}}}}";

const std::string jsonStringFromSpectrum =
    "\"skipSecondPass\" : false,"
    "\"continuumComponent\" : \"fromSpectrum\","
    "\"pdfCombination\" : \"bestChi2\","
    "\"tplRatioIsmFit\" : true,"
    "\"rules\" : \"balmerSingle\","
    "\"improveBalmerFit\" : true,"
    "\"lineRatioType\": \"tplRatio\","
    "\"continuumFit\" : { \"ignoreLineSupport\": false,"
    "\"negativeThreshold\": -5.0,"
    "\"count\" : 1,"
    "\"nullThreshold\": 3,"
    "\"ismFit\" : true,"
    "\"igmFit\" : true,"
    "\"fftProcessing\": true, "
    "\"priors\": { \"betaA\" : 1, \"betaTE\" : 1, \"betaZ\" : 1,"
    "\"catalogDirPath\" : \"\"}}}}}}}";

class fixture_LineModelSolveTest {
public:
  fixture_Context ctx;
  TScopeStack scopeStack;
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin =
      fixture_MeiskinCorrection().igmCorrectionMeiksin;
  std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti =
      fixture_CalzettiCorrection().ismCorrectionCalzetti;
  std::shared_ptr<CLSF> LSF =
      fixture_LSFGaussianConstantResolution(scopeStack).LSF;
  std::shared_ptr<CSpectrum> spc = fixture_SharedSpectrumExtended().spc;
  std::shared_ptr<CSpectrum> spcA = fixture_SharedMultiSpectrum().spcA;
  std::shared_ptr<CSpectrum> spcB = fixture_SharedMultiSpectrum().spcB;
  std::shared_ptr<CTemplateCatalog> catalog =
      fixture_sharedTemplateCatalog().catalog;
  std::shared_ptr<CPhotBandCatalog> photoBandCatalog =
      fixture_PhotoBandCatalog().photoBandCatalog;
  std::shared_ptr<CPhotometricData> photoData = fixture_PhotoData().photoData;
  std::shared_ptr<CLineCatalogsTplRatio> lineRatioTplCatalog =
      fixture_LineRatioTplCatalog().lineRatioTplCatalog;
  std::shared_ptr<CLineRatioCatalog> lineRatioCatalog =
      fixture_LineRatioCatalog().lineRatioCatalog;
  std::shared_ptr<CLineCatalog> lineCatalog = fixture_LineCatalog().lineCatalog;

  void fillCatalog() {
    catalog->Add(fixture_SharedGalaxyTemplate().tpl2);
    catalog->Add(fixture_SharedGalaxyTemplate().tpl2);
  }
};

class fixture_LineModelSolveTestTplFitRules
    : public fixture_LineModelSolveTest {
public:
  fixture_LineModelSolveTestTplFitRules() {
    fillCatalog();
    ctx.loadParameterStore(lambdaString + jsonString + jsonStringTplFitRules);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.setPhotoBandCatalog(photoBandCatalog);
    spc->SetPhotData(photoData);
    ctx.addSpectrum(spc, LSF);
    ctx.setLineRatioCatalogCatalog("galaxy", lineRatioTplCatalog);
    ctx.setLineCatalog("galaxy", "lineModelSolve", lineCatalog);
    ctx.initContext();
    lineRatioTplCatalog->addLineRatioCatalog(*lineRatioCatalog);
  }
};

class fixture_LineModelSolveTestTplFitTplRatio
    : public fixture_LineModelSolveTest {
public:
  fixture_LineModelSolveTestTplFitTplRatio() {
    fillCatalog();
    ctx.loadParameterStore(lambdaString + jsonString +
                           jsonStringTplFitTplRatio);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.setPhotoBandCatalog(photoBandCatalog);
    spc->SetPhotData(photoData);
    ctx.addSpectrum(spc, LSF);
    ctx.setLineRatioCatalogCatalog("galaxy", lineRatioTplCatalog);
    ctx.setLineCatalog("galaxy", "lineModelSolve", lineCatalog);
    ctx.initContext();
    lineRatioTplCatalog->addLineRatioCatalog(*lineRatioCatalog);
  }
};

class fixture_LineModelSolveTestNoContTplRatio
    : public fixture_LineModelSolveTest {
public:
  fixture_LineModelSolveTestNoContTplRatio() {
    fillCatalog();
    ctx.loadParameterStore(lambdaString + jsonStringS +
                           jsonStringNoContinuumTplRatio);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.setPhotoBandCatalog(photoBandCatalog);
    spc->SetPhotData(photoData);
    ctx.addSpectrum(spc, LSF);
    ctx.setLineRatioCatalogCatalog("galaxy", lineRatioTplCatalog);
    ctx.setLineCatalog("galaxy", "lineModelSolver", lineCatalog);
    ctx.initContext();
    lineRatioTplCatalog->addLineRatioCatalog(*lineRatioCatalog);
  }
};

class fixture_LineModelSolveTestMultiNoContTplRatio
    : public fixture_LineModelSolveTest {
public:
  fixture_LineModelSolveTestMultiNoContTplRatio() {
    fillCatalog();
    ctx.loadParameterStore(multiLambdaString + jsonStringS +
                           jsonStringTplFitTplRatio);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.setPhotoBandCatalog(photoBandCatalog);

    spcA->setObsID("A");
    spcB->setObsID("B");
    ctx.addSpectrum(spcA, LSF);
    ctx.addSpectrum(spcB, LSF);
    ctx.setLineRatioCatalogCatalog("galaxy", lineRatioTplCatalog);
    ctx.setLineCatalog("galaxy", "lineModelSolver", lineCatalog);
    ctx.initContext();
    lineRatioTplCatalog->addLineRatioCatalog(*lineRatioCatalog);
  }
};

class fixture_LineModelSolveTestFromSpectrum
    : public fixture_LineModelSolveTest {
public:
  fixture_LineModelSolveTestFromSpectrum() {
    fillCatalog();
    ctx.loadParameterStore(lambdaString + jsonString + jsonStringFromSpectrum);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.setPhotoBandCatalog(photoBandCatalog);
    spc->SetPhotData(photoData);
    ctx.addSpectrum(spc, LSF);
    ctx.setLineRatioCatalogCatalog("galaxy", lineRatioTplCatalog);
    ctx.setLineCatalog("galaxy", "lineModelSolve", lineCatalog);
    ctx.initContext();
    lineRatioTplCatalog->addLineRatioCatalog(*lineRatioCatalog);
  }
};

BOOST_AUTO_TEST_SUITE(lineModelSolve_test)

BOOST_FIXTURE_TEST_CASE(computeTplFitRules_test,
                        fixture_LineModelSolveTestTplFitRules) {

  CLineModelSolve lineModelSolve(Context.m_ScopeStack, "galaxy");
  BOOST_CHECK_NO_THROW(lineModelSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("galaxy", "redshiftSolver",
                                               "lineModelSolve");
  BOOST_CHECK(result_out.lock()->getType() == "CLineModelSolveResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "redshiftSolver", "lineModelSolve", "pdf");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "redshiftSolver", "lineModelSolve", "pdf_params");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  std::string resType = Context.GetResultStore()->GetCandidateResultType(
      "galaxy", "redshiftSolver", "lineModelSolve", "extrema_results",
      "model_parameters");
  BOOST_CHECK(resType == "TLineModelResult");

  std::shared_ptr<const TExtremaResult> res =
      Context.GetResultStore()->GetExtremaResult(
          "galaxy", "redshiftSolver", "lineModelSolve", "extrema_results",
          "model_parameters", 0);

  Float64 z = res->Redshift;
  BOOST_CHECK_CLOSE(z, 0.2596216267268967, 1e-6);

  ctx.reset();
}

BOOST_FIXTURE_TEST_CASE(computeTplFitTplRatio_test,
                        fixture_LineModelSolveTestTplFitTplRatio) {
  CLineModelSolve lineModelSolve(Context.m_ScopeStack, "galaxy");
  BOOST_CHECK_NO_THROW(lineModelSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("galaxy", "redshiftSolver",
                                               "lineModelSolve");
  BOOST_CHECK(result_out.lock()->getType() == "CLineModelSolveResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "redshiftSolver", "lineModelSolve", "pdf");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "redshiftSolver", "lineModelSolve", "pdf_params");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  std::string resType = Context.GetResultStore()->GetCandidateResultType(
      "galaxy", "redshiftSolver", "lineModelSolve", "extrema_results",
      "model_parameters");
  BOOST_CHECK(resType == "TLineModelResult");

  std::shared_ptr<const TExtremaResult> res =
      Context.GetResultStore()->GetExtremaResult(
          "galaxy", "redshiftSolver", "lineModelSolve", "extrema_results",
          "model_parameters", 0);

  Float64 z = res->Redshift;
  BOOST_CHECK_CLOSE(z, 0.2596216267268967, 1e-6);

  ctx.reset();
}

BOOST_FIXTURE_TEST_CASE(computeNoContTplRatio_test,
                        fixture_LineModelSolveTestNoContTplRatio) {

  bfs::path logFile = bfs::unique_path("/tmp/log_nocont_tplRatio");
  CLogFileHandler file_handler(logFile.c_str());
  file_handler.SetLevelMask(60);

  CLineModelSolve lineModelSolve(Context.m_ScopeStack, "galaxy");
  BOOST_CHECK_NO_THROW(lineModelSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("galaxy", "redshiftSolver",
                                               "lineModelSolve");
  BOOST_CHECK(result_out.lock()->getType() == "CLineModelSolveResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "redshiftSolver", "lineModelSolve", "pdf");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "redshiftSolver", "lineModelSolve", "pdf_params");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  std::string resType = Context.GetResultStore()->GetCandidateResultType(
      "galaxy", "redshiftSolver", "lineModelSolve", "extrema_results",
      "model_parameters");
  BOOST_CHECK(resType == "TLineModelResult");

  std::shared_ptr<const TExtremaResult> res =
      Context.GetResultStore()->GetExtremaResult(
          "galaxy", "redshiftSolver", "lineModelSolve", "extrema_results",
          "model_parameters", 0);

  Float64 z = res->Redshift;

  ctx.reset();
}

BOOST_FIXTURE_TEST_CASE(computeMultiNoContTplRatio_test,
                        fixture_LineModelSolveTestMultiNoContTplRatio) {
  bfs::path logFile = bfs::unique_path("/tmp/log_multi_nocont_tplRatio");
  CLogFileHandler file_handler(logFile.c_str());
  file_handler.SetLevelMask(60);

  CLineModelSolve lineModelSolve(Context.m_ScopeStack, "galaxy");
  BOOST_CHECK_NO_THROW(lineModelSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("galaxy", "redshiftSolver",
                                               "lineModelSolve");
  BOOST_CHECK(result_out.lock()->getType() == "CLineModelSolveResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "redshiftSolver", "lineModelSolve", "pdf");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "redshiftSolver", "lineModelSolve", "pdf_params");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  std::string resType = Context.GetResultStore()->GetCandidateResultType(
      "galaxy", "redshiftSolver", "lineModelSolve", "extrema_results",
      "model_parameters");
  BOOST_CHECK(resType == "TLineModelResult");

  std::shared_ptr<const TExtremaResult> res =
      Context.GetResultStore()->GetExtremaResult(
          "galaxy", "redshiftSolver", "lineModelSolve", "extrema_results",
          "model_parameters", 0);

  ctx.reset();
}

BOOST_FIXTURE_TEST_CASE(computeFromSpectrum_test,
                        fixture_LineModelSolveTestFromSpectrum) {
  CLineModelSolve lineModelSolve(Context.m_ScopeStack, "galaxy");
  BOOST_CHECK_NO_THROW(lineModelSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("galaxy", "redshiftSolver",
                                               "lineModelSolve");
  BOOST_CHECK(result_out.lock()->getType() == "CLineModelSolveResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "redshiftSolver", "lineModelSolve", "pdf");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "redshiftSolver", "lineModelSolve", "pdf_params");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  std::string resType = Context.GetResultStore()->GetCandidateResultType(
      "galaxy", "redshiftSolver", "lineModelSolve", "extrema_results",
      "model_parameters");
  BOOST_CHECK(resType == "TLineModelResult");

  std::shared_ptr<const TExtremaResult> res =
      Context.GetResultStore()->GetExtremaResult(
          "galaxy", "redshiftSolver", "lineModelSolve", "extrema_results",
          "model_parameters", 0);

  Float64 z = res->Redshift;
  BOOST_CHECK_CLOSE(z, 0.25969245809934272, 1e-6);

  ctx.reset();
}

BOOST_AUTO_TEST_SUITE_END()
