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
#include "RedshiftLibrary/method/linemodelsolve.h"
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
    "\"LSF\" : {\"LSFType\" : \"GaussianConstantResolution\", "
    "\"resolution\" : "
    "4300},"
    "\"extremaredshiftseparation\" : 0.01,"
    "\"objects\" : [\"galaxy\"],"
    "\"autocorrectinput\" : false,"
    "\"airvacuum_method\" : \"default\","
    "\"galaxy\" : {"
    "\"redshiftrange\" : [ 0.24, 0.3 ],"
    "\"redshiftstep\" : 0.0001,"
    "\"redshiftsampling\" : \"log\","
    "\"method\" : \"LineModelSolve\","
    "\"LineModelSolve\" : {"
    "\"linemodel\" : {"
    "\"continuumreestimation\" : \"no\","
    "\"velocityfit\" : true,"
    "\"emvelocityfitmin\" : 100,"
    "\"emvelocityfitmax\" : 700, "
    "\"emvelocityfitstep\" : 20,"
    "\"absvelocityfitmin\" : 150,"
    "\"absvelocityfitmax\" : 500, "
    "\"absvelocityfitstep\" : 50,"
    "\"extremacount\" : 5,"
    "\"extremacountB\" : 3,"
    "\"nsigmasupport\" : 8,"
    "\"haprior\" : 0.5,"
    "\"euclidnhaemittersStrength\" : 1.0,"
    "\"extremacutprobathreshold\" : -1,"
    "\"skipsecondpass\" : false,"
    "\"secondpasslcfittingmethod\" : -1,"
    "\"useloglambdasampling\": false,"
    "\"lyaforcefit\": false,"
    "\"lyaforcedisablefit\": false,"
    "\"stronglinesprior\" : 1.0,"
    "\"fittingmethod\": \"svd\","
    "\"linewidthtype\": \"combined\","
    "\"velocityemission\" : 100,"
    "\"velocityabsorption\": 100,"
    "\"linetypefilter\" : \"no\","
    "\"lineforcefilter\" : \"no\","
    "\"lyafit\": {"
    "\"asymfitmin\" : 0,"
    "\"asymfitmax\" : 4, \"asymfitstep\" : 1, "
    "\"widthfitmin\" : 1,"
    "\"widthfitmax\" : 4, \"widthfitstep\" : 1, "
    "\"deltafitmin\" : 0,"
    "\"deltafitmax\" : 0, \"deltafitstep\" : 1}, "
    "\"tplratio\": { \"priors\": {"
    "\"betaA\" : 1,    \"betaTE\" : 1, \"betaZ\" : 1, "
    "\"catalog_dirpath\" : \"\"}}, "
    "\"firstpass\": { \"fittingmethod\" : \"individual\", "
    "\"tplratio_ismfit\" : true,"
    "\"largegridstepratio\" : 6, "
    "\"multiplecontinuumfit_disable\": true},"
    "\"secondpass\" : {\"halfwindowsize\" : 0.001, "
    "\"continuumfit\" : \"refitfirstpass\"},";

const std::string jsonStringTplFitRules =
    "\"continuumcomponent\" : \"tplfit\","
    "\"pdfcombination\" : \"marg\","
    "\"tplratio_ismfit\" : true,"
    "\"rules\" : \"all\","
    "\"improveBalmerFit\" : true,"
    "\"lineRatioType\": \"rules\","
    "\"enablephotometry\" : true, "
    "\"photometry\" : {\"weight\" : 1},"
    "\"continuumfit\" : { \"ignorelinesupport\": false,"
    "\"negativethreshold\": -5.0,"
    "\"count\" : 1,"
    "\"nullthreshold\": 3,"
    "\"ismfit\" : true,"
    "\"igmfit\" : true,"
    "\"fftprocessing\": false, "
    "\"priors\": { \"betaA\" : 1, \"betaTE\" : 1, \"betaZ\" : 1,"
    "\"catalog_dirpath\" : \"\"}}}}}}";

const std::string jsonStringTplFitTplRatio =
    "\"continuumcomponent\" : \"tplfit\","
    "\"pdfcombination\" : \"marg\","
    "\"tplratio_ismfit\" : true,"
    "\"rules\" : \"all\","
    "\"improveBalmerFit\" : true,"
    "\"lineRatioType\": \"tplratio\","
    "\"enablephotometry\" : true, "
    "\"photometry\" : {\"weight\" : 1},"
    "\"continuumfit\" : { \"ignorelinesupport\": false,"
    "\"negativethreshold\": -5.0,"
    "\"count\" : 1,"
    "\"nullthreshold\": 3,"
    "\"ismfit\" : true,"
    "\"igmfit\" : true,"
    "\"fftprocessing\": false, "
    "\"priors\": { \"betaA\" : 1, \"betaTE\" : 1, \"betaZ\" : 1,"
    "\"catalog_dirpath\" : \"\"}}}}}}";

const std::string jsonStringFromSpectrum =
    "\"continuumcomponent\" : \"fromspectrum\","
    "\"pdfcombination\" : \"bestchi2\","
    "\"tplratio_ismfit\" : true,"
    "\"rules\" : \"balmerSingle\","
    "\"improveBalmerFit\" : true,"
    "\"lineRatioType\": \"tplratio\","
    "\"continuumfit\" : { \"ignorelinesupport\": false,"
    "\"negativethreshold\": -5.0,"
    "\"count\" : 1,"
    "\"nullthreshold\": 3,"
    "\"ismfit\" : true,"
    "\"igmfit\" : true,"
    "\"fftprocessing\": true, "
    "\"priors\": { \"betaA\" : 1, \"betaTE\" : 1, \"betaZ\" : 1,"
    "\"catalog_dirpath\" : \"\"}}}}}}";

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
    ctx.loadParameterStore(jsonString + jsonStringTplFitRules);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.setPhotoBandCatalog(photoBandCatalog);
    spc->SetPhotData(photoData);
    ctx.addSpectrum(spc, LSF);
    ctx.setLineRatioCatalogCatalog("galaxy", lineRatioTplCatalog);
    ctx.setLineCatalog("galaxy", "LineModelSolve", lineCatalog);
    ctx.initContext();
    lineRatioTplCatalog->addLineRatioCatalog(*lineRatioCatalog);
  }
};

class fixture_LineModelSolveTestTplFitTplRatio
    : public fixture_LineModelSolveTest {
public:
  fixture_LineModelSolveTestTplFitTplRatio() {
    fillCatalog();
    ctx.loadParameterStore(jsonString + jsonStringTplFitTplRatio);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.setPhotoBandCatalog(photoBandCatalog);
    spc->SetPhotData(photoData);
    ctx.addSpectrum(spc, LSF);
    ctx.setLineRatioCatalogCatalog("galaxy", lineRatioTplCatalog);
    ctx.setLineCatalog("galaxy", "LineModelSolve", lineCatalog);
    ctx.initContext();
    lineRatioTplCatalog->addLineRatioCatalog(*lineRatioCatalog);
  }
};

class fixture_LineModelSolveTestFromSpectrum
    : public fixture_LineModelSolveTest {
public:
  fixture_LineModelSolveTestFromSpectrum() {
    fillCatalog();
    ctx.loadParameterStore(jsonString + jsonStringFromSpectrum);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.setPhotoBandCatalog(photoBandCatalog);
    spc->SetPhotData(photoData);
    ctx.addSpectrum(spc, LSF);
    ctx.setLineRatioCatalogCatalog("galaxy", lineRatioTplCatalog);
    ctx.setLineCatalog("galaxy", "LineModelSolve", lineCatalog);
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
      Context.GetResultStore()->GetSolveResult("galaxy", "LineModelSolve");
  BOOST_CHECK(result_out.lock()->getType() == "CLineModelSolveResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "LineModelSolve", "pdf");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "LineModelSolve", "pdf_params");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  std::string resType = Context.GetResultStore()->GetCandidateResultType(
      "galaxy", "LineModelSolve", "extrema_results", "model_parameters");
  BOOST_CHECK(resType == "TLineModelResult");

  std::shared_ptr<const TExtremaResult> res =
      Context.GetResultStore()->GetExtremaResult(
          "galaxy", "LineModelSolve", "extrema_results", "model_parameters", 0);

  Float64 z = res->Redshift;
  BOOST_CHECK_CLOSE(z, 0.25886608043741388, 1e-6);

  ctx.reset();
}

BOOST_FIXTURE_TEST_CASE(computeTplFitTplRatio_test,
                        fixture_LineModelSolveTestTplFitTplRatio) {
  CLineModelSolve lineModelSolve(Context.m_ScopeStack, "galaxy");
  BOOST_CHECK_NO_THROW(lineModelSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("galaxy", "LineModelSolve");
  BOOST_CHECK(result_out.lock()->getType() == "CLineModelSolveResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "LineModelSolve", "pdf");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "LineModelSolve", "pdf_params");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  std::string resType = Context.GetResultStore()->GetCandidateResultType(
      "galaxy", "LineModelSolve", "extrema_results", "model_parameters");
  BOOST_CHECK(resType == "TLineModelResult");

  std::shared_ptr<const TExtremaResult> res =
      Context.GetResultStore()->GetExtremaResult(
          "galaxy", "LineModelSolve", "extrema_results", "model_parameters", 0);

  Float64 z = res->Redshift;
  BOOST_CHECK_CLOSE(z, 0.2591178788325017, 1e-6);

  ctx.reset();
}

BOOST_FIXTURE_TEST_CASE(computeFromSpectrum_test,
                        fixture_LineModelSolveTestFromSpectrum) {
  CLineModelSolve lineModelSolve(Context.m_ScopeStack, "galaxy");
  BOOST_CHECK_NO_THROW(lineModelSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("galaxy", "LineModelSolve");
  BOOST_CHECK(result_out.lock()->getType() == "CLineModelSolveResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "LineModelSolve", "pdf");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "LineModelSolve", "pdf_params");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  std::string resType = Context.GetResultStore()->GetCandidateResultType(
      "galaxy", "LineModelSolve", "extrema_results", "model_parameters");
  BOOST_CHECK(resType == "TLineModelResult");

  std::shared_ptr<const TExtremaResult> res =
      Context.GetResultStore()->GetExtremaResult(
          "galaxy", "LineModelSolve", "extrema_results", "model_parameters", 0);

  Float64 z = res->Redshift;
  BOOST_CHECK_CLOSE(z, 0.25905236255797504, 1e-6);

  ctx.reset();
}

BOOST_AUTO_TEST_SUITE_END()