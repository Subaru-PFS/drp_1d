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
#include "RedshiftLibrary/method/linemeassolve.h"
#include "RedshiftLibrary/processflow/context.h"
#include "tests/src/tool/inputContextLight.h"
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

const std::string jsonString =
    "{\"lambdaRange\" : [ 4680, 4712 ],"
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
    "\"lsf\" : {\"lsfType\" : \"gaussianConstantResolution\", "
    "\"resolution\" : "
    "4300},"
    "\"spectrumModels\" : [\"galaxy\"],"
    "\"autoCorrectInput\" : false,"
    "\"galaxy\" : {"
    "\"redshiftSampling\" : \"log\","
    "\"method\" : null ,"
    "\"linemeas_method\" : \"lineMeasSolver\","
    "\"lineMeasDzHalf\" : 0.0,"
    "\"lineMeasRedshiftStep\" : 0.0001,"
    "\"redshiftref\" : 0.25969245809934272,"
    "\"lineMeasSolver\" : {"
    "\"lineModel\" : {"
    "\"velocityEmission\" : 30.0,"
    "\"velocityAbsorption\" : 150.0,"
    "\"continuumComponent\" : \"noContinuum\","
    "\"lineTypeFilter\" : \"E\","
    "\"lineForceFilter\" : \"no\","
    "\"nSigmaSupport\" : 8,"
    "\"lineWidthType\" : \"combined\","
    "\"fittingMethod\" : \"hybrid\","
    "\"polynomialDegree\" : 2,"
    "\"velocityFit\" : false,"
    "\"ampOffsetFit\": \"true\","
    "\"lbdaOffsetFit\": \"true\","
    "\"lineRatioType\" : \"rules\","
    "\"rules\" : \"no\","
    "\"improveBalmerFit\" : true,"
    "\"lyaForceDisableFit\" : false "
    "}}}}";

// const std::string jsonString_lbfgsb =
//     "{\"lambdaRange\" : [ 4680, 4712 ],"
//     "\"smoothWidth\" : 0.0,"
//     "\"templateCatalog\" : {"
//     "\"continuumRemoval\" : {"
//     "\"method\" : \"zero\","
//     "\"medianKernelWidth\" : 75,"
//     "\"medianEvenReflection\" : true}},"
//     "\"continuumRemoval\" : {"
//     "\"method\" : \"irregularSamplingMedian\","
//     "\"medianKernelWidth\" : 400,"
//     "\"medianEvenReflection\" : true,"
//     "\"decompScales\" : 9},"
//     "\"lsf\" : {\"lsfType\" : \"gaussianConstantResolution\", "
//     "\"resolution\" : "
//     "4300},"
//     "\"spectrumModels\" : [\"galaxy\"],"
//     "\"autoCorrectInput\" : false,"
//     "\"galaxy\" : {"
//     "\"redshiftSampling\" : \"log\","
//     "\"method\" : null ,"
//     "\"linemeas_method\" : \"lineMeasSolver\","
//     "\"lineMeasDzHalf\" : 0.0,"
//     "\"lineMeasRedshiftStep\" : 0.0001,"
//     "\"redshiftref\" : 0.25969245809934272,"
//     "\"lineMeasSolver\" : {"
//     "\"lineModel\" : {"
//     "\"velocityEmission\" : 30.0,"
//     "\"velocityAbsorption\" : 150.0,"
//     "\"continuumComponent\" : \"noContinuum\","
//     "\"lineTypeFilter\" : \"E\","
//     "\"lineForceFilter\" : \"no\","
//     "\"nSigmaSupport\" : 14,"
//     "\"lineWidthType\" : \"combined\","
//     "\"fittingMethod\" : \"lbfgsb\","
//     "\"polynomialDegree\" : 2,"
//     "\"velocityFit\" : true,"
//     "\"ampOffsetFit\": \"true\","
//     "\"lbdaOffsetFit\": \"true\","
//     "\"emVelocityFitMin\" : 10,"
//     "\"emVelocityFitMax\" : 400,"
//     "\"absVelocityFitMin\" : 150,"
//     "\"absVelocityFitMax\" : 500,"
//     "\"lineRatioType\" : \"rules\","
//     "\"rules\" : \"no\","
//     "\"improveBalmerFit\" : true,"
//     "\"lyaForceDisableFit\" : false "
//     "}}}}";

const std::string jsonString_lbfgsb = "{"
"\"lambdaRange\": ["
"4680,"
"4712"
"],"
"\"smoothWidth\": 0.0,"
"\"templateCatalog\": {"
"\"continuumRemoval\": {"
"\"method\": \"zero\","
"\"medianKernelWidth\": 75,"
"\"medianEvenReflection\": true"
"}"
"},"
"\"continuumRemoval\": {"
"\"method\": \"irregularSamplingMedian\","
"\"medianKernelWidth\": 400,"
"\"medianEvenReflection\": true"
"},"
"\"lsf\": {"
"\"lsfType\": \"gaussianConstantResolution\","
"\"resolution\": 4300"
"},"
"\"spectrumModels\": ["
"\"galaxy\""
"],"
"\"autoCorrectInput\": false,"
"\"spectrumModel_galaxy\": {"
"\"stages\": null,"
"\"redshiftSampling\": \"log\","
"\"lineMeasDzHalf\": 0.0,"
"\"lineMeasRedshiftStep\": 0.0001,"
"\"lineMeasSolver\": {"
"\"lineModel\": {"
"\"nSigmaSupport\": 8,"
"\"lineWidthType\": \"combined\","
"\"fittingMethod\": \"hybrid\","
"\"polynomialDegree\": 2,"
"\"velocityFit\": false,"
"\"ampOffsetFit\": \"true\","
"\"lbdaOffsetFit\": \"true\","
"\"emVelocityFitMin\": 10,"
"\"emVelocityFitMax\": 400,"
"\"absVelocityFitMin\": 150,"
"\"absVelocityFitMax\": 500,"
"\"lineRatioType\": \"rules\","
"\"rules\": \"no\""
"}"
"}"
"}"
"}";

class fixture_LinemeasSolveTest {
public:
  fixture_Context ctx;

  fixture_LinemeasSolveTest() {
    fillCatalog();
    ctx.loadParameterStore(jsonString);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.setPhotoBandCatalog(photoBandCatalog);
    spc->SetPhotData(photoData);
    ctx.addSpectrum(spc, LSF);
    ctx.setLineRatioCatalogCatalog("galaxy", lineRatioTplCatalog);
    ctx.setLineCatalog("galaxy", "lineMeasSolver", lineCatalog);
    ctx.initContext();
    lineRatioTplCatalog->addLineRatioCatalog(*lineRatioCatalog);
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

#ifdef LBFGSBFITTER
class fixture_LinemeasSolveLbfgsbTest {
public:
  fixture_Context ctx;

  fixture_LinemeasSolveLbfgsbTest() {
    fillCatalog();
    ctx.loadParameterStore(jsonString_lbfgsb);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.setPhotoBandCatalog(photoBandCatalog);
    spc->SetPhotData(photoData);
    ctx.addSpectrum(spc, LSF);
    ctx.setLineRatioCatalogCatalog("galaxy", lineRatioTplCatalog);
    ctx.setLineCatalog("galaxy", "lineMeasSolver", lineCatalog);
    ctx.initContext();
    lineRatioTplCatalog->addLineRatioCatalog(*lineRatioCatalog);
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
#endif

BOOST_AUTO_TEST_SUITE(linemeasSolve_test)

BOOST_FIXTURE_TEST_CASE(compute_test, fixture_LinemeasSolveTest) {
  CLineMeasSolve lineMeasSolve(Context.m_ScopeStack, "galaxy");
  BOOST_CHECK_NO_THROW(lineMeasSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("galaxy", "lineMeasSolver");
  BOOST_CHECK(result_out.lock()->getType() == "CLineMeasSolveResult");

  std::shared_ptr<const CLineModelSolution> res =
      Context.GetResultStore()->GetLineModelSolution("galaxy", "lineMeasSolver",
                                                     "linemeas");

  Float64 snrOII = res->snrOII;
  BOOST_CHECK_CLOSE(snrOII, 16.196486940733053, 1e-6);

  Float64 lfOII = res->lfOII;
  BOOST_CHECK_CLOSE(lfOII, -15.658215485050579, 1e-6);

  ctx.reset();
}

#ifdef LBFGSBFITTER
BOOST_FIXTURE_TEST_CASE(compute_test_lbfgs, fixture_LinemeasSolveLbfgsbTest) {

  CLineMeasSolve lineMeasSolve(Context.m_ScopeStack, "galaxy");
  BOOST_CHECK_NO_THROW(lineMeasSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("galaxy", "lineMeasSolver");
  BOOST_CHECK(result_out.lock()->getType() == "CLineMeasSolveResult");

  std::shared_ptr<const CLineModelSolution> res =
      Context.GetResultStore()->GetLineModelSolution("galaxy", "lineMeasSolver",
                                                     "linemeas");

  Float64 snrOII = res->snrOII_DI;
  BOOST_CHECK_CLOSE(snrOII, 21.480993641608535, 1e-2);

  Float64 lfOII = res->lfOII;
  BOOST_CHECK_CLOSE(lfOII, -15.778872598441525, 1e-4);

  ctx.reset();
}
#endif

BOOST_AUTO_TEST_SUITE_END()
