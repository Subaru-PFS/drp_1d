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
    "{\"lambdarange\" : [ 4680, 4712 ],"
    "\"smoothWidth\" : 0.0,"
    "\"templateCatalog\" : {"
    "\"continuumRemoval\" : {"
    "\"method\" : \"zero\","
    "\"medianKernelWidth\" : 75,"
    "\"medianEvenReflection\" : true}},"
    "\"continuumRemoval\" : {"
    "\"method\" : \"IrregularSamplingMedian\","
    "\"medianKernelWidth\" : 400,"
    "\"medianEvenReflection\" : true,"
    "\"decompScales\" : 9},"
    "\"LSF\" : {\"LSFType\" : \"GaussianConstantResolution\", "
    "\"resolution\" : "
    "4300},"
    "\"objects\" : [\"galaxy\"],"
    "\"autocorrectinput\" : false,"
    "\"galaxy\" : {"
    "\"redshiftsampling\" : \"log\","
    "\"method\" : null ,"
    "\"linemeas_method\" : \"LineMeasSolve\","
    "\"linemeas_dzhalf\" : 0.0,"
    "\"linemeas_redshiftstep\" : 0.0001,"
    "\"redshiftref\" : 0.25969245809934272,"
    "\"LineMeasSolve\" : {"
    "\"linemodel\" : {"
    "\"velocityemission\" : 30.0,"
    "\"velocityabsorption\" : 150.0,"
    "\"continuumcomponent\" : \"nocontinuum\","
    "\"linetypefilter\" : \"E\","
    "\"lineforcefilter\" : \"no\","
    "\"nsigmasupport\" : 8,"
    "\"linewidthtype\" : \"combined\","
    "\"fittingmethod\" : \"hybrid\","
    "\"polynomialdegree\" : 2,"
    "\"velocityfit\" : false,"
    "\"lineRatioType\" : \"rules\","
    "\"rules\" : \"no\","
    "\"improveBalmerFit\" : true,"
    "\"lyaforcedisablefit\" : false "
    "}}}}";

const std::string jsonString_lbfgs =
    "{\"lambdarange\" : [ 4680, 4712 ],"
    "\"smoothWidth\" : 0.0,"
    "\"templateCatalog\" : {"
    "\"continuumRemoval\" : {"
    "\"method\" : \"zero\","
    "\"medianKernelWidth\" : 75,"
    "\"medianEvenReflection\" : true}},"
    "\"continuumRemoval\" : {"
    "\"method\" : \"IrregularSamplingMedian\","
    "\"medianKernelWidth\" : 400,"
    "\"medianEvenReflection\" : true,"
    "\"decompScales\" : 9},"
    "\"LSF\" : {\"LSFType\" : \"GaussianConstantResolution\", "
    "\"resolution\" : "
    "4300},"
    "\"objects\" : [\"galaxy\"],"
    "\"autocorrectinput\" : false,"
    "\"galaxy\" : {"
    "\"redshiftsampling\" : \"log\","
    "\"method\" : null ,"
    "\"linemeas_method\" : \"LineMeasSolve\","
    "\"linemeas_dzhalf\" : 0.0,"
    "\"linemeas_redshiftstep\" : 0.0001,"
    "\"redshiftref\" : 0.25969245809934272,"
    "\"LineMeasSolve\" : {"
    "\"linemodel\" : {"
    "\"velocityemission\" : 30.0,"
    "\"velocityabsorption\" : 150.0,"
    "\"continuumcomponent\" : \"nocontinuum\","
    "\"linetypefilter\" : \"E\","
    "\"lineforcefilter\" : \"no\","
    "\"nsigmasupport\" : 14,"
    "\"linewidthtype\" : \"combined\","
    "\"fittingmethod\" : \"lbfgs\","
    "\"polynomialdegree\" : 2,"
    "\"velocityfit\" : true,"
    "\"emvelocityfitmin\" : 10,"
    "\"emvelocityfitmax\" : 400,"
    "\"absvelocityfitmin\" : 150,"
    "\"absvelocityfitmax\" : 500,"
    "\"lineRatioType\" : \"rules\","
    "\"rules\" : \"no\","
    "\"improveBalmerFit\" : true,"
    "\"lyaforcedisablefit\" : false "
    "}}}}";

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
    ctx.setLineCatalog("galaxy", "LineMeasSolve", lineCatalog);
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

#ifdef LBFGSFITTER
class fixture_LinemeasSolveLbfgsTest {
public:
  fixture_Context ctx;

  fixture_LinemeasSolveLbfgsTest() {
    fillCatalog();
    ctx.loadParameterStore(jsonString_lbfgs);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.setPhotoBandCatalog(photoBandCatalog);
    spc->SetPhotData(photoData);
    ctx.addSpectrum(spc, LSF);
    ctx.setLineRatioCatalogCatalog("galaxy", lineRatioTplCatalog);
    ctx.setLineCatalog("galaxy", "LineMeasSolve", lineCatalog);
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
      Context.GetResultStore()->GetSolveResult("galaxy", "LineMeasSolve");
  BOOST_CHECK(result_out.lock()->getType() == "CLineMeasSolveResult");

  std::shared_ptr<const CLineModelSolution> res =
      Context.GetResultStore()->GetLineModelSolution("galaxy", "LineMeasSolve",
                                                     "linemeas");

  Float64 snrOII = res->snrOII;
  BOOST_CHECK_CLOSE(snrOII, 28.099083565445135, 1e-6);

  Float64 lfOII = res->lfOII;
  BOOST_CHECK_CLOSE(lfOII, -15.658705018780491, 1e-6);

  ctx.reset();
}

#ifdef LBFGSFITTER
BOOST_FIXTURE_TEST_CASE(compute_test_lbfgs, fixture_LinemeasSolveLbfgsTest) {

  CLineMeasSolve lineMeasSolve(Context.m_ScopeStack, "galaxy");
  BOOST_CHECK_NO_THROW(lineMeasSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("galaxy", "LineMeasSolve");
  BOOST_CHECK(result_out.lock()->getType() == "CLineMeasSolveResult");

  std::shared_ptr<const CLineModelSolution> res =
      Context.GetResultStore()->GetLineModelSolution("galaxy", "LineMeasSolve",
                                                     "linemeas");

  Float64 snrOII = res->snrOII;
  BOOST_CHECK_CLOSE(snrOII, 20.914260055007116, 1e-6);

  Float64 lfOII = res->lfOII;
  BOOST_CHECK_CLOSE(lfOII, -15.786954670784215, 1e-6);

  ctx.reset();
}
#endif

BOOST_AUTO_TEST_SUITE_END()