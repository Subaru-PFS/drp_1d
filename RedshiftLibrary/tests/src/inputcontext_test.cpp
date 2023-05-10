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
#include "RedshiftLibrary/processflow/inputcontext.h"
#include "tests/src/tool/inputContextLight.h"
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

const std::string jsonStringOneSpc = "{\"lambdarange\" : [ 4631, 4815 ],";

const std::string jsonStringMO =
    "{\"lambdarange\" : {\"1\" : [ 4631, 4815], \"2\" : [ 4631, 4815 ]},";

const std::string jsonString =
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
    "\"autocorrectinput\" : false,"
    "\"airvacuum_method\" : \"default\",";

const std::string jsonStringNoFFT =
    "\"galaxy\" : {"
    "\"redshiftrange\" : [ 2.84, 2.88 ],"
    "\"redshiftstep\" : 0.0001,"
    "\"redshiftsampling\" : \"log\","
    "\"method\" : \"TemplateFittingSolve\","
    "\"TemplateFittingSolve\" : {"
    "\"extremacount\" : 5,"
    "\"overlapThreshold\" : 1,"
    "\"spectrum\" : {\"component\" : \"raw\"},"
    "\"fftprocessing\" : false,"
    "\"interpolation\" : \"precomputedfinegrid\","
    "\"extinction\" : true,"
    "\"dustfit\" : true,"
    "\"pdfcombination\" : \"marg\","
    "\"enablephotometry\" : false}}}";

const std::string jsonStringFFT = "\"galaxy\" : {"
                                  "\"redshiftrange\" : [ 2.84, 2.88 ],"
                                  "\"redshiftstep\" : 0.0001,"
                                  "\"redshiftsampling\" : \"log\","
                                  "\"method\" : \"TemplateFittingSolve\","
                                  "\"TemplateFittingSolve\" : {"
                                  "\"extremacount\" : 5,"
                                  "\"overlapThreshold\" : 1,"
                                  "\"spectrum\" : {\"component\" : \"raw\"},"
                                  "\"fftprocessing\" : true,"
                                  "\"interpolation\" : \"precomputedfinegrid\","
                                  "\"extinction\" : true,"
                                  "\"dustfit\" : true,"
                                  "\"pdfcombination\" : \"marg\","
                                  "\"enablephotometry\" : false}}}";

const std::string jsonStringOrtho =
    "\"galaxy\" : {"
    "\"redshiftrange\" : [ 2.84, 2.88 ],"
    "\"redshiftstep\" : 0.0001,"
    "\"redshiftsampling\" : \"log\","
    "\"method\" : \"LineModelSolve\","
    "\"LineModelSolve\" : {"
    "\"linemodel\" : {"
    "\"continuumcomponent\" : \"tplfit\","
    "\"useloglambdasampling\": false,"
    "\"lyaforcefit\": false,"
    "\"lyaforcedisablefit\": false,"
    "\"lyafit\": {"
    "\"asymfitmin\" : 0,"
    "\"asymfitmax\" : 4, \"asymfitstep\" : 1, \"widthfitmin\" : 1,"
    "\"widthfitmax\" : 4, \"widthfitstep\" : 1, \"deltafitmin\" : 0,"
    "\"deltafitmax\" : 0, \"deltafitstep\" : 1}, "
    "\"continuumfit\" : { \"ignorelinesupport\": false,"
    "\"negativethreshold\": -5.0,"
    "\"nullthreshold\": 3,"
    "\"count\":1,"
    "\"fftprocessing\": false} ,"
    "\"firstpass\": { \"fittingmethod\" : \"individual\", "
    "\"multiplecontinuumfit_disable\": true},"
    "\"fittingmethod\": \"individual\","
    "\"ampoffsetfit\": \"false\","
    "\"lbdaoffsetfit\": \"false\","
    "\"linewidthtype\": \"combined\","
    "\"velocityemission\" : 100,"
    "\"velocityabsorption\": 100,"
    "\"lineRatioType\": \"tplratio\","
    "\"linetypefilter\" : \"no\","
    "\"lineforcefilter\" : \"no\"}}}}";

class fixture_inputcontextTest {
public:
  TScopeStack scopeStack;
  std::shared_ptr<CParameterStore> paramStoreNoFFT =
      fixture_ParamStore(jsonStringOneSpc + jsonString + jsonStringNoFFT,
                         scopeStack)
          .paramStore;
  std::shared_ptr<CParameterStore> paramStoreFFT =
      fixture_ParamStore(jsonStringOneSpc + jsonString + jsonStringFFT,
                         scopeStack)
          .paramStore;
  std::shared_ptr<CParameterStore> paramStoreMO =
      fixture_ParamStore(jsonStringMO + jsonString + jsonStringFFT, scopeStack)
          .paramStore;
  std::shared_ptr<CParameterStore> paramStoreOrtho =
      fixture_ParamStore(jsonStringOneSpc + jsonString + jsonStringOrtho,
                         scopeStack)
          .paramStore;
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin =
      fixture_MeiskinCorrection().igmCorrectionMeiksin;
  std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti =
      fixture_CalzettiCorrection().ismCorrectionCalzetti;
  std::shared_ptr<CLSF> LSF =
      fixture_LSFGaussianConstantResolution(scopeStack).LSF;
  std::shared_ptr<CSpectrum> spc = fixture_SharedSpectrumExtended().spc;
  std::shared_ptr<CSpectrum> spc2 = fixture_SharedSpectrumExtended().spc;
  std::shared_ptr<CTemplateCatalog> catalog =
      fixture_sharedTemplateCatalog().catalog;
  std::shared_ptr<CPhotBandCatalog> photoBandCatalog =
      fixture_PhotoBandCatalog().photoBandCatalog;
  std::shared_ptr<CLineCatalog> lineCatalog = fixture_LineCatalog().lineCatalog;
  std::shared_ptr<CLineCatalogsTplRatio> lineRatioTplCatalog =
      fixture_LineRatioTplCatalog().lineRatioTplCatalog;

  // list
  TFloat64List linLambdaList = fixture_SpectralAxis().linLambdaList;

  void setInputData(CInputContext &inputCtx) {
    spc->SetLSF(LSF);
    inputCtx.addSpectrum(spc);
    catalog->Add(fixture_SharedGalaxyTemplate().tpl);
    inputCtx.setTemplateCatalog(catalog);
    inputCtx.setPhotBandCatalog(photoBandCatalog);
    inputCtx.setLineRatioCatalogCatalog("galaxy", lineRatioTplCatalog);
    scopeStack.push_back("galaxy");
    inputCtx.setLineCatalog("galaxy", "LineModelSolve", lineCatalog);
    inputCtx.setfluxCorrectionCalzetti(ismCorrectionCalzetti);
    inputCtx.setfluxCorrectionMeiksin(igmCorrectionMeiksin);
  }

  void setInputDataMO(CInputContext &inputCtx) {
    setInputData(inputCtx);
    spc->setObsID("1");
    inputCtx.addSpectrum(spc2);
    spc2->setObsID("2");
  }
};

BOOST_FIXTURE_TEST_SUITE(inputContext_test, fixture_inputcontextTest)

BOOST_AUTO_TEST_CASE(getterSetter_test) {
  CInputContext inputCtx(paramStoreNoFFT);
  std::shared_ptr<CInputContext> inputCtxPtr =
      std::make_shared<CInputContext>(paramStoreNoFFT);

  spc->SetLSF(LSF);
  inputCtx.addSpectrum(spc);
  BOOST_CHECK(inputCtx.GetSpectrum() == spc);

  BOOST_CHECK(inputCtx.getRebinnedSpectra().size() == 0);

  inputCtx.setTemplateCatalog(catalog);
  BOOST_CHECK(inputCtx.GetTemplateCatalog() == catalog);

  inputCtx.setPhotBandCatalog(photoBandCatalog);
  BOOST_CHECK(inputCtx.GetPhotBandCatalog() == photoBandCatalog);

  inputCtx.setLineRatioCatalogCatalog("galaxy", lineRatioTplCatalog);
  scopeStack.push_back("galaxy");
  BOOST_CHECK(inputCtx.GetTemplateRatioCatalog(scopeStack[0]) ==
              lineRatioTplCatalog);

  inputCtx.setLineCatalog("galaxy", "LineModelSolve", lineCatalog);
  BOOST_CHECK(inputCtx.GetLineCatalog("galaxy", "LineModelSolve") ==
              lineCatalog);

  BOOST_CHECK(
      inputCtx.GetFilteredLineVector("galaxy", "LineModelSolve", "no", "no")
          .size() == fixture_LineCatalog().lineCatalogSize);

  BOOST_CHECK(inputCtx.m_ismcorrectionCalzetti == nullptr);
  inputCtx.setfluxCorrectionCalzetti(ismCorrectionCalzetti);
  BOOST_CHECK(inputCtx.m_ismcorrectionCalzetti == ismCorrectionCalzetti);

  BOOST_CHECK(inputCtx.m_igmcorrectionMeiksin == nullptr);
  inputCtx.setfluxCorrectionMeiksin(igmCorrectionMeiksin);
  BOOST_CHECK(inputCtx.m_igmcorrectionMeiksin == igmCorrectionMeiksin);

  // lambda range
  BOOST_CHECK(inputCtx.m_lambdaRanges.size() == 0);
  inputCtx.m_lambdaRanges.push_back(std::make_shared<TFloat64Range>(1, 3));
  BOOST_CHECK(inputCtx.getLambdaRange()->GetBegin() == 1);
  BOOST_CHECK(inputCtx.getLambdaRange()->GetEnd() == 3);

  BOOST_CHECK(inputCtx.m_clampedLambdaRanges.size() == 0);
  inputCtx.m_clampedLambdaRanges.push_back(
      std::make_shared<TFloat64Range>(5, 6));
  BOOST_CHECK(inputCtx.getClampedLambdaRange(0)->GetBegin() == 5);
  BOOST_CHECK(inputCtx.getClampedLambdaRange(0)->GetEnd() == 6);

  BOOST_CHECK(inputCtx.m_rebinnedClampedLambdaRanges.size() == 0);
  inputCtx.m_rebinnedClampedLambdaRanges.push_back(
      std::make_shared<TFloat64Range>(15, 20));
  BOOST_CHECK(inputCtx.getClampedLambdaRange(1)->GetBegin() == 15);
  BOOST_CHECK(inputCtx.getClampedLambdaRange(1)->GetEnd() == 20);

  inputCtx.m_constClampedLambdaRanges.push_back(
      std::make_shared<TFloat64Range>(5, 6));
  BOOST_CHECK(inputCtx.getClampedLambdaRanges().size() == 1);
  BOOST_CHECK(inputCtx.getClampedLambdaRanges()[0]->GetBegin() ==
              inputCtx.getClampedLambdaRange(0)->GetBegin());
  BOOST_CHECK(inputCtx.getClampedLambdaRanges()[0]->GetEnd() ==
              inputCtx.getClampedLambdaRange(0)->GetEnd());

  inputCtx.m_constRebinnedClampedLambdaRanges.push_back(
      std::make_shared<TFloat64Range>(15, 20));
  BOOST_CHECK(inputCtx.getRebinnedClampedLambdaRanges().size() == 1);
  BOOST_CHECK(inputCtx.getRebinnedClampedLambdaRanges()[0]->GetBegin() ==
              inputCtx.getClampedLambdaRange(1)->GetBegin());
  BOOST_CHECK(inputCtx.getRebinnedClampedLambdaRanges()[0]->GetEnd() ==
              inputCtx.getClampedLambdaRange(1)->GetEnd());
}

BOOST_AUTO_TEST_CASE(initAndReset_test) {
  CInputContext inputCtx(paramStoreNoFFT);
  setInputData(inputCtx);
  inputCtx.Init();
  BOOST_CHECK(inputCtx.GetSpectrum(0) != nullptr);
  BOOST_CHECK(inputCtx.GetTemplateCatalog() != nullptr);
  BOOST_CHECK(inputCtx.GetLineCatalog("galaxy", "LineModelSolve") != nullptr);
  BOOST_CHECK(inputCtx.GetTemplateRatioCatalog(scopeStack[0]) != nullptr);
  BOOST_CHECK(inputCtx.GetPhotBandCatalog() != nullptr);

  inputCtx.resetSpectrumSpecific();

  BOOST_CHECK(inputCtx.getSpectra().size() == 0);
  BOOST_CHECK(inputCtx.GetTemplateCatalog() == nullptr);
  BOOST_CHECK(inputCtx.GetLineCatalog("galaxy", "LineModelSolve") == nullptr);
  BOOST_CHECK(inputCtx.GetTemplateRatioCatalog(scopeStack[0]) == nullptr);
  BOOST_CHECK(inputCtx.GetPhotBandCatalog() == nullptr);

  // rebined spectrum (FFT)
  CInputContext inputCtx2(paramStoreFFT);
  setInputData(inputCtx2);
  inputCtx2.Init();
  BOOST_CHECK(inputCtx2.GetRebinnedSpectrum(0) != nullptr);
  inputCtx2.resetSpectrumSpecific();
  BOOST_CHECK(inputCtx2.getRebinnedSpectra().size() == 0);

  // multi-obs
  CInputContext inputCtx3(paramStoreMO);
  setInputDataMO(inputCtx3);
  BOOST_CHECK(inputCtx3.GetSpectrum(0, 0)->getObsID() == "1");
  BOOST_CHECK(inputCtx3.GetSpectrum(0, 1)->getObsID() == "2");
  inputCtx3.Init();
  BOOST_CHECK(inputCtx3.GetRebinnedSpectrum(0) != nullptr);
  BOOST_CHECK(inputCtx3.GetRebinnedSpectrum(1) != nullptr);
  inputCtx3.resetSpectrumSpecific();
  BOOST_CHECK(inputCtx2.getRebinnedSpectra().size() == 0);
}

BOOST_AUTO_TEST_CASE(rebinInputs_test) {
  CInputContext inputCtx(paramStoreNoFFT);
  setInputData(inputCtx);
  inputCtx.m_categories = paramStoreNoFFT->GetList<std::string>("objects");
  inputCtx.m_lambdaRanges.push_back(std::make_shared<TFloat64Range>(
      paramStoreNoFFT->Get<TFloat64Range>("lambdarange")));

  inputCtx.RebinInputs();
  BOOST_CHECK(inputCtx.m_use_LogLambaSpectrum == false);
  BOOST_CHECK(inputCtx.getRebinnedSpectra().size() == 0);

  CInputContext inputCtx2(paramStoreFFT);
  setInputData(inputCtx2);
  inputCtx2.m_categories = paramStoreFFT->GetList<std::string>("objects");
  inputCtx2.m_lambdaRanges.push_back(std::make_shared<TFloat64Range>(
      paramStoreFFT->Get<TFloat64Range>("lambdarange")));
  inputCtx2.m_rebinnedClampedLambdaRanges.push_back(
      std::make_shared<TFloat64Range>(
          paramStoreFFT->Get<TFloat64Range>("lambdarange")));

  inputCtx2.RebinInputs();
  BOOST_CHECK(inputCtx2.m_use_LogLambaSpectrum == true);
  BOOST_CHECK(inputCtx2.GetRebinnedSpectrum() != nullptr);

  inputCtx2.m_lambdaRanges[0]->Set(1000, 2000);
  BOOST_CHECK_THROW(inputCtx2.RebinInputs(), GlobalException);

  inputCtx2.resetSpectrumSpecific();
  spc->SetSpectralAxis(CSpectrumSpectralAxis(linLambdaList));
  setInputData(inputCtx2);
  inputCtx2.m_lambdaRanges.push_back(std::make_shared<TFloat64Range>(
      paramStoreFFT->Get<TFloat64Range>("lambdarange")));
  inputCtx2.m_rebinnedClampedLambdaRanges.push_back(
      std::make_shared<TFloat64Range>(
          paramStoreFFT->Get<TFloat64Range>("lambdarange")));
  inputCtx2.RebinInputs();
  BOOST_CHECK(inputCtx2.m_logGridStep == 0.0001);
  BOOST_CHECK(inputCtx2.GetSpectrum(1)->GetSpectralAxis().IsLogSampled() ==
              true);
}
BOOST_AUTO_TEST_CASE(OrthogonalizeTemplates_test) {
  CInputContext inputCtx(paramStoreOrtho);
  setInputData(inputCtx);

  inputCtx.m_categories = paramStoreNoFFT->GetList<std::string>("objects");
  inputCtx.m_lambdaRanges.push_back(std::make_shared<TFloat64Range>(
      paramStoreNoFFT->Get<TFloat64Range>("lambdarange")));
  BOOST_CHECK(catalog->GetTemplateCount("galaxy", 1, 0) == 0);

  //
  Context.LoadParameterStore(jsonStringOneSpc + jsonString + jsonStringOrtho);
  Context.setLineCatalog("galaxy", "LineModelSolve", lineCatalog);
  inputCtx.OrthogonalizeTemplates();
  BOOST_CHECK(catalog->GetTemplateCount("galaxy", 1, 0) == 1);
}

BOOST_AUTO_TEST_SUITE_END()
