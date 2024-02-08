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

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/spectrum/logrebinning.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "tests/src/tool/inputContextLight.h"

using namespace NSEpic;

const std::string jsonString =
    "{\"lambdaRange\" : [4680.0, 4712.0], "
    "\"smoothWidth\" : 0.5, "
    "\"redshiftStep\" : 0.001, "
    "\"star\" : { \"redshiftRange\" : [ 2.84, 2.88 ], "
    "\"stages\": [\"redshiftSolver\"],"
    "\"redshiftSolver\": {"
    "\"method\" : \"lineModelSolve\", "
    "\"lineModelSolve\" : {\"lineModel\" : { \"firstPass\" : { "
    "\"largeGridStepRatio\" : 1 }}}}}, "
    "\"galaxy\" : { \"redshiftRange\" : [ 2.84, 2.88 ], \"redshiftStep\" : "
    "0.0001, "
    "\"stages\": [\"redshiftSolver\"],"
    "\"redshiftSolver\": {"
    "\"method\" : \"lineModelSolve\","
    "\"lineModelSolve\" : {\"lineModel\" : { \"firstPass\" : { "
    "\"largeGridStepRatio\" : 1 },"
    "\"continuumFit\" : {\"fftProcessing\" : true }}}}, "
    "\"continuumRemoval\" : { \"medianKernelWidth\" : 74.0, "
    "\"medianEvenReflection\" : false, "
    "\"method\" : \"irregularSamplingMedian\"}}}";
class fixture_logRebinningTest {
public:
  TScopeStack scopeStack;
  std::shared_ptr<CParameterStore> paramStore =
      fixture_ParamStore(jsonString, scopeStack).paramStore;
  std::shared_ptr<CInputContext> ctx_logSampled =
      fixture_InputContext(jsonString, paramStore).ctx;
  std::shared_ptr<CInputContext> ctx_notLogSampled =
      fixture_InputContext2(jsonString, paramStore).ctx;
  std::shared_ptr<CTemplate> tplStar_logSampled =
      fixture_SharedStarTemplate().tpl;
  std::shared_ptr<CTemplate> tplStar_notLogSampled =
      fixture_SharedStarNotLogTemplate().tpl;
  std::shared_ptr<CTemplateCatalog> catalog =
      fixture_sharedTemplateCatalog().catalog;

  // fill catalog for those tests
  // one template 'star' not rebinned
  // 2 templates 'galaxy' : one not rebinned on already
  // rebinned
  void fillCatalog() {
    catalog->Add(fixture_SharedStarNotLogTemplate().tpl);
    catalog->Add(fixture_SharedGalaxyTemplate().tpl);
    catalog->m_logsampling = 1;
    catalog->Add(tplStar_logSampled);
  }
};

BOOST_FIXTURE_TEST_SUITE(LogRebinning_test, fixture_logRebinningTest)

Float64 precision = 1e-12;

BOOST_AUTO_TEST_CASE(Constructor_test) {
  CSpectrumLogRebinning logRebinning(*ctx_logSampled);
  BOOST_CHECK(logRebinning.m_logGridStep == ctx_logSampled->m_logGridStep);

  CSpectrumLogRebinning logRebinningNotLog(*ctx_notLogSampled);
  BOOST_CHECK(logRebinningNotLog.m_logGridStep ==
              ctx_notLogSampled->m_logGridStep);
}

BOOST_AUTO_TEST_CASE(setupRebinning_test) {
  // use context just for create CSpectrumLogRebinning object
  CSpectrumLogRebinning logRebinning(*ctx_logSampled);

  // logSampled spectrum
  CSpectrumAxis lbdaAxis({exp(2), exp(3), exp(4)});
  CSpectrumFluxAxis fluxAxis({0, 0, 0});
  CSpectrum spc_1(lbdaAxis, fluxAxis);
  TFloat64Range lbdaRange(exp(3), exp(6));

  logRebinning.m_logGridStep = 0.1;
  BOOST_CHECK_THROW(logRebinning.setupRebinning(spc_1, lbdaRange),
                    GlobalException);

  logRebinning.m_logGridStep = 1.;
  logRebinning.setupRebinning(spc_1, lbdaRange);
  BOOST_CHECK_CLOSE(logRebinning.m_lambdaRange_ref.GetBegin(), exp(3.),
                    precision);
  BOOST_CHECK_CLOSE(logRebinning.m_lambdaRange_ref.GetEnd(), exp(6), precision);

  // no logSampled spectrum
  CSpectrumAxis lbdaAxis2({2, 4, 6});
  spc_1.SetSpectralAxis(lbdaAxis2);
  TFloat64Range lbdaRange2(exp(3), exp(15));

  logRebinning.m_logGridStep = 1.;
  logRebinning.setupRebinning(spc_1, lbdaRange2);
  BOOST_CHECK_CLOSE(logRebinning.m_lambdaRange_ref.GetBegin(), exp(3),
                    precision);
  BOOST_CHECK_CLOSE(logRebinning.m_lambdaRange_ref.GetEnd(), exp(15),
                    precision);

  logRebinning.m_logGridStep = 1.1;
  logRebinning.setupRebinning(spc_1, lbdaRange2);
  BOOST_CHECK_CLOSE(logRebinning.m_lambdaRange_ref.GetBegin(), exp(3),
                    precision);
  BOOST_CHECK_CLOSE(logRebinning.m_lambdaRange_ref.GetEnd(), exp(14),
                    precision);

  logRebinning.m_logGridStep = 1.9;
  logRebinning.setupRebinning(spc_1, lbdaRange2);
  BOOST_CHECK_CLOSE(logRebinning.m_lambdaRange_ref.GetBegin(), exp(3),
                    precision);
  BOOST_CHECK_CLOSE(logRebinning.m_lambdaRange_ref.GetEnd(), exp(14.4),
                    precision);
}

BOOST_AUTO_TEST_CASE(computeTargetLogSpectralAxis_test) {
  // use context just for create CSpectrumLogRebinning object
  CSpectrumLogRebinning logRebinning(*ctx_logSampled);

  TFloat64Range lbdaRange(1, 10000);
  TFloat64List tgtRef = lbdaRange.SpreadOverLogEpsilon(1.);
  logRebinning.m_logGridStep = 1.;
  CSpectrumSpectralAxis tgtAxis =
      logRebinning.computeTargetLogSpectralAxis(lbdaRange, 10);

  BOOST_CHECK(tgtAxis.GetSamplesCount() == 10);
  for (std::size_t i = 0; i < tgtAxis.GetSamplesCount(); i++)
    BOOST_CHECK_CLOSE(tgtAxis.GetSamplesVector()[i], tgtRef[i], 1e-6);

  BOOST_CHECK_THROW(logRebinning.computeTargetLogSpectralAxis(lbdaRange, 8),
                    GlobalException);
}

BOOST_AUTO_TEST_CASE(loglambdaRebinSpectrum_test) {
  // create logSampled spectrum
  CSpectrumLogRebinning logRebinning(*ctx_logSampled);

  BOOST_CHECK_THROW(
      logRebinning.loglambdaRebinSpectrum(*ctx_logSampled->GetSpectrum(), "no"),
      GlobalException);

  // create not logSampled spectrum
  CSpectrumLogRebinning logRebinningNotLog(*ctx_notLogSampled);
  ctx_notLogSampled->GetSpectrum()->SetName("spc_notLog");

  std::shared_ptr<CSpectrum> spcLogRebinning =
      logRebinningNotLog.loglambdaRebinSpectrum(
          *ctx_notLogSampled->GetSpectrum(), "no");
  BOOST_CHECK(spcLogRebinning->GetName() == "spc_notLog");
  BOOST_CHECK(spcLogRebinning->GetSpectralAxis().IsLogSampled() == true);
  TFloat64Range lbdaRange(4680.4680234007774, 4711.9324472744638);
  BOOST_CHECK_CLOSE(
      spcLogRebinning->GetSpectralAxis().GetSamplesVector().front(),
      lbdaRange.GetBegin(), 1e-8);
  BOOST_CHECK_CLOSE(
      spcLogRebinning->GetSpectralAxis().GetSamplesVector().back(),
      lbdaRange.GetEnd(), 1e-8);

  ctx_notLogSampled->GetSpectrum()->SetSpectralAndFluxAxes(
      CSpectrumSpectralAxis(TFloat64List{1212, 1212.4, 1213}),
      CSpectrumFluxAxis(TFloat64List{0, 0, 0}));
  BOOST_CHECK_THROW(logRebinningNotLog.loglambdaRebinSpectrum(
                        *ctx_notLogSampled->GetSpectrum(), "no"),
                    GlobalException);
}

BOOST_AUTO_TEST_CASE(isRebinningNeeded_test) {
  TFloat64Range lbdaRange(4680.282, 4712.085);
  CSpectrumLogRebinning logRebinning(*ctx_logSampled);
  logRebinning.m_logGridStep = 0.8;

  // 1st case
  bool res = logRebinning.isRebinningNeeded(tplStar_logSampled, lbdaRange);
  BOOST_CHECK(res == true);

  // change logGridStep to match with tpl step to skip 1st case
  logRebinning.m_logGridStep =
      tplStar_logSampled->GetSpectralAxis().GetlogGridStep();

  // 2nd case
  lbdaRange.SetBegin(4680.);
  res = logRebinning.isRebinningNeeded(tplStar_logSampled, lbdaRange);
  BOOST_CHECK(res == true);

  // 3rd case
  lbdaRange.Set(4680.282, 4713);
  res = logRebinning.isRebinningNeeded(tplStar_logSampled, lbdaRange);
  BOOST_CHECK(res == true);

  // last case
  lbdaRange.Set(4680.482, 4712.085);
  res = logRebinning.isRebinningNeeded(tplStar_logSampled, lbdaRange);
  BOOST_CHECK(res == true);

  // not rebin needed
  lbdaRange.Set(4680.282, 4712.085);
  res = logRebinning.isRebinningNeeded(tplStar_logSampled, lbdaRange);
  BOOST_CHECK(res == false);
}

BOOST_AUTO_TEST_CASE(checkTemplateAlignment_test) {
  // create logSampled spectrum
  CSpectrumLogRebinning logRebinning(*ctx_logSampled);

  TFloat64Range lbdaRange(4680.282, 4712.085);

  bool res = logRebinning.checkTemplateAlignment(tplStar_logSampled, lbdaRange);
  BOOST_CHECK(res == true);

  lbdaRange.SetBegin(4680.482);
  res = logRebinning.checkTemplateAlignment(tplStar_logSampled, lbdaRange);
  BOOST_CHECK(res == false);
}

BOOST_AUTO_TEST_CASE(inferTemplateRebinningSetup_test) {
  // create not logSampled spectrum
  CSpectrumLogRebinning logRebinning(*ctx_notLogSampled);

  TFloat64Range lbdaRange;
  Float64 zmin_new = 0.01;
  Float64 zmax_new = 0.03;
  Float64 log_zmin_new_p1 = log(zmin_new + 1.);
  Float64 log_zmax_new_p1 = log(zmax_new + 1.);
  Int32 nb_z = Int32(ceil((log_zmax_new_p1 - log_zmin_new_p1) /
                          ctx_notLogSampled->getLogGridStep()));
  zmax_new =
      exp(log_zmin_new_p1 + nb_z * ctx_notLogSampled->getLogGridStep()) - 1.;
  TFloat64Range zrange(zmin_new, zmax_new);

  Int32 loglambda_count_tpl =
      logRebinning.inferTemplateRebinningSetup(zrange, lbdaRange);

  Float64 loglbdamin = log(4680.0 / (1.0 + zrange.GetEnd()));
  Float64 loglbdamax = log(4712.0 / (1.0 + zrange.GetBegin()));
  BOOST_CHECK(
      loglambda_count_tpl ==
      std::round((loglbdamax - loglbdamin) / ctx_notLogSampled->m_logGridStep) +
          1);
}

BOOST_AUTO_TEST_CASE(loglambdaRebinTemplate_test) {
  // create logSampled spectrum
  CSpectrumLogRebinning logRebinning(*ctx_logSampled);

  TFloat64Range lbdaRange;
  Float64 zmin_new = 2.84;
  Float64 zmax_new = 2.88;
  Float64 log_zmin_new_p1 = log(zmin_new + 1.);
  Float64 log_zmax_new_p1 = log(zmax_new + 1.);
  Int32 nb_z = Int32(ceil((log_zmax_new_p1 - log_zmin_new_p1) /
                          ctx_logSampled->getLogGridStep()));
  zmax_new =
      exp(log_zmin_new_p1 + nb_z * ctx_logSampled->getLogGridStep()) - 1.;
  TFloat64Range zrange(zmin_new, zmax_new);

  Int32 loglambda_count_tpl =
      logRebinning.inferTemplateRebinningSetup(zrange, lbdaRange);

  std::shared_ptr<NSEpic::CTemplate> tpl = logRebinning.loglambdaRebinTemplate(
      tplStar_notLogSampled, lbdaRange, loglambda_count_tpl);
  BOOST_CHECK(tpl->GetSpectralAxis().IsLogSampled() == true);
  BOOST_CHECK_CLOSE(tpl->GetSpectralAxis().GetSamplesVector().front(),
                    lbdaRange.GetBegin(), 1e-8);
  BOOST_CHECK_CLOSE(tpl->GetSpectralAxis().GetSamplesVector().back(),
                    lbdaRange.GetEnd(), 1e-8);

  // template not overlap
  lbdaRange.Set(4679., 4712);
  BOOST_CHECK_THROW(logRebinning.loglambdaRebinTemplate(
                        tplStar_logSampled, lbdaRange, loglambda_count_tpl),
                    GlobalException);

  lbdaRange.Set(4680.282, 4713);
  BOOST_CHECK_THROW(logRebinning.loglambdaRebinTemplate(
                        tplStar_logSampled, lbdaRange, loglambda_count_tpl),
                    GlobalException);
}

BOOST_AUTO_TEST_CASE(logRebinTemplateCatalog_test) {
  // create logSampled spectrum
  fillCatalog();
  ctx_logSampled->setTemplateCatalog(catalog);
  CSpectrumLogRebinning logRebinning(*ctx_logSampled);

  BOOST_CHECK_NO_THROW(logRebinning.logRebinTemplateCatalog("star"));
  BOOST_CHECK_NO_THROW(logRebinning.logRebinTemplateCatalog("galaxy"));
}

BOOST_AUTO_TEST_SUITE_END()