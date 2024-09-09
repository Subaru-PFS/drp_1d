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
#include "RedshiftLibrary/operator/powerlaw.h"
#include "RedshiftLibrary/processflow/context.h"
#include "tests/src/powerlaw/doublevardata.h"
#include "tests/src/powerlaw/simplevardata.h"
#include "tests/src/tool/inputContextLight.h"
#include <boost/test/unit_test.hpp>

#include <random>
using namespace NSEpic;

class PowerLaw_fixture {
public:
  void Init(std::string jsonString,
            TList<std::shared_ptr<CSpectrum>> const &spc);
  void addNoiseAxis(CSpectrumFluxAxis &fluxAxis);
  CSpectrumSpectralAxis createSpectralAxis(Float64 lambdaMin, Float64 lambdaMax,
                                           Int16 step);
  CSpectrumSpectralAxis
  createSpectralAxisRest(CSpectrumSpectralAxis const &spectralAxis, Float64 z);
  CSpectrumFluxAxis createFluxAxis(CSpectrumSpectralAxis const &spectralAxis,
                                   Float64 a1, Float64 b1, Float64 b2 = NAN,
                                   Float64 xc = INFINITY);
  Float64 computea2(Float64 a1, Float64 b1, Float64 b2, Float64 xc);
  void applyIsmIgmOnSpectrum(CSpectrumSpectralAxis const &spectralAxisRest,
                             CSpectrumFluxAxis const &fluxAxis, Float64 z,
                             std::shared_ptr<CSpectrum> &spc);
};

void PowerLaw_fixture::Init(std::string jsonString,
                            TList<std::shared_ptr<CSpectrum>> const &spcList) {
  std::shared_ptr<CScopeStack> scopeStack = std::make_shared<CScopeStack>();

  Context.LoadParameterStore(jsonString);
  Context.setfluxCorrectionMeiksin(
      fixture_MeiskinCorrection().igmCorrectionMeiksin);
  Context.setfluxCorrectionCalzetti(
      fixture_CalzettiCorrection().ismCorrectionCalzetti);

  std::shared_ptr<CLSF> LSF =
      fixture_LSFGaussianConstantResolution(scopeStack).LSF;
  for (Int16 spectrumIdx = 0; spectrumIdx < spcList.size(); spectrumIdx++) {
    spcList[spectrumIdx]->SetLSF(LSF);
    Context.addSpectrum(spcList[spectrumIdx]);
  }

  // TODO see if could be removed
  std::shared_ptr<CTemplateCatalog> catalog =
      fixture_sharedTemplateCatalog().catalog;
  Context.setTemplateCatalog(catalog);

  Context.Init();
}

CSpectrumSpectralAxis PowerLaw_fixture::createSpectralAxis(Float64 lambdaMin,
                                                           Float64 lambdaMax,
                                                           Int16 step) {
  // Build lambda
  std::vector<Float64> spectralValues;
  spectralValues.reserve(std::ceil(lambdaMax - lambdaMin / step));
  for (Float64 lambda = lambdaMin; lambda <= lambdaMax; lambda += step) {
    spectralValues.push_back(lambda);
  }
  return CSpectrumSpectralAxis(spectralValues);
}

CSpectrumFluxAxis
PowerLaw_fixture::createFluxAxis(CSpectrumSpectralAxis const &spectralAxis,
                                 Float64 a1, Float64 b1, Float64 b2,
                                 Float64 xc) {
  // Build flux
  Int32 n = spectralAxis.GetSamplesVector().size();
  Float64 a2 = 0;
  if (!std::isnan(b2))
    a2 = computea2(a1, b1, b2, xc);
  TFloat64List fluxValues(n);
  for (Int32 i = 0; i < n; i += 1) {
    if (spectralAxis.GetSamplesVector()[i] < xc)
      fluxValues[i] = a1 * pow(spectralAxis[i], b1);
    else
      fluxValues[i] = a2 * pow(spectralAxis[i], b2);
  }
  return CSpectrumFluxAxis(fluxValues);
}

CSpectrumSpectralAxis PowerLaw_fixture::createSpectralAxisRest(
    CSpectrumSpectralAxis const &spectralAxis, Float64 z) {
  TList<Float64> samplesRest(spectralAxis.GetSamplesCount());
  for (Int32 sampleIdx = 0; sampleIdx < spectralAxis.GetSamplesCount();
       sampleIdx++) {
    samplesRest[sampleIdx] = spectralAxis.GetSamplesVector()[sampleIdx] /
                             (1 + z); // <-> spectralAxis/(z+1)
  }
  CSpectrumSpectralAxis spectralAxisRest = CSpectrumSpectralAxis(samplesRest);
  return spectralAxisRest;
}

void PowerLaw_fixture::applyIsmIgmOnSpectrum(
    CSpectrumSpectralAxis const &spectralAxisRest,
    CSpectrumFluxAxis const &fluxAxis, Float64 z,
    std::shared_ptr<CSpectrum> &spc) {
  CTemplate spectrumTemplate1("", "", spectralAxisRest, fluxAxis);
  spectrumTemplate1.InitIsmIgmConfig(z, Context.getFluxCorrectionCalzetti(),
                                     Context.getFluxCorrectionMeiksin());
  spectrumTemplate1.ApplyMeiksinCoeff(0);
  spectrumTemplate1.ApplyDustCoeff(0);
  spc->SetFluxAxis(spectrumTemplate1.GetFluxAxis());
}

Float64 PowerLaw_fixture::computea2(Float64 a1, Float64 b1, Float64 b2,
                                    Float64 xc) {
  return a1 * std::pow(xc, b1 - b2);
}

void PowerLaw_fixture::addNoiseAxis(CSpectrumFluxAxis &fluxAxis) {
  // Build error
  // Calculate the mean of flux and divide by 10
  Int32 n = fluxAxis.GetSamplesVector().size();
  Float64 meanFlux = std::reduce(fluxAxis.GetSamplesVector().begin(),
                                 fluxAxis.GetSamplesVector().end()) /
                     n;

  // Build noise from a random number generator
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(meanFlux / 5, meanFlux / 10);
  std::vector<Float64> noiseValues(n);
  for (Int32 i = 0; i < n; ++i) {
    noiseValues[i] = dis(gen);
  }

  CSpectrumNoiseAxis noiseAxis(noiseValues);
  fluxAxis.setError(noiseAxis);
}

BOOST_FIXTURE_TEST_SUITE(powerLawOperator_test, PowerLaw_fixture)

const std::string jsonStringEnd =
    "\"smoothWidth\" : 0.0,"
    "\"spectrumModels\" : [\"galaxy\"],"
    "\"autoCorrectInput\" : false,"
    "\"continuumRemoval\": {"
    "\"method\": \"zero\", "
    "\"medianKernelWidth\" : 40.0, "
    "\"medianEvenReflection\" : false},"
    "\"lsf\" : {\"lsfType\" : \"gaussianConstantResolution\", \"resolution\" "
    ": "
    "4300},"
    "\"templateCatalog\": {"
    "\"continuumRemoval\": {"
    "\"method\": \"zero\", "
    "\"medianKernelWidth\" : 40.0, "
    "\"medianEvenReflection\" : false"
    "}}}";

const std::string jsonStringBegin1 =
    "{\"lambdaRange\" : [ 4000, 7000 ],"; // So that xc is in range
const std::string jsonStringBegin2 =
    "{\"lambdaRange\" : [ 2700, 4200 ],"; // To be in the range of ism / igm
                                          // fixture
const std::string jsonString1 = jsonStringBegin1 + jsonStringEnd;
const std::string jsonString2 = jsonStringBegin2 + jsonStringEnd;

const Float64 nullThreshold = 2;
Int32 nMinSamples = 100;
const Float64 xc = 5400;

BOOST_AUTO_TEST_CASE(basicfit_simple_without_extinction) {
  // We consider z = 0 here
  Float64 a = 1.5e-16;
  Float64 b = -0.5;

  CSpectrumSpectralAxis spectralAxis = createSpectralAxis(
      3999, 6498,
      10); // Checks that OK is spectral axis larger than lambdaRange
  CSpectrumFluxAxis fluxAxis = createFluxAxis(spectralAxis, a, b);
  addNoiseAxis(fluxAxis);

  // Initialize power law operator
  std::shared_ptr<CSpectrum> spc =
      std::make_shared<CSpectrum>(spectralAxis, fluxAxis);
  Init(jsonString1, {spc});
  COperatorPowerLaw operatorPowerLaw;
  operatorPowerLaw.m_nLogSamplesMin = nMinSamples;
  TPowerLawResult result =
      operatorPowerLaw.BasicFit(0, false, false, nullThreshold, "simple");

  // Accepts a 1% error for calculated coefs
  BOOST_TEST(result.coefs.first.a == a, boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.first.b == b, boost::test_tools::tolerance(0.01));
  BOOST_CHECK_EQUAL(result.coefs.second.a, 0);
  BOOST_CHECK_EQUAL(result.coefs.second.b, 0);
  Context.reset();

  Init(jsonString1, {spc});
  COperatorPowerLaw operatorPowerLaw2;
  operatorPowerLaw2.m_nLogSamplesMin = nMinSamples;
  TPowerLawResult result2 = operatorPowerLaw2.BasicFit(
      0, false, false, nullThreshold, "simpleWeighted");

  // Accepts a 1% error for calculated coefs
  BOOST_TEST(result2.coefs.first.a == a, boost::test_tools::tolerance(0.01));
  BOOST_TEST(result2.coefs.first.b == b, boost::test_tools::tolerance(0.01));
  BOOST_CHECK_EQUAL(result2.coefs.second.a, 0);
  BOOST_CHECK_EQUAL(result2.coefs.second.b, 0);
  Context.reset();
}

BOOST_AUTO_TEST_CASE(basicfit_simple_var) {
  // We consider z = 0 here
  Float64 a = 1.5e-16;
  Float64 b = -0.5;
  nMinSamples = 10;

  TList<Float64> lambda = fixture_powerLawSimpleVar().lambda;
  TList<Float64> flux = fixture_powerLawSimpleVar().flux;
  TList<Float64> noise = fixture_powerLawSimpleVar().noise;

  CSpectrumSpectralAxis spectralAxis(lambda);
  CSpectrumFluxAxis fluxAxis(flux);
  CSpectrumNoiseAxis noiseAxis(noise);
  fluxAxis.setError(noiseAxis);

  // Initialize power law operator
  std::shared_ptr<CSpectrum> spc =
      std::make_shared<CSpectrum>(spectralAxis, fluxAxis);
  Init(jsonString1, {spc});
  COperatorPowerLaw operatorPowerLaw;
  operatorPowerLaw.m_nLogSamplesMin = nMinSamples;
  TPowerLawResult result =
      operatorPowerLaw.BasicFit(0, false, false, nullThreshold, "simple");

  // Accepts a 1% error compared to the results found with python notebook
  BOOST_TEST(result.coefs.first.a == 1.452346083570847e-16,
             boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.first.b == -0.4962537047556169,
             boost::test_tools::tolerance(0.01));
  BOOST_CHECK_EQUAL(result.coefs.second.a, 0);
  BOOST_CHECK_EQUAL(result.coefs.second.b, 0);
  BOOST_TEST(result.coefs.first.sigmaa == 5.634124732229189e-16,
             boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.first.sigmab == 0.4533861348305359,
             boost::test_tools::tolerance(0.01));
  Context.reset();

  Init(jsonString1, {spc});
  COperatorPowerLaw operatorPowerLaw2;
  operatorPowerLaw2.m_nLogSamplesMin = nMinSamples;
  TPowerLawResult result2 = operatorPowerLaw2.BasicFit(
      0, false, false, nullThreshold, "simpleWeighted");

  // Accepts a 1% error compared to the results found with python notebook
  BOOST_TEST(result2.coefs.first.a == 1.4963406789134072e-16,
             boost::test_tools::tolerance(0.01));
  BOOST_TEST(result2.coefs.first.b == -0.49983673257891786,
             boost::test_tools::tolerance(0.01));
  BOOST_CHECK_EQUAL(result2.coefs.second.a, 0);
  BOOST_CHECK_EQUAL(result2.coefs.second.b, 0);
  BOOST_TEST(result2.coefs.first.sigmaa == 1.0296994218492315e-17,
             boost::test_tools::tolerance(0.01));
  BOOST_TEST(result2.coefs.first.sigmab == 0.008048716415326033,
             boost::test_tools::tolerance(0.01));
  Context.reset();
}

BOOST_AUTO_TEST_CASE(basicfit_double_without_extinction) {
  nMinSamples = 10;
  Float64 a1 = 1.5e-16;
  Float64 b1 = -0.5;
  Float64 b2 = -0.2;
  Float64 a2 = computea2(a1, b1, b2, xc);
  CSpectrumSpectralAxis spectralAxis = createSpectralAxis(3999, 6498, 10);
  CSpectrumFluxAxis fluxAxis = createFluxAxis(spectralAxis, a1, b1, b2, xc);
  addNoiseAxis(fluxAxis);

  // Initialize power law operator
  std::shared_ptr<CSpectrum> spc =
      std::make_shared<CSpectrum>(spectralAxis, fluxAxis);
  Init(jsonString1, {spc});
  COperatorPowerLaw operatorPowerLaw;
  operatorPowerLaw.m_nLogSamplesMin = nMinSamples;

  TPowerLawResult result =
      operatorPowerLaw.BasicFit(0, false, false, nullThreshold, "full");

  // Accepts a 1% error for calculated coefs
  BOOST_TEST(result.coefs.first.a == a1, boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.first.b == b1, boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.second.a == a2,
             boost::test_tools::tolerance(10.)); // a2 has big errors
  BOOST_TEST(result.coefs.second.b == b2, boost::test_tools::tolerance(0.01));
  Context.reset();
}

BOOST_AUTO_TEST_CASE(basicfit_double_with_var) {
  nMinSamples = 10;
  Float64 a1 = 1.5e-16;
  Float64 b1 = -0.5;
  Float64 b2 = -0.2;
  Float64 a2 = computea2(a1, b1, b2, xc);

  TList<Float64> lambda = fixture_powerLawDoubleVar().lambda;
  TList<Float64> flux = fixture_powerLawDoubleVar().flux;
  TList<Float64> noise = fixture_powerLawDoubleVar().noise;

  CSpectrumSpectralAxis spectralAxis(lambda);
  CSpectrumFluxAxis fluxAxis(flux);
  CSpectrumNoiseAxis noiseAxis(noise);
  fluxAxis.setError(noiseAxis);

  // Initialize power law operator
  std::shared_ptr<CSpectrum> spc =
      std::make_shared<CSpectrum>(spectralAxis, fluxAxis);
  Init(jsonString1, {spc});
  COperatorPowerLaw operatorPowerLaw;
  operatorPowerLaw.m_nLogSamplesMin = nMinSamples;

  TPowerLawResult result =
      operatorPowerLaw.BasicFit(0, false, false, nullThreshold, "full");

  // Accepts a 1% error compared to the results found with python notebook
  BOOST_TEST(result.coefs.first.a == 1.3593026417419627e-16,
             boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.first.b == -0.4883729145733192,
             boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.second.a == 1.1911841337173753e-17,
             boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.second.b == -0.20508628168193455,
             boost::test_tools::tolerance(0.01));

  BOOST_TEST(result.coefs.first.sigmaa == 6.6531274177064696e-18,
             boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.first.sigmab == 0.005769929304253846,
             boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.second.sigmaa == 6.868548415411009e-19,
             boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.second.sigmab == 0.0066342065076361876,
             boost::test_tools::tolerance(0.01));
  Context.reset();
}

BOOST_AUTO_TEST_CASE(basicfit_simple_with_extinction) {
  Float64 z = 2;

  // Build lambda obs & lambda rest
  CSpectrumSpectralAxis spectralAxis = createSpectralAxis(2700, 4200, 5);
  // CSpectrumSpectralAxis spectralAxisRest =
  //     createSpectralAxis(900, 1400, 1); // <-> spectralAxis/(z+1)

  TList<Float64> samplesRest(spectralAxis.GetSamplesCount());
  for (Int32 sampleIdx = 0; sampleIdx < spectralAxis.GetSamplesCount();
       sampleIdx++) {
    samplesRest[sampleIdx] = spectralAxis.GetSamplesVector()[sampleIdx] /
                             3; // <-> spectralAxis/(z+1)
  }
  CSpectrumSpectralAxis spectralAxisRest = CSpectrumSpectralAxis(samplesRest);

  Float64 a = 1.5e-16;
  Float64 b = -0.5;
  CSpectrumFluxAxis fluxAxis = createFluxAxis(spectralAxisRest, a, b);
  addNoiseAxis(fluxAxis);

  // Initialize power law operator
  std::shared_ptr<CSpectrum> spc =
      std::make_shared<CSpectrum>(spectralAxis, fluxAxis);
  Init(jsonString2, {spc});
  COperatorPowerLaw operatorPowerLaw;
  operatorPowerLaw.m_nLogSamplesMin = nMinSamples;

  // Without extinction, with redshift
  TPowerLawResult result =
      operatorPowerLaw.BasicFit(z, false, false, nullThreshold, "simple");

  // Accepts a 1% error for calculated coefs
  BOOST_TEST(result.coefs.first.a == a, boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.first.b == b, boost::test_tools::tolerance(0.01));
  BOOST_CHECK_EQUAL(result.coefs.second.a, 0);
  BOOST_CHECK_EQUAL(result.coefs.second.b, 0);
  Context.reset();

  // With extinction

  // Initialize another power law operator
  std::shared_ptr<CSpectrum> spc2 =
      std::make_shared<CSpectrum>(spectralAxis, fluxAxis);
  Init(jsonString2, {spc2});
  COperatorPowerLaw operatorPowerLaw2;
  operatorPowerLaw2.m_nLogSamplesMin = nMinSamples;

  // Creates a template to be able to apply ism / igm
  CTemplate spectrumTemplate("", "", spectralAxisRest, fluxAxis);
  spectrumTemplate.InitIsmIgmConfig(z, Context.getFluxCorrectionCalzetti(),
                                    Context.getFluxCorrectionMeiksin());
  spectrumTemplate.ApplyMeiksinCoeff(0);
  spectrumTemplate.ApplyDustCoeff(0);

  // Adapts spectrum flux
  spc2->SetFluxAxis(spectrumTemplate.GetFluxAxis());

  result = operatorPowerLaw2.BasicFit(z, true, true, nullThreshold, "simple");

  // accepts a 1% error for calculated coefs
  BOOST_TEST(result.coefs.first.a == a, boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.first.b == b, boost::test_tools::tolerance(0.01));
  BOOST_CHECK_EQUAL(result.coefs.second.a, 0);
  BOOST_CHECK_EQUAL(result.coefs.second.b, 0);
  Context.reset();
}

BOOST_AUTO_TEST_CASE(basicfit_multiobs) {
  // Updates the xc value so it is in the ism / igm correction range
  Float64 xc2 = 1200;
  Float64 z = 2;
  Float64 a1 = 1.5e-16;
  Float64 b1 = -0.5;
  Float64 b2 = -0.2;
  Float64 a2 = computea2(a1, b1, b2, xc);

  // Build lambda obs & lambda rest
  CSpectrumSpectralAxis spectralAxis = createSpectralAxis(2700, 4200, 5);
  CSpectrumSpectralAxis spectralAxis2 = createSpectralAxis(2701, 4201, 5);

  CSpectrumSpectralAxis spectralAxisRest =
      createSpectralAxisRest(spectralAxis, z);
  CSpectrumSpectralAxis spectralAxisRest2 =
      createSpectralAxisRest(spectralAxis2, z);

  CSpectrumFluxAxis fluxAxis =
      createFluxAxis(spectralAxisRest, a1, b1, b2, xc2);
  CSpectrumFluxAxis fluxAxis2 =
      createFluxAxis(spectralAxisRest2, a1, b1, b2, xc2);

  addNoiseAxis(fluxAxis);
  addNoiseAxis(fluxAxis2);

  // Same values in both spectrum
  std::shared_ptr<CSpectrum> spc =
      std::make_shared<CSpectrum>(spectralAxis, fluxAxis);
  std::shared_ptr<CSpectrum> spc2 =
      std::make_shared<CSpectrum>(spectralAxis2, fluxAxis2);

  // Ism Igm
  applyIsmIgmOnSpectrum(spectralAxisRest, fluxAxis, z, spc);
  applyIsmIgmOnSpectrum(spectralAxisRest2, fluxAxis2, z, spc2);

  Init(jsonString2, {spc, spc2});
  COperatorPowerLaw operatorPowerLaw({}, xc2);

  operatorPowerLaw.m_nLogSamplesMin = nMinSamples;
  TPowerLawResult result =
      operatorPowerLaw.BasicFit(z, true, true, nullThreshold, "full");

  // Accepts a 1% error for calculated coefs
  BOOST_TEST(result.coefs.first.a == a1, boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.first.b == b1, boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.second.a == a2,
             boost::test_tools::tolerance(10.)); // a2 has big errors
  BOOST_TEST(result.coefs.second.b == b2, boost::test_tools::tolerance(0.01));
  Context.reset();

  // Disjoint samples
  CSpectrumSpectralAxis spectralAxis1 = createSpectralAxis(2700, 3100, 5);
  CSpectrumSpectralAxis spectralAxisRest1 =
      createSpectralAxisRest(spectralAxis1, z);
  CSpectrumFluxAxis fluxAxis1 =
      createFluxAxis(spectralAxisRest1, a1, b1, b2, xc2);
  addNoiseAxis(fluxAxis1);

  spectralAxis2 = createSpectralAxis(3150, 4200, 5);
  spectralAxisRest2 = createSpectralAxisRest(spectralAxis2, z);
  fluxAxis2 = createFluxAxis(spectralAxisRest2, a1, b1, b2, xc2);
  addNoiseAxis(fluxAxis2);

  spc = std::make_shared<CSpectrum>(spectralAxis1, fluxAxis1);
  spc2 = std::make_shared<CSpectrum>(spectralAxis2, fluxAxis2);

  Init(jsonString2, {spc, spc2});

  applyIsmIgmOnSpectrum(spectralAxisRest1, fluxAxis1, z, spc);
  applyIsmIgmOnSpectrum(spectralAxisRest2, fluxAxis2, z, spc2);

  COperatorPowerLaw operatorPowerLaw2{{}, xc2};
  operatorPowerLaw2.m_nLogSamplesMin = nMinSamples;
  result = operatorPowerLaw2.BasicFit(z, true, true, nullThreshold, "full");

  // Accepts a 1% error for calculated coefs
  BOOST_TEST(result.coefs.first.a == a1, boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.first.b == b1, boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.second.a == a2,
             boost::test_tools::tolerance(10.)); // a2 has big errors
  BOOST_TEST(result.coefs.second.b == b2, boost::test_tools::tolerance(0.01));
  Context.reset();

  // Too small samples (left)
  spectralAxis1 = createSpectralAxis(2700, 3100, 50);
  spectralAxisRest1 = createSpectralAxisRest(spectralAxis1, z);
  fluxAxis1 = createFluxAxis(spectralAxisRest1, a1, b1, b2, xc2);
  addNoiseAxis(fluxAxis1);

  spectralAxis2 = createSpectralAxis(3650, 4200, 5);
  spectralAxisRest2 = createSpectralAxisRest(spectralAxis2, z);
  fluxAxis2 = createFluxAxis(spectralAxisRest2, a1, b1, b2, xc2);
  addNoiseAxis(fluxAxis2);

  spc = std::make_shared<CSpectrum>(spectralAxis1, fluxAxis1);
  spc2 = std::make_shared<CSpectrum>(spectralAxis2, fluxAxis2);

  Init(jsonString2, {spc, spc2});

  applyIsmIgmOnSpectrum(spectralAxisRest1, fluxAxis1, z, spc);
  applyIsmIgmOnSpectrum(spectralAxisRest2, fluxAxis2, z, spc2);

  COperatorPowerLaw operatorPowerLaw3{{}, xc2};
  operatorPowerLaw3.m_nLogSamplesMin = nMinSamples;

  result = operatorPowerLaw3.BasicFit(z, true, true, nullThreshold, "full");

  // Accepts a 1% error for calculated coefs
  BOOST_TEST(result.coefs.first.a == 0);
  BOOST_TEST(result.coefs.first.b == 0);
  BOOST_TEST(result.coefs.second.a == a2,
             boost::test_tools::tolerance(10.)); // a2 has big errors
  BOOST_TEST(result.coefs.second.b == b2, boost::test_tools::tolerance(0.01));
  Context.reset();
}
BOOST_AUTO_TEST_SUITE_END()
