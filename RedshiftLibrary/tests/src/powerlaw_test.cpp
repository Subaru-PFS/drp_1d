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
#include "tests/src/tool/inputContextLight.h"
#include <boost/test/unit_test.hpp>
#include <random>
using namespace NSEpic;

class PowerLaw_fixture {
public:
  COperatorPowerLaw Init(std::string jsonString,
                         std::shared_ptr<CSpectrum> spc);
};

COperatorPowerLaw PowerLaw_fixture::Init(std::string jsonString,
                                         std::shared_ptr<CSpectrum> spc) {
  std::shared_ptr<CScopeStack> scopeStack = std::make_shared<CScopeStack>();

  Context.LoadParameterStore(jsonString);
  // TODO here how to better set flux correction meiksin / calzetti : use real
  // ones or create custom
  Context.setfluxCorrectionMeiksin(
      fixture_MeiskinCorrection().igmCorrectionMeiksin);
  Context.setfluxCorrectionCalzetti(
      fixture_CalzettiCorrection().ismCorrectionCalzetti);

  std::shared_ptr<CLSF> LSF =
      fixture_LSFGaussianConstantResolution(scopeStack).LSF;
  spc->SetLSF(LSF);
  Context.addSpectrum(spc);

  // TODO see if could be removed
  std::shared_ptr<CTemplateCatalog> catalog =
      fixture_sharedTemplateCatalog().catalog;
  Context.setTemplateCatalog(catalog);

  Context.Init();
  COperatorPowerLaw operatorPowerlaw;
  return operatorPowerlaw;
}

BOOST_FIXTURE_TEST_SUITE(powerLawOperator_test, PowerLaw_fixture)

const std::string jsonString2 =
    "{\"lambdaRange\" : [ 4000, 6500 ],"
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

BOOST_AUTO_TEST_CASE(basicfit_simple_without_extinction) {
  // Build lambda
  Float64 lambdaMin = 4000;
  Float64 lambdaMax = 6500;
  std::vector<Float64> spectralValues;
  for (Float64 lambda = lambdaMin; lambda <= lambdaMax; lambda += 10) {
    spectralValues.push_back(lambda);
  }
  CSpectrumSpectralAxis spectralAxis(spectralValues);

  // Build flux
  Int32 n = spectralValues.size();
  Float64 a = 1.5e-16;
  Float64 b = -0.5;
  std::vector<Float64> fluxValues(n);
  for (Int32 i = 0; i < n; i += 1) {
    fluxValues[i] = a * pow(spectralAxis[i], b);
  }
  CSpectrumFluxAxis fluxAxis(fluxValues);

  // Build error
  // Calculate the mean of flux and divide by 10
  Float64 meanFlux =
      std::accumulate(fluxValues.begin(), fluxValues.end(), 0.0) /
      fluxValues.size();

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

  // Spectrum
  std::shared_ptr<CSpectrum> spc =
      std::make_shared<CSpectrum>(spectralAxis, fluxAxis);

  COperatorPowerLaw operatorPowerLaw = Init(jsonString2, spc);

  Float64 nullThreshold = 2;
  Int32 nMinSamples = 100;
  TPowerLawResult result = operatorPowerLaw.BasicFit(
      0, false, false, nullThreshold, nMinSamples, "simple");

  // accepts a 1% error for calculated coefs
  BOOST_TEST(result.coefs.first.a == a, boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.first.b == b, boost::test_tools::tolerance(0.01));
  BOOST_CHECK_EQUAL(result.coefs.second.a, 0);
  BOOST_CHECK_EQUAL(result.coefs.second.b, 0);
  Context.reset();

  result = operatorPowerLaw.BasicFit(0, false, false, nullThreshold,
                                     nMinSamples, "simpleWeighted");

  // accepts a 1% error for calculated coefs
  BOOST_TEST(result.coefs.first.a == a, boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.first.b == b, boost::test_tools::tolerance(0.01));
  BOOST_CHECK_EQUAL(result.coefs.second.a, 0);
  BOOST_CHECK_EQUAL(result.coefs.second.b, 0);
  Context.reset();
}

BOOST_AUTO_TEST_CASE(basicfit_double_without_extinction) {
  // TODO see what to do whith this
  Float64 xc = 5400;
  // Build lambda
  Float64 lambdaMin = 4000;
  Float64 lambdaMax = 6500;
  std::vector<Float64> spectralValues;
  for (Float64 lambda = lambdaMin; lambda <= lambdaMax; lambda += 1) {
    spectralValues.push_back(lambda);
  }
  CSpectrumSpectralAxis spectralAxis(spectralValues);

  // Build flux
  Int32 n = spectralValues.size();
  Float64 a1 = 1.5e-16;
  Float64 b1 = -0.5;
  Float64 b2 = -0.2;
  Float64 a2 = a1 * std::pow(xc, b1 - b2);
  TFloat64List fluxValues(n);
  for (Int32 i = 0; i < n; i += 1) {
    if (spectralValues[i] < xc)
      fluxValues[i] = a1 * pow(spectralAxis[i], b1);
    else
      fluxValues[i] = a2 * pow(spectralAxis[i], b2);
  }
  CSpectrumFluxAxis fluxAxis(fluxValues);

  // Build error
  // Calculate the mean of flux and divide by 10
  Float64 meanFlux =
      std::accumulate(fluxValues.begin(), fluxValues.end(), 0.0) /
      fluxValues.size();

  // Build noise from a random number generator
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(meanFlux / 5, meanFlux / 10);
  TFloat64List noiseValues(n);
  for (Int32 i = 0; i < n; ++i) {
    noiseValues[i] = dis(gen);
  }

  CSpectrumNoiseAxis noiseAxis(noiseValues);
  fluxAxis.setError(noiseAxis);

  // Spectrum
  std::shared_ptr<CSpectrum> spc =
      std::make_shared<CSpectrum>(spectralAxis, fluxAxis);

  COperatorPowerLaw operatorPowerLaw = Init(jsonString2, spc);

  Float64 nullThreshold = 2;
  Int32 nMinSamples = 100;
  TPowerLawResult result = operatorPowerLaw.BasicFit(
      0, false, false, nullThreshold, nMinSamples, "full");

  // accepts a 1% error for calculated coefs
  BOOST_TEST(result.coefs.first.a == a1, boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.first.b == b1, boost::test_tools::tolerance(0.01));
  BOOST_TEST(result.coefs.second.a == a2,
             boost::test_tools::tolerance(10.)); // a2 has big errors
  BOOST_TEST(result.coefs.second.b == b2, boost::test_tools::tolerance(0.01));
  Context.reset();
}

BOOST_AUTO_TEST_SUITE_END()
