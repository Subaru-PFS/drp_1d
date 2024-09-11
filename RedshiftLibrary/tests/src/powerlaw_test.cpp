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
    "{\"lambdaRange\" : [ 4000, 6500 ],"; // So that xc is in range
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

// BOOST_AUTO_TEST_CASE(basicfit_simple_var) {
//   // We consider z = 0 here
//   Float64 a = 1.5e-16;
//   Float64 b = -0.5;
//   nMinSamples = 10;

//   TList<Float64> lambda {4000, 4010, 4020, 4030, 4040, 4050, 4060, 4070,
//   4080, 4090, 4100, 4110, 4120, 4130, 4140, 4150, 4160, 4170, 4180, 4190,
//   4200, 4210, 4220, 4230, 4240, 4250, 4260, 4270, 4280, 4290, 4300, 4310,
//   4320, 4330, 4340, 4350, 4360, 4370, 4380, 4390, 4400, 4410, 4420, 4430,
//   4440, 4450, 4460, 4470, 4480, 4490, 4500, 4510, 4520, 4530, 4540, 4550,
//   4560, 4570, 4580, 4590, 4600, 4610, 4620, 4630, 4640, 4650, 4660, 4670,
//   4680, 4690, 4700, 4710, 4720, 4730, 4740, 4750, 4760, 4770, 4780, 4790,
//   4800, 4810, 4820, 4830, 4840, 4850, 4860, 4870, 4880, 4890, 4900, 4910,
//   4920, 4930, 4940, 4950, 4960, 4970, 4980, 4990};
//   // TList<Float64> flux
//   {2.371708245126284e-18, 2.368749156955745e-18, 2.365801117036919e-18, 2.3628640567896592e-18,
//   2.3599379082283223e-18, 2.3570226039551584e-18, 2.354118077153792e-18, 2.351224261582788e-18,
//   2.3483410915693107e-18, 2.3454685020028577e-18, 2.342606428329091e-18, 2.339754806543738e-18,
//   2.336913573186584e-18, 2.3340826653355393e-18, 2.3312620206007845e-18, 2.3284515771189983e-18,
//   2.3256512735476584e-18, 2.3228610490594156e-18, 2.3200808433365473e-18, 2.3173105965654786e-18,
//   2.3145502494313784e-18, 2.3117997431128268e-18, 2.3090590192765472e-18, 2.3063280200722125e-18,
//   2.303606688127317e-18, 2.3008949665421113e-18, 2.298192798884608e-18, 2.2955001291856497e-18,
//   2.2928169019340366e-18, 2.2901430620717222e-18, 2.2874785549890698e-18, 2.284823326520165e-18,
//   2.282177322938192e-18, 2.2795404909508683e-18, 2.2769127776959362e-18, 2.2742941307367104e-18,
//   2.2716844980576853e-18, 2.2690838280601936e-18, 2.266492069558123e-18, 2.2639091717736827e-18,
//   2.261335084333227e-18, 2.2587697572631283e-18, 2.2562131409857007e-18, 2.2536651863151784e-18,
//   2.251125844453741e-18, 2.2485950669875843e-18, 2.2460728058830496e-18, 2.243559013482789e-18,
//   2.241053642501988e-18, 2.2385566460246264e-18, 2.2360679774997898e-18, 2.2335875907380243e-18,
//   2.231115439907736e-18, 2.2286514795316342e-18, 2.226195664483218e-18, 2.2237479499833035e-18,
//   2.2213082915965963e-18, 2.2188766452283024e-18, 2.2164529671207812e-18, 2.214037213850238e-18,
//   2.211629342323457e-18, 2.2092293097745713e-18, 2.206837073761874e-18, 2.204452592164665e-18,
//   2.2020758231801378e-18, 2.199706725320299e-18, 2.1973452574089286e-18, 2.194991378578572e-18,
//   2.192645048267573e-18, 2.190306226217134e-18, 2.1879748724684182e-18, 2.18565094735968e-18,
//   2.1833344115234323e-18, 2.1810252258836445e-18, 2.1787233516529753e-18, 2.176428750330035e-18,
//   2.1741413836966814e-18, 2.171861213815347e-18, 2.1695882030263937e-18, 2.1673223139455044e-18,
//   2.1650635094610964e-18, 2.1628117527317716e-18, 2.16056700718379e-18, 2.158329236508578e-18,
//   2.1560984046602585e-18, 2.153874475853214e-18, 2.151657414559676e-18, 2.149447185507339e-18,
//   2.1472437536770057e-18, 2.1450470843002562e-18, 2.1428571428571426e-18, 2.1406738950739127e-18,
//   2.1384973069207536e-18, 2.136327344609569e-18, 2.1341639745917706e-18, 2.132007163556104e-18,
//   2.129856878426493e-18, 2.1277130863599075e-18, 2.125575754744259e-18, 2.1234448511963156e-18};
//   // TList<Float64> flux
//   {2.351192061643426e-18, 2.8328275828168832e-18, 2.5572574789614242e-18, 2.25861506281451e-18,
//   2.3837104119670393e-18, 2.3525819102169956e-18, 2.1689410543108362e-18, 2.6270016807720085e-18,
//   2.3182003042934535e-18, 2.4458571781384045e-18, 2.2193156738371962e-18, 2.3448588800909612e-18,
//   2.2472801894308287e-18, 2.5433897957513505e-18, 2.5718541131665013e-18, 2.347783811827672e-18,
//   2.1631639049369797e-18, 2.494061415523371e-18, 2.3615233619261285e-18, 2.5452563143023004e-18,
//   2.435168916463962e-18, 2.6630005515191404e-18, 2.1715876501076054e-18, 2.0261267850228935e-18,
//   2.2921274889575585e-18, 2.2016205637329056e-18, 2.250246103086722e-18, 2.2912578280548543e-18,
//   2.0232125579401678e-18, 2.4397138224374727e-18, 2.5707711641581053e-18, 2.4407962230616223e-18,
//   2.4393479246167503e-18, 2.1176623735987592e-18, 2.4135809988980303e-18, 1.935975848683607e-18,
//   2.2066361161802933e-18, 2.28574628547126e-18, 2.196744345988345e-18, 2.459841319181855e-18,
//   2.4590139144327713e-18, 2.4173671625909085e-18, 2.3811172011942147e-18, 2.4042361394102205e-18,
//   2.5396236987383866e-18, 2.1988751442526796e-18, 2.1926382673164417e-18, 2.341982714013016e-18,
//   2.2940847023554022e-18, 2.3162845901750082e-18, 1.946604592705617e-18, 2.1672282198051887e-18,
//   2.193691171328001e-18, 2.48583344118838e-18, 1.6602565161446682e-18, 2.0681134791198003e-18,
//   2.047895644705224e-18, 2.526786496718541e-18, 2.15017046395677e-18, 2.3918625615929724e-18,
//   2.1427643768481926e-18, 2.185374342493257e-18, 2.4747660136434003e-18, 2.01541097430307e-18,
//   2.2876451410224693e-18, 2.092537065353695e-18, 2.2323767410210064e-18, 2.1387190804852368e-18,
//   2.479080088940028e-18, 1.8495615016050605e-18, 2.6757941114300463e-18, 2.02081534601526e-18,
//   1.8076272246798338e-18, 2.0394733299300884e-18, 2.5491889164424137e-18, 1.7342416694597616e-18,
//   2.4945268416265332e-18, 2.2210686066571884e-18, 2.374369420727408e-18, 1.959212583560078e-18,
//   1.944737426977144e-18, 2.058857872588191e-18, 1.9887495725499225e-18, 2.376161052147815e-18,
//   2.4233897948444134e-18, 1.8967288627836677e-18, 2.0805645288393606e-18, 2.0801956242083287e-18,
//   2.0730476687161096e-18, 2.2372908056250317e-18, 2.1318774898474037e-18, 2.154702333232344e-18,
//   2.1495334527883394e-18, 2.0969084124909334e-18, 2.366396217511721e-18, 2.006283009037979e-18,
//   2.277464361842206e-18, 1.8505148521259174e-18, 2.2820181297625883e-18, 2.1228732948454444e-18};
//   TList<Float64> flux
//   {2.338548211551598e-18, 2.375291402061136e-18, 2.340032792689709e-18, 2.363798670022807e-18,
//   2.3474173844228473e-18, 2.3877991647644225e-18, 2.345087611192833e-18, 2.324290747959023e-18,
//   2.343842987310314e-18, 2.310661520610388e-18, 2.3319967127943927e-18, 2.3565398402880014e-18,
//   2.3509968896295147e-18, 2.329518607333574e-18, 2.3541335088026202e-18, 2.3475276814361116e-18,
//   2.3005152329980288e-18, 2.3147287648324446e-18, 2.335714576806453e-18, 2.326193607343992e-18,
//   2.2890612537311866e-18, 2.296551966890712e-18, 2.2953918462435495e-18, 2.299365145900833e-18,
//   2.2969121631911147e-18, 2.302571003485119e-18, 2.3054647732100008e-18, 2.2854365853708265e-18,
//   2.321117949416075e-18, 2.2999629192588634e-18, 2.2878663036007676e-18, 2.2593851096372157e-18,
//   2.2924112921269782e-18, 2.2505091125539328e-18, 2.300229264861786e-18, 2.262898344215314e-18,
//   2.272674725267458e-18, 2.2474693023834234e-18, 2.2527238303145746e-18, 2.273455005066231e-18,
//   2.2462464826804264e-18, 2.2947547471984237e-18, 2.2753393248487003e-18, 2.272493774598556e-18,
//   2.2568742423225944e-18, 2.2595579999071467e-18, 2.2287044526270677e-18, 2.2481341313291648e-18,
//   2.256318679328381e-18, 2.2750815372406043e-18, 2.2594837471493213e-18, 2.235833610314262e-18,
//   2.2579597046887296e-18, 2.241073185706498e-18, 2.219525378750373e-18, 2.186992797999218e-18,
//   2.2305560819250526e-18, 2.2238783330090177e-18, 2.163744583976687e-18, 2.2272847834375233e-18,
//   2.1967671224621006e-18, 2.172911809810659e-18, 2.2263955376047785e-18, 2.22744397677875e-18,
//   2.2201206449959538e-18, 2.23915682334975e-18, 2.1465435783346552e-18, 2.216337921300234e-18,
//   2.215277625569155e-18, 2.185779423215783e-18, 2.199637434112605e-18, 2.2147003322936443e-18,
//   2.1495537767731402e-18, 2.193390104799245e-18, 2.1769216452177974e-18, 2.187729421790446e-18,
//   2.1779454943142484e-18, 2.1554108503334496e-18, 2.169813154669971e-18, 2.1433257598476266e-18,
//   2.1626080044330038e-18, 2.166706368713327e-18, 2.139454377095029e-18, 2.1458561251660883e-18,
//   2.1585780685556736e-18, 2.1464903893535957e-18, 2.17992128529597e-18, 2.1448712623692013e-18,
//   2.1646232918076888e-18, 2.1485761436701006e-18, 2.1498029360032138e-18, 2.141885099306405e-18,
//   2.1286886259866316e-18, 2.16107921055554e-18, 2.104462009999916e-18, 2.1227167042812475e-18,
//   2.084066437162797e-18, 2.1534181831266105e-18, 2.09249467928985e-18, 2.1414592087394346e-18};
//   TList<Float64> noise(lambda.size(), 1e-40);

//   CSpectrumSpectralAxis spectralAxis(lambda);
//   CSpectrumFluxAxis fluxAxis(flux);
//   CSpectrumNoiseAxis noiseAxis(noise);
//   fluxAxis.setError(noiseAxis);

//   // Initialize power law operator
//   std::shared_ptr<CSpectrum> spc =
//       std::make_shared<CSpectrum>(spectralAxis, fluxAxis);
//   Init(jsonString1, {spc});
//   COperatorPowerLaw operatorPowerLaw;
//   operatorPowerLaw.m_nLogSamplesMin = nMinSamples;
//   TPowerLawResult result =
//       operatorPowerLaw.BasicFit(0, false, false, nullThreshold, "simple");

//   // Accepts a 1% error for calculated coefs
//   BOOST_TEST(result.coefs.first.a == a, boost::test_tools::tolerance(0.01));
//   BOOST_TEST(result.coefs.first.b == b, boost::test_tools::tolerance(0.01));
//   BOOST_CHECK_EQUAL(result.coefs.second.a, 0);
//   BOOST_CHECK_EQUAL(result.coefs.second.b, 0);
//   BOOST_TEST(result.coefs.first.sigmaa == 1.3635324582905924e-17,
//   boost::test_tools::tolerance(0.01)); BOOST_TEST(result.coefs.first.sigmab
//   == 0.010959492066037258, boost::test_tools::tolerance(0.01));
//   Context.reset();

//   Init(jsonString1, {spc});
//   COperatorPowerLaw operatorPowerLaw2;
//   operatorPowerLaw2.m_nLogSamplesMin = nMinSamples;
//   TPowerLawResult result2 = operatorPowerLaw2.BasicFit(
//       0, false, false, nullThreshold, "simpleWeighted");

//   // Accepts a 1% error for calculated coefs
//   BOOST_TEST(result2.coefs.first.a == a, boost::test_tools::tolerance(0.01));
//   BOOST_TEST(result2.coefs.first.b == b, boost::test_tools::tolerance(0.01));
//   BOOST_CHECK_EQUAL(result2.coefs.second.a, 0);
//   BOOST_CHECK_EQUAL(result2.coefs.second.b, 0);
//   BOOST_TEST(result.coefs.first.sigmaa == 1.6956935444526124e-17,
//   boost::test_tools::tolerance(0.01)); BOOST_TEST(result.coefs.first.sigmab
//   == 0.01382943476639579, boost::test_tools::tolerance(0.01));
//   Context.reset();
// }

BOOST_AUTO_TEST_CASE(basicfit_double_without_extinction) {
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
