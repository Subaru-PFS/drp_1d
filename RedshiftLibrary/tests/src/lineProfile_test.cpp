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
#include "RedshiftLibrary/line/lineprofileASYM.h"
#include "RedshiftLibrary/line/lineprofileASYMFIT.h"
#include "RedshiftLibrary/line/lineprofileLOR.h"
#include "RedshiftLibrary/line/lineprofileSYM.h"
#include "tests/src/tool/inputContextLight.h"

#include <boost/test/unit_test.hpp>
#include <cmath>

using namespace NSEpic;

class fixture_LineProfileTest {
public:
  fixture_LineProfileTest() {
    igmCorrectionMeiksin->convolveByLSF(LSF, TFloat64Range(1000, 12500));
  }

  TScopeStack scopeStack;
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin =
      fixture_MeiskinCorrection().igmCorrectionMeiksin;
  std::shared_ptr<CLSF> LSF =
      fixture_LSFGaussianConstantResolution(scopeStack).LSF;
  Int32 idxCount = fixture_MeiskinCorrection().idxCount;
};

BOOST_AUTO_TEST_SUITE(lineProfile_test)

BOOST_FIXTURE_TEST_CASE(lineprofile_test, fixture_LineProfileTest) {
  // default constructor for params
  TSymIgmParams igmParams;
  igmParams = TSymIgmParams();
  BOOST_CHECK(igmParams.m_igmidx == -1 && std::isnan(igmParams.m_redshift));

  // virtual functions
  CLineProfileLOR profileLOR = CLineProfileLOR(4.0);

  TAsymParams params = profileLOR.GetAsymParams();
  BOOST_CHECK(std::isnan(params.alpha) && std::isnan(params.delta) &&
              std::isnan(params.sigma));

  Float64 delta = profileLOR.GetDelta();
  BOOST_CHECK(delta == 0);

  igmParams = profileLOR.GetSymIgmParams();
  BOOST_CHECK(igmParams.m_igmidx == -1 && std::isnan(igmParams.m_redshift));

  bool res = profileLOR.isAsym();
  BOOST_CHECK(res == false);

  res = profileLOR.isAsymFit();
  BOOST_CHECK(res == false);

  res = profileLOR.isAsymFixed();
  BOOST_CHECK(res == false);

  res = profileLOR.isSymIgm();
  BOOST_CHECK(res == false);

  BOOST_CHECK(profileLOR.isSymIgmFit() == false);
  BOOST_CHECK_NO_THROW(profileLOR.SetSymIgmFit());
  BOOST_CHECK(profileLOR.isSymIgmFit() == false);
  BOOST_CHECK_NO_THROW(profileLOR.SetSymIgmFit(false));
  BOOST_CHECK(profileLOR.isSymIgmFit() == false);
  BOOST_CHECK_NO_THROW(profileLOR.SetSymIgmFixed());
  BOOST_CHECK(profileLOR.isSymIgmFit() == false);

  profileLOR.SetAsymParams(params);
  profileLOR.SetSymIgmParams(igmParams);
  profileLOR.resetParams();

  BOOST_CHECK_THROW(profileLOR.getIGMIdxCount(), GlobalException);
}

BOOST_AUTO_TEST_CASE(lineprofileSYM_test) {
  // constructor
  CLineProfileSYM profileSYM = CLineProfileSYM(4.0);
  CLineProfileSYM profileSYM2 = CLineProfileSYM(4.0, TProfile::NONE);

  // copy
  CLineProfileSYM profileSYM5(profileSYM);
  BOOST_CHECK(profileSYM5.GetNSigmaSupport() == 4.0);

  // move
  CLineProfileSYM profileSYM6(std::move(profileSYM2));

  // clone
  CLineProfile_ptr profileSYM3{
      std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM(2.))};
  CLineProfile_ptr profileSYM4(profileSYM3->Clone());
  BOOST_CHECK(profileSYM4->GetNSigmaSupport() == 2.0);
  BOOST_CHECK(profileSYM4->GetName() == TProfile::SYM);
  BOOST_CHECK(profileSYM4 != profileSYM3);

  Float64 x0 = 10.;
  Float64 sigma = 2.;

  // GetLineProfileVal
  Float64 val = profileSYM.GetLineProfileVal(10.0, x0, sigma);
  BOOST_CHECK(val == 1.0);

  // GetLineProfileDerivX0
  val = profileSYM.GetLineProfileDerivX0(10.0, x0, sigma);
  BOOST_CHECK(val == 0.0);

  val = profileSYM.GetLineProfileDerivX0(10.0 - sigma, x0, sigma);
  Float64 val2 = profileSYM.GetLineProfileDerivX0(10.0 + sigma, x0, sigma);
  BOOST_CHECK(val == -val2);

  // GetLineProfileDerivSigma
  val = profileSYM.GetLineProfileDerivSigma(10.0, x0, sigma);
  BOOST_CHECK(val == 0.0);

  val = profileSYM.GetLineProfileDerivSigma(10.0 - sigma, x0, sigma);
  val2 = profileSYM.GetLineProfileDerivSigma(10.0 + sigma, x0, sigma);
  BOOST_CHECK(val == val2);

  // GetLineFlux
  val = profileSYM.GetLineFlux(x0, sigma, 1.);
  BOOST_CHECK(val == sigma * sqrt(2 * M_PI));

  val2 = profileSYM.GetLineFlux(x0, sigma, 2.);
  BOOST_CHECK(val2 == 2 * val);

  val = profileSYM.GetLineFlux(x0, 0., 1.);
  BOOST_CHECK(val == 0.);
}

BOOST_FIXTURE_TEST_CASE(lineprofileSYMIGM_test, fixture_LineProfileTest) {
  // constructor
  CLineProfileSYMIGM profileSYMIGM =
      CLineProfileSYMIGM(igmCorrectionMeiksin, 4.0);
  CLineProfileSYM profileSYM = CLineProfileSYM(4.0);

  // clone
  CLineProfile_ptr profileSYMIGM2{std::unique_ptr<CLineProfileSYMIGM>(
      new CLineProfileSYMIGM(igmCorrectionMeiksin, 2.))};
  CLineProfile_ptr profileSYMIGM3(profileSYMIGM2->Clone());
  BOOST_CHECK(profileSYMIGM3->GetNSigmaSupport() == 2.0);
  BOOST_CHECK(profileSYMIGM3->GetName() == TProfile::SYMIGM);
  BOOST_CHECK(profileSYMIGM3 != profileSYMIGM2);

  // isSymIgm
  BOOST_CHECK(profileSYMIGM.isSymIgm() == true);

  // isSymIgmFit
  BOOST_CHECK(profileSYMIGM.isSymIgmFit() == false);
  BOOST_CHECK_NO_THROW(profileSYMIGM.SetSymIgmFit());
  BOOST_CHECK(profileSYMIGM.isSymIgmFit() == true);
  BOOST_CHECK_NO_THROW(profileSYMIGM.SetSymIgmFit(false));
  BOOST_CHECK(profileSYMIGM.isSymIgmFit() == false);
  BOOST_CHECK_NO_THROW(profileSYMIGM.SetSymIgmFit());
  BOOST_CHECK_NO_THROW(profileSYMIGM.SetSymIgmFixed());
  BOOST_CHECK(profileSYMIGM.isSymIgmFit() == false);

  // CheckMeiksinInit
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> nullIgm;
  CLineProfileSYMIGM profileSYMIGM_badIgm = CLineProfileSYMIGM(nullIgm, 4.0);
  BOOST_CHECK_THROW(profileSYMIGM_badIgm.CheckMeiksinInit(), GlobalException);

  // igmIdx, redshift
  TSymIgmParams params = {1, 2.1};
  profileSYMIGM.SetSymIgmParams(params);
  TSymIgmParams res = profileSYMIGM.GetSymIgmParams();
  BOOST_CHECK(res.m_igmidx == 1);
  BOOST_CHECK(res.m_redshift == 2.1);

  profileSYMIGM.resetParams();
  res = profileSYMIGM.GetSymIgmParams();
  BOOST_CHECK(res.m_igmidx == 1);
  BOOST_CHECK(res.m_redshift == 2.1);

  profileSYMIGM.SetSymIgmFit();
  profileSYMIGM.resetParams();
  res = profileSYMIGM.GetSymIgmParams();
  BOOST_CHECK(res.m_igmidx == -1);
  BOOST_CHECK(std::isnan(res.m_redshift));
  profileSYMIGM.SetSymIgmFit(false);

  Float64 x0 = 3769.7; //
  Float64 sigma = 2.;
  Float64 corr = 0;

  // getIGMCorrection
  corr = profileSYMIGM.getIGMCorrection(3769);
  BOOST_CHECK(corr = 1);

  profileSYMIGM.SetSymIgmParams(params);
  corr = profileSYMIGM.getIGMCorrection(3769);
  BOOST_CHECK(corr != 1); // lambda < RESTLAMBDA_LYA
  corr = profileSYMIGM.getIGMCorrection(3770);
  BOOST_CHECK(corr == 1); // lambda > RESTLAMBDA_LYA

  // getIGMIdxCount
  Int32 idxCount = profileSYMIGM.getIGMIdxCount();
  BOOST_CHECK(idxCount == idxCount);

  // GetLineProfileVal
  Float64 val = profileSYMIGM.GetLineProfileVal(x0, x0, sigma);
  BOOST_CHECK(val == 1.0);

  val = profileSYMIGM.GetLineProfileVal(x0 + 2., x0, sigma);
  Float64 val2 = profileSYM.GetLineProfileVal(x0 + 2., x0, sigma);
  BOOST_CHECK(val == val2);

  val = profileSYMIGM.GetLineProfileVal(x0 - 2., x0, sigma);
  val2 = profileSYM.GetLineProfileVal(x0 - 2., x0, sigma);
  BOOST_CHECK(val != val2);

  // GetLineProfileDerivX0
  val = profileSYMIGM.GetLineProfileDerivX0(x0, x0, sigma);
  BOOST_CHECK(val == 0.0);

  // GetLineProfileDerivSigma
  val = profileSYMIGM.GetLineProfileDerivSigma(x0, x0, sigma);
  BOOST_CHECK(val == 0.0);

  // GetLineFlux
  val = profileSYMIGM.GetLineFlux(x0, sigma, 1);
  BOOST_CHECK_CLOSE(val, 297.44547483232265, 1e-12);

  BOOST_CHECK_THROW(profileSYMIGM.GetLineFlux(x0 - 2., 0.001, 1),
                    GlobalException);

  val = profileSYMIGM.GetLineFlux(x0 + 50., sigma, 1);
  val2 = profileSYM.GetLineFlux(x0 + 50., sigma, 1);
  BOOST_CHECK(val == val2);
}

BOOST_AUTO_TEST_CASE(lineprofileLOR_test) {
  // constructor
  CLineProfileLOR profileLOR = CLineProfileLOR(4.0);

  CLineProfile_ptr profileLOR3{
      std::unique_ptr<CLineProfileLOR>(new CLineProfileLOR(2.))};
  CLineProfile_ptr profileLOR4(profileLOR3->Clone());
  BOOST_CHECK(profileLOR4->GetNSigmaSupport() == 4.0);
  BOOST_CHECK(profileLOR4->GetName() == TProfile::LOR);
  BOOST_CHECK(profileLOR4 != profileLOR3);

  Float64 x0 = 10.;
  Float64 sigma = 2.;

  // GetLineProfileVal
  Float64 val = profileLOR.GetLineProfileVal(10.0, x0, sigma);
  BOOST_CHECK(val == 1.0);

  // GetLineProfileDerivX0
  BOOST_CHECK_THROW(profileLOR.GetLineProfileDerivX0(10.0, x0, sigma),
                    GlobalException);

  // GetLineProfileDerivSigma
  BOOST_CHECK_THROW(profileLOR.GetLineProfileDerivSigma(10.0, x0, sigma),
                    GlobalException);

  // GetLineFlux
  val = profileLOR.GetLineFlux(x0, sigma, 1.);
  BOOST_CHECK(val == sigma * M_PI);

  Float64 val2 = profileLOR.GetLineFlux(x0, sigma, 2.);
  BOOST_CHECK(val2 == 2 * val);

  val = profileLOR.GetLineFlux(x0, 0., 1.);
  BOOST_CHECK(val == 0.);
}

BOOST_AUTO_TEST_CASE(lineprofileASYM_test) {
  // constructor with sigma - alpha - delta : {1.0, 4.5, 0.} & method none
  CLineProfileASYM profileASYM = CLineProfileASYM(4.0);

  TAsymParams params = {1.2, 4.6, 0.3};
  CLineProfileASYM profileASYM2 = CLineProfileASYM(4.0, params, "mean");
  CLineProfileASYM profileASYM2B = CLineProfileASYM(4.0, params, "mode");
  // alpha = 0 & delta = 0 -> profile sym
  params = {1., 0.0, 0.0};
  CLineProfileASYM profileSYM = CLineProfileASYM(4.0, params, "mean");
  CLineProfileASYM profileSYM2 = CLineProfileASYM(4.0, params, "mode");

  params = {1.2, 4.6, 0.3};
  CLineProfileASYM profileASYM_badDelta = CLineProfileASYM(4.0, params, "none");
  params = TAsymParams();
  CLineProfileASYM profileASYM_badParam = CLineProfileASYM(4.0, params, "none");

  // clone ASYM
  // nsigma support with default val = sigma * 1.0 * 1.0
  CLineProfile_ptr profileASYM3{
      std::unique_ptr<CLineProfileASYM>(new CLineProfileASYM(2.))};
  CLineProfile_ptr profileASYM4(profileASYM3->Clone());
  BOOST_CHECK(profileASYM4->GetNSigmaSupport() == 2.);
  BOOST_CHECK(profileASYM4->GetName() == TProfile::ASYM);
  BOOST_CHECK(profileASYM4 != profileASYM3);

  // convert ASYM to ASYMFIT
  CLineProfile_ptr profileASYMFIT = profileASYM2.cloneToASYMFIT();
  BOOST_CHECK(profileASYMFIT->GetNSigmaSupport() ==
              profileASYM2.GetNSigmaSupport());
  BOOST_CHECK(profileASYMFIT->GetName() == TProfile::ASYMFIT);

  CLineProfileASYM profileASYM6 = CLineProfileASYM(2.);
  CLineProfileASYMFIT profileASYMFIT2 = CLineProfileASYMFIT(profileASYM6);
  BOOST_CHECK(profileASYM6.GetNSigmaSupport() ==
              profileASYMFIT2.GetNSigmaSupport());

  // convert ASYMFIT to ASYM
  CLineProfile_ptr profileASYM7 = profileASYMFIT2.cloneToASYM();
  BOOST_CHECK(profileASYM7->GetNSigmaSupport() ==
              profileASYMFIT2.GetNSigmaSupport());
  BOOST_CHECK(profileASYM7->GetName() == TProfile::ASYM);

  CLineProfileASYMFIT profileASYMFIT3 = CLineProfileASYMFIT(2.);
  CLineProfileASYM profileASYM5 = CLineProfileASYM(profileASYMFIT3);
  BOOST_CHECK(profileASYM5.GetNSigmaSupport() ==
              profileASYMFIT3.GetNSigmaSupport());

  // clone ASYMFIT with sigma - alpha - delta : {2.0, 2.0, 0.}
  // nsigma support with default val = sigma * 2.0 * 2.5
  CLineProfile_ptr profileASYMFIT4{
      std::unique_ptr<CLineProfileASYMFIT>(new CLineProfileASYMFIT(2.))};
  CLineProfile_ptr profileASYMFIT5(profileASYMFIT4->Clone());
  BOOST_CHECK(profileASYMFIT5->GetNSigmaSupport() == 10.);
  BOOST_CHECK(profileASYMFIT5->GetName() == TProfile::ASYMFIT);
  BOOST_CHECK(profileASYMFIT5 != profileASYMFIT4);

  // Specfic method for ASYMFIT
  params = TAsymParams();
  profileASYMFIT3.SetAsymParams(params);
  BOOST_CHECK(std::isnan(profileASYMFIT3.GetDelta()));

  params = {1.2, 4.6, 0.3};
  profileASYMFIT3.SetAsymParams(params);
  BOOST_CHECK(profileASYMFIT3.GetDelta() == 0.3);

  profileASYMFIT3.resetParams();
  params = profileASYMFIT3.GetAsymParams();
  BOOST_CHECK(params.alpha == 0);
  BOOST_CHECK(params.delta == 0);
  BOOST_CHECK(params.sigma == 2);

  // GetAsymParams
  params = profileASYM.GetAsymParams();
  BOOST_CHECK(params.alpha == ASYM_DEFAULT_PARAMS.alpha);
  BOOST_CHECK(params.delta == ASYM_DEFAULT_PARAMS.delta);
  BOOST_CHECK(params.sigma == ASYM_DEFAULT_PARAMS.sigma);

  // GetDelta
  Float64 delta = profileASYM.GetDelta();
  BOOST_CHECK(delta == ASYM_DEFAULT_PARAMS.delta);

  // GetCenteringMethod
  std::string method = profileASYM2.GetCenteringMethod();
  BOOST_CHECK(method == "mean");

  // ASYM & not ASYMFIT
  BOOST_CHECK(profileASYM.isAsym() == true);
  BOOST_CHECK(profileASYM.isAsymFixed() == true);
  BOOST_CHECK(profileASYM.isAsymFit() == false);

  // ASYMFIT & not ASYMFIXED
  BOOST_CHECK(profileASYMFIT->isAsym() == true);
  BOOST_CHECK(profileASYMFIT->isAsymFixed() == false);
  BOOST_CHECK(profileASYMFIT->isAsymFit() == true);

  // isValid
  BOOST_CHECK(profileASYM.isValid() == 1);
  BOOST_CHECK(profileASYM_badParam.isValid() == 0);

  Float64 x0 = 10.;
  Float64 sigma = 2.;
  Float64 xsurc;

  // GetXSurc
  Float64 val = profileASYM.GetXSurc(0., sigma, xsurc);
  BOOST_CHECK(val == 0);
  BOOST_CHECK(xsurc == 0);

  Float64 muz = 4.6 / std::sqrt(1. + 4.6 * 4.6) * sqrt(2. / M_PI);
  val = profileASYM2.GetXSurc(0., sigma, xsurc);
  BOOST_CHECK_CLOSE(val, muz, 1e-12);

  val = profileSYM.GetXSurc(0., sigma, xsurc);
  BOOST_CHECK(xsurc == 0.);

  val = profileSYM2.GetXSurc(0., sigma, xsurc);
  BOOST_CHECK(xsurc == 0.);

  BOOST_CHECK_THROW(profileASYM_badDelta.GetXSurc(0., sigma, xsurc),
                    GlobalException);

  // GetLineProfileVal
  val = profileASYM.GetLineProfileVal(10.0, x0, sigma);
  BOOST_CHECK(val == 1.0);

  val = profileSYM.GetLineProfileVal(10.0, x0, sigma);
  BOOST_CHECK(val == 1.0);

  val = profileSYM2.GetLineProfileVal(10.0, x0, sigma);
  BOOST_CHECK(val == 1.0);

  BOOST_CHECK_THROW(profileASYM_badParam.GetLineProfileVal(10.0, x0, sigma),
                    GlobalException);

  // GetLineProfileDerivX0
  val = profileSYM.GetLineProfileDerivX0(10.0, x0, sigma);
  BOOST_CHECK(val == 0.0);

  val = profileSYM.GetLineProfileDerivX0(10.0 - sigma, x0, sigma);
  Float64 val2 = profileSYM.GetLineProfileDerivX0(10.0 + sigma, x0, sigma);
  BOOST_CHECK(val == -val2);

  BOOST_CHECK_THROW(profileASYM_badParam.GetLineProfileDerivX0(10.0, x0, sigma),
                    GlobalException);

  // GetLineProfileDerivSigma
  val = profileSYM.GetLineProfileDerivSigma(10.0, x0, sigma);
  BOOST_CHECK(val == 0.0);

  val = profileSYM.GetLineProfileDerivSigma(10.0 - sigma, x0, sigma);
  val2 = profileSYM.GetLineProfileDerivSigma(10.0 + sigma, x0, sigma);
  BOOST_CHECK_CLOSE(val, val2, 1e-12);

  BOOST_CHECK_THROW(
      profileASYM_badParam.GetLineProfileDerivSigma(10.0, x0, sigma),
      GlobalException);

  // GetLineFlux
  val = profileASYM.GetLineFlux(x0, sigma, 1.);
  BOOST_CHECK(val == sigma * sqrt(2 * M_PI));

  val2 = profileASYM.GetLineFlux(x0, sigma, 2.);
  BOOST_CHECK(val2 == 2 * val);

  val = profileASYM.GetLineFlux(x0, 0., 1.);
  BOOST_CHECK(val == 0.);
}

BOOST_AUTO_TEST_SUITE_END()