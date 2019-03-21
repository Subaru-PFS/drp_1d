#include <RedshiftLibrary/linemodel/element.h>
#include <RedshiftLibrary/linemodel/multiline.h>
#include <RedshiftLibrary/ray/ray.h>

#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(test_element)

BOOST_AUTO_TEST_CASE(Instance){

  CRay ray = CRay("O2",0.1, 1, CRay::SYM, 2, 0.2, 0.3, 0.4 ,0.5 , 0.6, 0.7, "group", 0.8);
  std::vector<CRay> rs;
  rs.push_back(ray);
  std::vector<Float64> nominalAmplitudes = std::vector<Float64> ();
  nominalAmplitudes.push_back(0.8);
  std::vector<UInt32> catalogIndexes;
  catalogIndexes.push_back(1);
  catalogIndexes.push_back(0);

  CMultiLine element = CMultiLine(rs,  "fixed",  0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);
  BOOST_CHECK_CLOSE( 1.0, element.GetVelocityEmission(), 0.01 );
  BOOST_CHECK_CLOSE( 1.1, element.GetVelocityAbsorption(), 0.01 );
  BOOST_CHECK(element.GetElementTypeTag() == "CMultiLine" );
  BOOST_CHECK(element.GetSize()==2);

  BOOST_CHECK_THROW(CMultiLine(rs,  "foobar",  0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes),
		    std::runtime_error);

  element.SetVelocityEmission(2.0);
  BOOST_CHECK_CLOSE( 2.0, element.GetVelocityEmission(), 0.01 );

  element.SetVelocityAbsorption(2.1);
  BOOST_CHECK_CLOSE( 2.1, element.GetVelocityAbsorption(), 0.01 );

  element.SetAsymfitWidthCoeff(2.2);
  BOOST_CHECK_CLOSE( 2.2, element.GetAsymfitWidthCoeff(), 0.01 );

  element.SetAsymfitAlphaCoeff(2.3);
  BOOST_CHECK_CLOSE( 2.3, element.GetAsymfitAlphaCoeff(), 0.01 );

  element.SetAsymfitDelta(2.4);
  BOOST_CHECK_CLOSE( 2.4, element.GetAsymfitDelta(), 0.01 );

  element.SetSumCross(2.5);
  BOOST_CHECK_CLOSE( 2.5, element.GetSumCross(), 0.01 );

  element.SetSumGauss(2.6);
  BOOST_CHECK_CLOSE( 2.6, element.GetSumGauss(), 0.01 );

  BOOST_CHECK_CLOSE( 0., element.GetFitAmplitude(), 0.01 );
}

BOOST_AUTO_TEST_CASE(GetLineWidth){

    CRay ray = CRay("Halpha",6564.61, 2, CRay::SYM, 2,1.0, 0.5);
    std::vector<CRay> rs;
    rs.push_back(ray);
    std::vector<Float64> nominalAmplitudes = std::vector<Float64> ();
    nominalAmplitudes.push_back(0.8);
    std::vector<UInt32> catalogIndexes;
    catalogIndexes.push_back(1);
    catalogIndexes.push_back(0);

    CMultiLine elementID = CMultiLine(rs,  "instrumentdriven",  0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);
    CMultiLine elementfixed = CMultiLine(rs,  "fixed",  0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);
    CMultiLine elementcombined = CMultiLine(rs,  "combined",  0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);
    CMultiLine elementVD = CMultiLine(rs,  "velocitydriven",  0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);
    CMultiLine elementNip = CMultiLine(rs,  "nispsim2016",  0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);

    BOOST_CHECK_CLOSE( 3346.06, elementID.GetLineWidth(10000., 1., true, CRay::NONE), 0.001);
    BOOST_CHECK_CLOSE( 3346.06, elementID.GetLineWidth(10000., 1., false, CRay::NONE),0.001);
    BOOST_CHECK_CLOSE( 1.2,  elementfixed.GetLineWidth(10000., 1., true, CRay::NONE),0.001);
    BOOST_CHECK_CLOSE( 1.2, elementfixed.GetLineWidth(10000., 1., false, CRay::NONE),0.001);
    BOOST_CHECK_CLOSE( 3346.06, elementcombined.GetLineWidth(10000., 1., true, CRay::NONE),0.001);
    BOOST_CHECK_CLOSE( 3346.06, elementcombined.GetLineWidth(10000., 1., false, CRay::NONE),0.001);
    BOOST_CHECK_CLOSE( 0.0333333, elementVD.GetLineWidth(10000., 1., true, CRay::NONE), 0.001);
    BOOST_CHECK_CLOSE( 0.0366667, elementVD.GetLineWidth(10000., 1., false, CRay::NONE), 0.001);
    BOOST_CHECK_CLOSE( 6.61532, elementNip.GetLineWidth(10000., 1., true, CRay::NONE), 0.001);
    BOOST_CHECK_CLOSE( 6.61534, elementNip.GetLineWidth(10000., 1., false, CRay::NONE), 0.001);
    BOOST_CHECK_CLOSE( 600., elementNip.GetLineWidth(10000., 1., false, CRay::EXTINCT), 0.001);

}

BOOST_AUTO_TEST_CASE(GetLineProfile){
  CRay ray = CRay("Halpha",6564.61, 2, CRay::SYM, 2,1.0, 0.5);
  std::vector<CRay> rs;
  rs.push_back(ray);
  std::vector<Float64> nominalAmplitudes = std::vector<Float64> ();
  nominalAmplitudes.push_back(0.8);
  std::vector<UInt32> catalogIndexes;
  catalogIndexes.push_back(1);
  catalogIndexes.push_back(0);
  CMultiLine element = CMultiLine(rs,  "combined",  0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);

  BOOST_CHECK_CLOSE(0.237755, element.GetLineProfile(CRay::SYM,6564.61, 6568., 2. ), 0.001);
  BOOST_CHECK_CLOSE(0.944159, element.GetLineProfile(CRay::SYMXL,6564.61, 6568., 2. ), 0.001);
  BOOST_CHECK_CLOSE(0.2581961, element.GetLineProfile(CRay::LOR,6564.61, 6568., 2. ), 0.001);
  BOOST_CHECK_CLOSE(0.373055, element.GetLineProfile(CRay::ASYM,6564.61, 6565., 2. ), 0.001);
  BOOST_CHECK_CLOSE(0.0628983, element.GetLineProfile(CRay::ASYM2,6564.61, 6568., 2. ), 0.001);

  // Theses value are for asym profile with mean shift
  BOOST_CHECK_CLOSE(0.2893982, element.GetLineProfile(CRay::ASYMFIT,6564.61, 6568., 2. ), 0.001);
  BOOST_CHECK_CLOSE(0.2893982, element.GetLineProfile(CRay::ASYMFIXED,6564.61, 6568., 2. ), 0.001);

  // Theses tests failed since 56e4b5fdc567f65370075474d74f6a79ed3093c4
  // c.f. asym profile with mean shift
  //BOOST_CHECK_CLOSE(0.0628983, element.GetLineProfile(CRay::ASYMFIT,6564.61, 6568., 2. ), 0.001);
  //BOOST_CHECK_CLOSE(0.0628983, element.GetLineProfile(CRay::ASYMFIXED,6564.61, 6568., 2. ), 0.001);

}

BOOST_AUTO_TEST_CASE(GetLineProfileDerivSigma){
  CRay ray = CRay("Halpha",6564.61, 2, CRay::SYM, 2,1.0, 0.5);
  std::vector<CRay> rs;
  rs.push_back(ray);
  std::vector<Float64> nominalAmplitudes = std::vector<Float64> ();
  nominalAmplitudes.push_back(0.8);
  std::vector<UInt32> catalogIndexes;
  catalogIndexes.push_back(1);
  catalogIndexes.push_back(0);
  CMultiLine element = CMultiLine(rs,  "velocitydriven",  0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);

  BOOST_CHECK_CLOSE(0.34153872866337925, element.GetLineProfileDerivSigma(CRay::SYM, 6564.61, 6568., 2. ), 0.001);
  BOOST_CHECK_CLOSE(0.010850371757731672, element.GetLineProfileDerivSigma(CRay::SYMXL, 6564.61, 6568., 2. ), 0.001);
  BOOST_CHECK_CLOSE(0.24081246138668605, element.GetLineProfileDerivSigma(CRay::ASYM, 6564.61, 6565., 2. ), 0.001);
  BOOST_CHECK_CLOSE(0.067426590554372501, element.GetLineProfileDerivSigma(CRay::ASYM2, 6564.61, 6568., 2. ), 0.001);
  BOOST_CHECK_CLOSE(0.067426590554372501, element.GetLineProfileDerivSigma(CRay::ASYMFIT, 6564.61, 6568., 2. ), 0.001);
  BOOST_CHECK_CLOSE(0.067426590554372501, element.GetLineProfileDerivSigma(CRay::ASYMFIXED, 6564.61, 6568., 2. ), 0.001);
}


BOOST_AUTO_TEST_CASE(GetNSigmaSupport){
  CRay ray = CRay("Halpha",6564.61, 2, CRay::SYM, 2,1.0, 0.5);
  std::vector<CRay> rs;
  rs.push_back(ray);
  std::vector<Float64> nominalAmplitudes = std::vector<Float64> ();
  nominalAmplitudes.push_back(0.8);
  std::vector<UInt32> catalogIndexes;
  catalogIndexes.push_back(1);
  catalogIndexes.push_back(0);
  CMultiLine element = CMultiLine(rs,  "nispsim2016",  0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);


  BOOST_CHECK_CLOSE(8., element.GetNSigmaSupport(CRay::SYM), 0.001);
  BOOST_CHECK_CLOSE(16., element.GetNSigmaSupport(CRay::LOR), 0.001);
  BOOST_CHECK_CLOSE(8., element.GetNSigmaSupport(CRay::ASYM), 0.001);
  BOOST_CHECK_CLOSE(16., element.GetNSigmaSupport(CRay::ASYM2), 0.001);
  BOOST_CHECK_CLOSE(40., element.GetNSigmaSupport(CRay::SYMXL), 0.001);
  BOOST_CHECK_CLOSE(40., element.GetNSigmaSupport(CRay::ASYMFIT), 0.001);
  BOOST_CHECK_CLOSE(40., element.GetNSigmaSupport(CRay::ASYMFIXED), 0.001);

}



BOOST_AUTO_TEST_SUITE_END()
