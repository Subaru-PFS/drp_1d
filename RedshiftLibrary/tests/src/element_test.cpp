#include <RedshiftLibrary/linemodel/element.h>
#include <RedshiftLibrary/linemodel/multiline.h>
#include <RedshiftLibrary/ray/ray.h>
#include "RedshiftLibrary/ray/lineprofile.h"

#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(test_element)

BOOST_AUTO_TEST_CASE(Instance){
  std::shared_ptr<CLineProfile> profilesym{std::make_shared<CLineProfileSYM>()};
  CRay ray = CRay("O2",0.1, 1, profilesym, 2, 0.2, 0.3, 0.4 ,0.5 , 0.6, 0.7, "group", 0.8);
  std::vector<CRay> rs;
  rs.push_back(ray);
  std::vector<Float64> nominalAmplitudes = std::vector<Float64> ();
  nominalAmplitudes.push_back(0.8);
  std::vector<UInt32> catalogIndexes;
  catalogIndexes.push_back(1);
  catalogIndexes.push_back(0);

  BOOST_CHECK_THROW(CMultiLine(rs,  "foobar",  8.0, 0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes),
		    std::runtime_error);

    BOOST_CHECK_THROW(CMultiLine(rs,  "fixed",  8.0, 0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes),
		    std::runtime_error);
  /*
  CMultiLine element = CMultiLine(rs,  "fixed",  8.0, 0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);
  BOOST_CHECK_CLOSE( 1.0, element.GetVelocityEmission(), 0.01 );
  BOOST_CHECK_CLOSE( 1.1, element.GetVelocityAbsorption(), 0.01 );
  BOOST_CHECK(element.GetElementTypeTag() == "CMultiLine" );
  BOOST_CHECK(element.GetSize()==2);
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
  */
}

BOOST_AUTO_TEST_CASE(GetLineWidth){
    std::shared_ptr<CLineProfile> profilesym{std::make_shared<CLineProfileSYM>()};
    CRay ray = CRay("Halpha",6564.61, 2, profilesym, 2,1.0, 0.5);
    std::vector<CRay> rs;
    rs.push_back(ray);
    std::vector<Float64> nominalAmplitudes = std::vector<Float64> ();
    nominalAmplitudes.push_back(0.8);
    std::vector<UInt32> catalogIndexes;
    catalogIndexes.push_back(1);
    catalogIndexes.push_back(0);

    CMultiLine elementID = CMultiLine(rs,  "instrumentdriven", 8.0, 0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);
    //CMultiLine elementfixed = CMultiLine(rs,  "fixed", 8.0, 0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);
    CMultiLine elementcombined = CMultiLine(rs,  "combined", 8.0, 0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);
    CMultiLine elementVD = CMultiLine(rs,  "velocitydriven", 8.0, 0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);
    //CMultiLine elementNip = CMultiLine(rs,  "nispsim2016", 8.0, 0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);

    BOOST_CHECK_CLOSE( 3346.06, elementID.GetLineWidth(10000., 1., true), 0.001);
    BOOST_CHECK_CLOSE( 3346.06, elementID.GetLineWidth(10000., 1., false),0.001);
    /*
    BOOST_CHECK_CLOSE( 1.2,  elementfixed.GetLineWidth(10000., 1., true),0.001);
    BOOST_CHECK_CLOSE( 1.2, elementfixed.GetLineWidth(10000., 1., false),0.001);
    */
    BOOST_CHECK_CLOSE( 3346.06, elementcombined.GetLineWidth(10000., 1., true),0.001);
    BOOST_CHECK_CLOSE( 3346.06, elementcombined.GetLineWidth(10000., 1., false),0.001);
    BOOST_CHECK_CLOSE( 0.0333564, elementVD.GetLineWidth(10000., 1., true), 0.001);
    BOOST_CHECK_CLOSE( 0.0366920, elementVD.GetLineWidth(10000., 1., false), 0.001);
   /*
    BOOST_CHECK_CLOSE( 6.61532, elementNip.GetLineWidth(10000., 1., true), 0.001);
    BOOST_CHECK_CLOSE( 6.61534, elementNip.GetLineWidth(10000., 1., false), 0.001);
    */
    //BOOST_CHECK_CLOSE( 600., elementNip.GetLineWidth(10000., 1., false, CRay::EXTINCT), 0.001);

}

BOOST_AUTO_TEST_CASE(GetLineProfile){
  std::shared_ptr<CLineProfile> profilesym{std::make_shared<CLineProfileSYM>()};
  CRay ray = CRay("Halpha",6564.61, 2, profilesym, 2,1.0, 0.5);
  std::vector<CRay> rs;
  rs.push_back(ray);
  std::vector<Float64> nominalAmplitudes = std::vector<Float64> ();
  nominalAmplitudes.push_back(0.8);
  std::vector<UInt32> catalogIndexes;
  catalogIndexes.push_back(1);
  catalogIndexes.push_back(0);
  CMultiLine element = CMultiLine(rs,  "combined", 8.0, 0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);

  std::shared_ptr<CLineProfile> profilelor{std::make_shared<CLineProfileLOR>()};
  std::shared_ptr<CLineProfile> profileasym{std::make_shared<CLineProfileASYM>()};
  std::shared_ptr<CLineProfile> profileasymfit{std::make_shared<CLineProfileASYMFIT>()};
  std::shared_ptr<CLineProfile> profileasymfixed{std::make_shared<CLineProfileASYMFIXED>()};

  BOOST_CHECK_CLOSE(0.237755, profilesym->GetLineProfile(6564.61, 6568., 2. ), 0.001);
  BOOST_CHECK_CLOSE(0.2581961, profilelor->GetLineProfile(6564.61, 6568., 2. ), 0.001);
  BOOST_CHECK_CLOSE(0.373055, profileasym->GetLineProfile(6564.61, 6565., 2. ), 0.001);
  // Theses value are for asym profile with mean shift
  BOOST_CHECK_CLOSE(0.781894, profileasymfit->GetLineProfile(6564.61, 6568., 2. ), 0.001);
  BOOST_CHECK_CLOSE(0.781894, profileasymfixed->GetLineProfile(6564.61, 6568., 2. ), 0.001);

}

BOOST_AUTO_TEST_CASE(GetLineProfileDerivSigma){
  
  std::shared_ptr<CLineProfile> profilesym{std::make_shared<CLineProfileSYM>()};

  CRay ray = CRay("Halpha",6564.61, 2, profilesym, 2,1.0, 0.5);
  std::vector<CRay> rs;
  rs.push_back(ray);
  std::vector<Float64> nominalAmplitudes = std::vector<Float64> ();
  nominalAmplitudes.push_back(0.8);
  std::vector<UInt32> catalogIndexes;
  catalogIndexes.push_back(1);
  catalogIndexes.push_back(0);
  CMultiLine element = CMultiLine(rs,  "velocitydriven", 8.0, 0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);

  std::shared_ptr<CLineProfile> profileasym{std::make_shared<CLineProfileASYM>()};
  std::shared_ptr<CLineProfile> profileasymfit{std::make_shared<CLineProfileASYMFIT>()};
  std::shared_ptr<CLineProfile> profileasymfixed{std::make_shared<CLineProfileASYMFIXED>()};

  BOOST_CHECK_CLOSE(0.34153872866337925, profilesym->GetLineProfileDerivSigma(6564.61, 6568., 2. ), 0.001);
  BOOST_CHECK_CLOSE(0.24081246138668605, profileasym->GetLineProfileDerivSigma( 6564.61, 6565., 2. ), 0.001);
  BOOST_CHECK_CLOSE(0.6909365681559104, profileasymfit->GetLineProfileDerivSigma( 6564.61, 6568., 2. ), 0.001);
  BOOST_CHECK_CLOSE(0.6909365681559104, profileasymfixed->GetLineProfileDerivSigma(6564.61, 6568., 2. ), 0.001);
}


BOOST_AUTO_TEST_CASE(GetNSigmaSupport){

  std::shared_ptr<CLineProfile> profilesym{std::make_shared<CLineProfileSYM>(8.0)};
  std::shared_ptr<CLineProfile> profilelor{std::make_shared<CLineProfileLOR>(16.)};
  std::shared_ptr<CLineProfile> profileasym{std::make_shared<CLineProfileASYM>(8.)};
  std::shared_ptr<CLineProfile> profileasymfixed{std::make_shared<CLineProfileASYMFIXED>(40.)};
  std::shared_ptr<CLineProfile> profileasymfit{std::make_shared<CLineProfileASYMFIT>(40.)};

  CRay ray = CRay("Halpha",6564.61, 2, profilesym, 2,1.0, 0.5);
  std::vector<CRay> rs;
  rs.push_back(ray);
  std::vector<Float64> nominalAmplitudes = std::vector<Float64> ();
  nominalAmplitudes.push_back(0.8);
  std::vector<UInt32> catalogIndexes;
  catalogIndexes.push_back(1);
  catalogIndexes.push_back(0);
  CMultiLine element = CMultiLine(rs,  "nispsim2016", 8.0, 0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);
  
  std::cout << "hey"; 
  std::cout << profilesym->GetNSigmaSupport()<<"\n";
  std::cout << profilelor->GetNSigmaSupport()<<"\n";
  std::cout << profileasym->GetNSigmaSupport()<<"\n";
  std::cout << profileasymfit->GetNSigmaSupport()<<"\n";
  std::cout << profileasymfixed->GetNSigmaSupport()<<"\n";

  BOOST_CHECK_CLOSE(8., profilesym->GetNSigmaSupport(), 0.001);
  BOOST_CHECK_CLOSE(16., profilelor->GetNSigmaSupport(), 0.001);
  BOOST_CHECK_CLOSE(8., profileasym->GetNSigmaSupport(), 0.001);
  BOOST_CHECK_CLOSE(40., profileasymfit->GetNSigmaSupport(), 0.001);
  BOOST_CHECK_CLOSE(40., profileasymfixed->GetNSigmaSupport(), 0.001);

}


BOOST_AUTO_TEST_SUITE_END()
