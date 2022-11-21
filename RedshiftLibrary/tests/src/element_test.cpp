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
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/line/line.h"
#include "RedshiftLibrary/line/lineprofile.h"
#include "RedshiftLibrary/linemodel/element.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"

#include <boost/test/unit_test.hpp>

#include <iostream>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(test_element)

BOOST_AUTO_TEST_CASE(Instance) {
  CLineProfile_ptr profilesym =
      std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM());
  CLine line = CLine("O2", 0.1, 1, std::move(profilesym), 2, 0.2, 0.3, 0.4, 0.5,
                     0.6, 0.7, "group", 0.8);
  std::vector<CLine> rs;
  rs.push_back(line);
  TFloat64List nominalAmplitudes = TFloat64List();
  nominalAmplitudes.push_back(0.8);
  TInt32List catalogIndexes;
  catalogIndexes.push_back(1);
  catalogIndexes.push_back(0);

  BOOST_CHECK_THROW(CLineModelElement(rs, "foobar", 1.0, 1.1, nominalAmplitudes,
                                      1.2, catalogIndexes),
                    GlobalException);

  BOOST_CHECK_THROW(CLineModelElement(rs, "fixed", 1.0, 1.1, nominalAmplitudes,
                                      1.2, catalogIndexes),
                    GlobalException);
  /*
  CLineModelElement element = CLineModelElement(rs,  "fixed",  8.0,
  0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes); BOOST_CHECK_CLOSE( 1.0,
  element.GetVelocityEmission(), 0.01 ); BOOST_CHECK_CLOSE( 1.1,
  element.GetVelocityAbsorption(), 0.01 );
  BOOST_CHECK(element.GetElementTypeTag() == "CLineModelElement" );
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

BOOST_AUTO_TEST_CASE(GetLineWidth) {
  CLineProfile_ptr profilesym{
      std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM())};
  CLine line = CLine("Halpha", 6564.61, 2, std::move(profilesym), 2, 1.0, 0.5);
  std::vector<CLine> rs;
  rs.push_back(line);
  TFloat64List nominalAmplitudes = TFloat64List();
  nominalAmplitudes.push_back(0.8);
  TInt32List catalogIndexes;
  catalogIndexes.push_back(1);
  catalogIndexes.push_back(0);

  CLineModelElement elementID = CLineModelElement(
      rs, "instrumentdriven", 1.0, 1.1, nominalAmplitudes, 1.2, catalogIndexes);
  // CLineModelElement elementfixed = CLineModelElement(rs,  "fixed", 1.0, 1.1,
  // nominalAmplitudes, 1.2,catalogIndexes);
  CLineModelElement elementcombined = CLineModelElement(
      rs, "combined", 1.0, 1.1, nominalAmplitudes, 1.2, catalogIndexes);
  CLineModelElement elementVD = CLineModelElement(
      rs, "velocitydriven", 1.0, 1.1, nominalAmplitudes, 1.2, catalogIndexes);
  // CLineModelElement elementNip = CLineModelElement(rs,
  // "nispsim2016", 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);

  // setLSF on multiLines
  std::string lsfType = "GaussianConstantResolution"; // TBC
  TScopeStack scopeStack;
  std::shared_ptr<CParameterStore> store =
      std::make_shared<CParameterStore>(scopeStack);
  store->Set("LSF.resolution", 0.9);
  std::shared_ptr<TLSFArguments> args =
      std::make_shared<TLSFGaussianConstantResolutionArgs>(store);
  std::shared_ptr<CLSF> lsf = LSFFactory.Create(lsfType, args);

  elementID.SetLSF(lsf);
  elementcombined.SetLSF(lsf);
  elementVD.SetLSF(lsf);

  BOOST_CHECK_CLOSE(3346.06, elementID.GetLineWidth(10000., 1., true), 0.001);
  BOOST_CHECK_CLOSE(3346.06, elementID.GetLineWidth(10000., 1., false), 0.001);
  /*
  BOOST_CHECK_CLOSE( 1.2,  elementfixed.GetLineWidth(10000., 1., true),0.001);
  BOOST_CHECK_CLOSE( 1.2, elementfixed.GetLineWidth(10000., 1., false),0.001);
  */
  BOOST_CHECK_CLOSE(3346.06, elementcombined.GetLineWidth(10000., 1., true),
                    0.001);
  BOOST_CHECK_CLOSE(3346.06, elementcombined.GetLineWidth(10000., 1., false),
                    0.001);
  BOOST_CHECK_CLOSE(0.0333564, elementVD.GetLineWidth(10000., 1., true), 0.001);
  BOOST_CHECK_CLOSE(0.0366920, elementVD.GetLineWidth(10000., 1., false),
                    0.001);
  /*
   BOOST_CHECK_CLOSE( 6.61532, elementNip.GetLineWidth(10000., 1., true),
   0.001); BOOST_CHECK_CLOSE( 6.61534, elementNip.GetLineWidth(10000., 1.,
   false), 0.001);
   */
  // BOOST_CHECK_CLOSE( 600., elementNip.GetLineWidth(10000., 1., false,
  // CLine::EXTINCT), 0.001);
}

BOOST_AUTO_TEST_CASE(GetLineProfileVal) {
  CLineProfile_ptr profilesym{
      std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM())};
  CLine line = CLine("Halpha", 6564.61, 2, profilesym->Clone(), 2, 1.0, 0.5);
  std::vector<CLine> rs;
  rs.push_back(line);
  TFloat64List nominalAmplitudes = TFloat64List();
  nominalAmplitudes.push_back(0.8);
  TInt32List catalogIndexes;
  catalogIndexes.push_back(1);
  catalogIndexes.push_back(0);

  Float64 nsigmasupport = 8.;
  Float64 resolution = 0.9;
  TAsymParams _asymParams = {1., 4.5, 0.};
  TAsymParams _asymFixedParams = {2., 2., 0.};
  TAsymParams _asymFitParams = {2., 2., 0.};

  CLineModelElement element = CLineModelElement(
      rs, "combined", 1.0, 1.1, nominalAmplitudes, 1.2, catalogIndexes);

  CLineProfile_ptr profilelor{
      std::unique_ptr<CLineProfileLOR>(new CLineProfileLOR(nsigmasupport))};
  CLineProfile_ptr profileasym{std::unique_ptr<CLineProfileASYM>(
      new CLineProfileASYM(nsigmasupport, _asymParams, "none"))};
  CLineProfile_ptr profileasymfit{std::unique_ptr<CLineProfileASYMFIT>(
      new CLineProfileASYMFIT(nsigmasupport, _asymFixedParams, "mean"))};
  CLineProfile_ptr profileasymfixed{std::unique_ptr<CLineProfileASYM>(
      new CLineProfileASYM(nsigmasupport, _asymFitParams, "mean"))};

  BOOST_CHECK_CLOSE(0.237755, profilesym->GetLineProfileVal(6564.61, 6568., 2.),
                    0.001);
  BOOST_CHECK_CLOSE(0.2581961,
                    profilelor->GetLineProfileVal(6564.61, 6568., 2.), 0.001);
  BOOST_CHECK_CLOSE(0.373055,
                    profileasym->GetLineProfileVal(6564.61, 6565., 2.), 0.001);
  // Theses value are for asym profile with mean shift
  BOOST_CHECK_CLOSE(
      0.781894, profileasymfit->GetLineProfileVal(6564.61, 6568., 2.), 0.001);
  BOOST_CHECK_CLOSE(
      0.781894, profileasymfixed->GetLineProfileVal(6564.61, 6568., 2.), 0.001);
}

BOOST_AUTO_TEST_CASE(GetLineProfileDerivSigma) {

  CLineProfile_ptr profilesym{
      std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM())};

  CLine line = CLine("Halpha", 6564.61, 2, profilesym->Clone(), 2, 1.0, 0.5);
  std::vector<CLine> rs;
  rs.push_back(line);
  TFloat64List nominalAmplitudes = TFloat64List();
  nominalAmplitudes.push_back(0.8);
  TInt32List catalogIndexes;
  catalogIndexes.push_back(1);
  catalogIndexes.push_back(0);
  CLineModelElement element = CLineModelElement(
      rs, "velocitydriven", 1.0, 1.1, nominalAmplitudes, 1.2, catalogIndexes);

  Float64 nsigmasupport = 8.;
  Float64 resolution = 0.9;
  TAsymParams _asymParams = {1., 4.5, 0.};
  TAsymParams _asymFixedParams = {2., 2., 0.};
  TAsymParams _asymFitParams = {2., 2., 0.};

  CLineProfile_ptr profileasym{std::unique_ptr<CLineProfileASYM>(
      new CLineProfileASYM(nsigmasupport, _asymParams, "none"))};
  CLineProfile_ptr profileasymfit{std::unique_ptr<CLineProfileASYMFIT>(
      new CLineProfileASYMFIT(nsigmasupport, _asymFixedParams, "mean"))};
  CLineProfile_ptr profileasymfixed{std::unique_ptr<CLineProfileASYM>(
      new CLineProfileASYM(nsigmasupport, _asymFitParams, "mean"))};

  BOOST_CHECK_CLOSE(0.34153872866337925,
                    profilesym->GetLineProfileDerivSigma(6564.61, 6568., 2.),
                    0.001);
  BOOST_CHECK_CLOSE(0.24081246138668605,
                    profileasym->GetLineProfileDerivSigma(6564.61, 6565., 2.),
                    0.001);
  BOOST_CHECK_CLOSE(
      0.6909365681559104,
      profileasymfit->GetLineProfileDerivSigma(6564.61, 6568., 2.), 0.001);
  BOOST_CHECK_CLOSE(
      0.6909365681559104,
      profileasymfixed->GetLineProfileDerivSigma(6564.61, 6568., 2.), 0.001);
}

BOOST_AUTO_TEST_CASE(GetNSigmaSupport) {

  CLineProfile_ptr profilelor{
      std::unique_ptr<CLineProfileLOR>(new CLineProfileLOR(16.))};
  CLineProfile_ptr profilesym{
      std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM(8.))};
  CLineProfile_ptr profileasymfit{
      std::unique_ptr<CLineProfileASYMFIT>(new CLineProfileASYMFIT(40.))};
  CLineProfile_ptr profileasym{
      std::unique_ptr<CLineProfileASYM>(new CLineProfileASYM(8.))};

  /*
    CLine line = CLine("Halpha",6564.61, 2, profilesym, 2,1.0, 0.5);
    std::vector<CLine> rs;
    rs.push_back(line);
    TFloat64List nominalAmplitudes = TFloat64List ();
    nominalAmplitudes.push_back(0.8);
    TInt32List catalogIndexes;
    catalogIndexes.push_back(1);
    catalogIndexes.push_back(0);
    CLineModelElement element = CLineModelElement(rs,  "nispsim2016", 8.0,
    0.9, 1.0, 1.1, nominalAmplitudes, 1.2,catalogIndexes);
    */
  /*std::cout << profilesym->GetNSigmaSupport()<<"\n";
  std::cout << profilelor->GetNSigmaSupport()<<"\n";
  std::cout << profileasym->GetNSigmaSupport()<<"\n";
  std::cout << profileasymfit->GetNSigmaSupport()<<"\n";
  std::cout << profileasymfixed->GetNSigmaSupport()<<"\n";
*/
  BOOST_CHECK_CLOSE(8., profilesym->GetNSigmaSupport(), 0.001);
  BOOST_CHECK_CLOSE(32., profilelor->GetNSigmaSupport(), 0.001);
  BOOST_CHECK_CLOSE(8., profileasym->GetNSigmaSupport(), 0.001);
  BOOST_CHECK_CLOSE(200., profileasymfit->GetNSigmaSupport(), 0.001);
}

BOOST_AUTO_TEST_SUITE_END()
