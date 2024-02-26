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
#include <iostream>

#include <boost/test/unit_test.hpp>

#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/line/line.h"
#include "RedshiftLibrary/line/lineprofile.h"
#include "RedshiftLibrary/linemodel/element.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(test_element)

BOOST_AUTO_TEST_CASE(Instance) {
  CLineProfile_ptr profilesym =
      std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM());
  CLine line("O2", 0.1, CLine::EType::nType_Emission, std::move(profilesym),
             CLine::EForce::nForce_Strong, 0., false, "group", 0.8, "em", 0,
             "id0");
  CLineProfile_ptr profileAsym =
      std::unique_ptr<CLineProfileASYMFIT>(new CLineProfileASYMFIT());
  CLine lineAsym("Lya", 0.1, CLine::EType::nType_Emission,
                 std::move(profileAsym), CLine::EForce::nForce_Strong, 0.,
                 false, "group", 0.8, "em", 1, "id1");

  CLineVector rs;
  rs.push_back(std::move(line));
  rs.push_back(std::move(lineAsym));

  TLineModelElementParam_ptr fdata =
      std::make_shared<TLineModelElementParam>(rs, 1.0, 1.1);

  BOOST_CHECK_THROW(CLineModelElement(fdata, "foobar"), GlobalException);

  CLineModelElement element = CLineModelElement(fdata, "combined");
  BOOST_CHECK_CLOSE(1.0, element.getVelocityEmission(), 0.01);
  BOOST_CHECK_CLOSE(1.1, element.getVelocityAbsorption(), 0.01);
  BOOST_CHECK(element.GetElementType() == CLine::EType::nType_Emission);
  BOOST_CHECK(element.IsEmission());
  BOOST_CHECK(element.GetSize() == 2);
  element.SetVelocityEmission(2.0);
  BOOST_CHECK_CLOSE(2.0, element.getVelocityEmission(), 0.01);

  element.SetVelocityAbsorption(2.1);
  BOOST_CHECK_CLOSE(2.1, element.getVelocityAbsorption(), 0.01);

  TAsymParams const params = {2.2, 2.3, 2.4};
  element.SetAsymfitParams(params);
  TAsymParams const rparam = fdata->GetAsymfitParams(0);

  BOOST_CHECK_CLOSE(2.2, rparam.sigma, 0.01);
  BOOST_CHECK_CLOSE(2.3, rparam.alpha, 0.01);
  BOOST_CHECK_CLOSE(2.4, rparam.delta, 0.01);

  element.SetSumCross(2.5);
  BOOST_CHECK_CLOSE(2.5, element.GetSumCross(), 0.01);

  element.SetSumGauss(2.6);
  BOOST_CHECK_CLOSE(2.6, element.GetSumGauss(), 0.01);

  BOOST_CHECK(std::isnan(element.GetFittedAmplitude(0)));
}

BOOST_AUTO_TEST_CASE(GetLineWidth) {
  CLineProfile_ptr profilesym{
      std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM())};
  CLine line("Halpha", 6564.61, CLine::EType::nType_Emission,
             std::move(profilesym), CLine::EForce::nForce_Strong, 0., false,
             "group", 1.0, "velgroup", 1, "id1");

  line.setNominalAmplitude(0.8);
  CLineVector rs;
  rs.push_back(std::move(line));

  TLineModelElementParam_ptr fdata_id =
      std::make_shared<TLineModelElementParam>(rs, 1.0, 1.1);
  TLineModelElementParam_ptr fdata_c =
      std::make_shared<TLineModelElementParam>(rs, 1.0, 1.1);
  TLineModelElementParam_ptr fdata_vd =
      std::make_shared<TLineModelElementParam>(rs, 1.0, 1.1);

  CLineModelElement elementID = CLineModelElement(fdata_id, "instrumentDriven");
  CLineModelElement elementcombined = CLineModelElement(fdata_c, "combined");
  CLineModelElement elementVD = CLineModelElement(fdata_vd, "velocityDriven");

  // setLSF on multiLines
  std::string lsfType = "gaussianConstantResolution"; // TBC
  std::shared_ptr<CScopeStack> scopeStack = std::make_shared<CScopeStack>();
  std::shared_ptr<CParameterStore> store =
      std::make_shared<CParameterStore>(scopeStack);
  store->Set("lsf.resolution", 0.9);
  std::shared_ptr<TLSFArguments> args =
      std::make_shared<TLSFGaussianConstantResolutionArgs>(store);
  std::shared_ptr<CLSF> lsf = LSFFactory.Create(lsfType, args);

  elementID.SetLSF(lsf);
  elementcombined.SetLSF(lsf);
  elementVD.SetLSF(lsf);

  BOOST_CHECK_CLOSE(4718.09, elementID.GetLineWidth(10000., true), 0.001);
  BOOST_CHECK_CLOSE(4718.09, elementID.GetLineWidth(10000., false), 0.001);

  BOOST_CHECK_CLOSE(4718.09, elementcombined.GetLineWidth(10000., true), 0.001);
  BOOST_CHECK_CLOSE(4718.09, elementcombined.GetLineWidth(10000., false),
                    0.001);
  BOOST_CHECK_CLOSE(0.0333564, elementVD.GetLineWidth(10000., true), 0.001);
  BOOST_CHECK_CLOSE(0.0366920, elementVD.GetLineWidth(10000., false), 0.001);
}

BOOST_AUTO_TEST_CASE(GetLineProfileVal) {
  CLineProfile_ptr profilesym{
      std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM())};
  CLine line("Halpha", 6564.61, CLine::EType::nType_Emission,
             profilesym->Clone(), CLine::EForce::nForce_Strong, 0., false,
             "group", 1.0, "velgroup", 1, "id1");

  line.setNominalAmplitude(0.8);
  CLineVector rs;
  rs.push_back(std::move(line));

  Float64 nsigmasupport = 8.;
  Float64 resolution = 0.9;
  TAsymParams _asymParams = {1., 4.5, 0.};
  TAsymParams _asymFixedParams = {2., 2., 0.};
  TAsymParams _asymFitParams = {2., 2., 0.};
  TLineModelElementParam_ptr fdata =
      std::make_shared<TLineModelElementParam>(rs, 1.0, 1.1);

  CLineModelElement element = CLineModelElement(fdata, "combined");

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
  CLine line("Halpha", 6564.61, CLine::EType::nType_Emission,
             profilesym->Clone(), CLine::EForce::nForce_Strong, 0., false,
             "group", 1.0, "velgroup", 1, "id1");

  line.setNominalAmplitude(0.8);
  CLineVector rs;
  rs.push_back(std::move(line));
  TFloat64List nominalAmplitudes = TFloat64List();
  nominalAmplitudes.push_back(0.8);
  std::shared_ptr<TLineModelElementParam> fdata =
      std::make_shared<TLineModelElementParam>(rs, 1.0, 1.1);

  CLineModelElement element = CLineModelElement(fdata, "velocityDriven");

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

  BOOST_CHECK_CLOSE(8., profilesym->GetNSigmaSupport(), 0.001);
  BOOST_CHECK_CLOSE(32., profilelor->GetNSigmaSupport(), 0.001);
  BOOST_CHECK_CLOSE(8., profileasym->GetNSigmaSupport(), 0.001);
  BOOST_CHECK_CLOSE(200., profileasymfit->GetNSigmaSupport(), 0.001);
}

BOOST_AUTO_TEST_SUITE_END()
