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
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/spectrum/LSFConstantResolution.h"
#include "RedshiftLibrary/spectrum/LSFConstantWidth.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"
#include "RedshiftLibrary/spectrum/LSFVariableWidth.h"
#include "RedshiftLibrary/spectrum/LSF_NISPSIM_2016.h"
#include "RedshiftLibrary/spectrum/LSF_NISPVSSPSF_201707.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(LSF_test)

bool correctMessage(const GlobalException &ex) {
  BOOST_CHECK_EQUAL(ex.what(), std::string("Size do not match "));
  return true;
}

Float64 precision = 1e-12;

BOOST_AUTO_TEST_CASE(LSF_ConstantWidth) {
  Int32 n = 10;
  CSpectrumFluxAxis FluxAxis(n, 1);
  TFloat64List lbdaList(n);
  for (Int32 i = 0; i < n; i++) {
    lbdaList[i] = i * 0.4;
  }
  CSpectrumSpectralAxis SpectralAxis(lbdaList);
  std::string lsfType = "GaussianConstantWidth";
  Float64 width = 1.09;
  TScopeStack scopeStack;
  std::shared_ptr<CParameterStore> store =
      std::make_shared<CParameterStore>(scopeStack);
  store->Set("LSF.width", width);
  std::shared_ptr<TLSFArguments> args =
      std::make_shared<TLSFGaussianConstantWidthArgs>(store);
  std::shared_ptr<CLSF> LSF = LSFFactory.Create(lsfType, args);

  // Test constructor with spectralAxis, fluxAxis and LSF
  CSpectrum object_CSpectrum = CSpectrum(SpectralAxis, FluxAxis, LSF);
  Float64 lambda = 7000.;
  BOOST_CHECK(object_CSpectrum.GetLSF()->IsValid() == true);
  BOOST_CHECK(object_CSpectrum.GetLSF()->checkAvailability(lambda) == true);
  BOOST_CHECK(object_CSpectrum.GetLSF()->getSpectralRange() ==
              TFloat64Range(LSF_MIN_LAMBDA, LSF_MAX_LAMBDA));
  BOOST_CHECK_CLOSE(object_CSpectrum.GetLSF()->GetWidth(lambda), 1.09,
                    precision);
  // Test assignment copy constructor
  CSpectrum object_CSpectrum1 = object_CSpectrum;

  BOOST_CHECK(object_CSpectrum1.GetLSF()->IsValid() == true);
  BOOST_CHECK_CLOSE(object_CSpectrum1.GetLSF()->GetWidth(lambda), 1.09,
                    precision);

  // Test copy constructor
  CSpectrum object_CSpectrum1_bis(object_CSpectrum1);

  BOOST_CHECK(object_CSpectrum1_bis.GetLSF()->IsValid() == true);
  BOOST_CHECK_CLOSE(object_CSpectrum1_bis.GetLSF()->GetWidth(lambda), 1.09,
                    precision);

  // Test constructor with spectralAxis and fluxAxis
  CSpectrum object_CSpectrum2(SpectralAxis, FluxAxis);
  BOOST_CHECK(object_CSpectrum2.GetLSF() == nullptr);

  // Test default constructor
  CSpectrum object_CSpectrum3;
  BOOST_CHECK(object_CSpectrum3.GetLSF() == nullptr);
}

BOOST_AUTO_TEST_CASE(LSF_constantWidth_test) {
  TScopeStack scopeStack;
  std::shared_ptr<CParameterStore> store =
      std::make_shared<CParameterStore>(scopeStack);
  std::string lsfType = "GaussianConstantWidth";

  // TEST OK constructor 1
  Float64 width = 1.0;
  Float64 lambda = 7000.;
  store->Set("LSF.width", width);
  std::shared_ptr<TLSFArguments> args =
      std::make_shared<TLSFGaussianConstantWidthArgs>(store);
  std::shared_ptr<CLSF> LSF = LSFFactory.Create(lsfType, args);

  BOOST_CHECK(LSF->IsValid() == true);
  BOOST_CHECK(LSF->GetWidth(lambda) == width);

  Float64 x = 8000.;
  Float64 x0 = 7999.;

  Float64 res_1 = LSF->GetProfileVal(x, x0);

  Float64 sigma = LSF->GetWidth(x0);
  Float64 res_2 = LSF->GetProfileVal(x, x0, sigma);

  Float64 res_ref = exp(-0.5 * ((x - x0) / sigma) * ((x - x0) / sigma));

  BOOST_CHECK(res_1 == res_ref);
  BOOST_CHECK(res_2 == res_ref);

  Float64 const z = 2.2;
  Float64 const lambda0_rest = 7200.0;
  Float64 const lambda0_obs = lambda0_rest * (1.0 + z);
  Float64 const sigma_obs = LSF->GetWidth(lambda0_obs);
  Float64 const sigmaSupport = LSF->GetProfile().GetNSigmaSupport();
  Float64 lbdastep_rest = 1.; // value in angstrom based on calibration-igm
                              // files
  Int32 Nhalf =
      std::round(sigmaSupport * sigma_obs / (1.0 + z) / lbdastep_rest);
  Int32 len = 2 * Nhalf + 1;
  // change to observedframe
  TFloat64List lambdas_obs(len);
  for (Int32 i = 0; i < len; i++)
    lambdas_obs[i] = (lambda0_rest + (i - Nhalf) * lbdastep_rest) * (1.0 + z);

  TFloat64List res_3 =
      LSF->getNormalizedProfileVector(lambdas_obs, lambda0_obs);

  // TEST OK constructor 2
  std::shared_ptr<TLSFArguments> args_2 =
      std::make_shared<TLSFGaussianConstantWidthArgs>(width);
  std::shared_ptr<CLSF> LSF_2 = LSFFactory.Create(lsfType, args_2);

  BOOST_CHECK(LSF_2->IsValid() == true);
  BOOST_CHECK(LSF_2->GetWidth(lambda) == width);

  // TEST KO
  width = 0.0;
  store->Set("LSF.width", width);
  args = std::make_shared<TLSFGaussianConstantWidthArgs>(store);
  std::shared_ptr<CLSF> LSF_3 = LSFFactory.Create(lsfType, args);
  BOOST_CHECK(LSF_3->IsValid() == false);
  BOOST_CHECK(LSF_3->GetWidth(lambda) == 0.0);

  Float64 res_4 = LSF_3->GetProfileVal(0., 0.);
  BOOST_CHECK(std::isnan(res_4));

  BOOST_CHECK_THROW(LSF_3->getNormalizedProfileVector(lambdas_obs, lambda0_obs),
                    GlobalException);
}

BOOST_AUTO_TEST_CASE(LSF_GaussianConstantResolution_test) {
  TScopeStack scopeStack;
  std::shared_ptr<CParameterStore> store =
      std::make_shared<CParameterStore>(scopeStack);
  std::string lsfType = "GaussianConstantResolution";

  // TEST OK constructor 1
  Float64 resolution = 100.0;
  Float64 lambda = 7000.;
  const Float64 instrumentResolutionEmpiricalFactor = 230.0 / 325.0 / 2.35;
  store->Set("LSF.resolution", resolution);
  std::shared_ptr<TLSFArguments> args =
      std::make_shared<TLSFGaussianConstantResolutionArgs>(store);
  std::shared_ptr<CLSF> LSF = LSFFactory.Create(lsfType, args);

  Float64 width_ref = lambda / resolution * instrumentResolutionEmpiricalFactor;

  BOOST_CHECK(LSF->IsValid() == true);
  BOOST_CHECK(LSF->GetWidth(lambda) == width_ref);
  BOOST_CHECK(LSF->GetWidth(lambda, true) == width_ref);

  Float64 res_ref = 7200. / 22. * instrumentResolutionEmpiricalFactor;
  Float64 res_1 = CLSFGaussianConstantResolution::computeResolution(7200, 22.);
  BOOST_CHECK(res_1 == res_ref);

  // TEST OK constructor 2
  std::shared_ptr<TLSFArguments> args_2 =
      std::make_shared<TLSFGaussianConstantResolutionArgs>(resolution);
  std::shared_ptr<CLSF> LSF_2 = LSFFactory.Create(lsfType, args_2);

  BOOST_CHECK(LSF_2->IsValid() == true);
  BOOST_CHECK(LSF_2->GetWidth(lambda) == width_ref);

  // TEST KO
  resolution = 0.5;
  store->Set("LSF.resolution", resolution);
  args = std::make_shared<TLSFGaussianConstantResolutionArgs>(store);
  std::shared_ptr<CLSF> LSF_3 = LSFFactory.Create(lsfType, args);
  BOOST_CHECK(LSF_3->IsValid() == false);
  BOOST_CHECK(LSF_3->GetWidth(0.0) == 0.0);
}

BOOST_AUTO_TEST_CASE(LSF_GaussianNISPVSSPSF201707_test) {
  TScopeStack scopeStack;
  std::shared_ptr<CParameterStore> store =
      std::make_shared<CParameterStore>(scopeStack);
  std::string lsfType = "GaussianNISPVSSPSF201707";

  // TEST OK
  Float64 sourcesize = 0.3;
  Float64 lambda = 7000.;
  store->Set("LSF.sourcesize", sourcesize);
  std::shared_ptr<TLSFArguments> args =
      std::make_shared<TLSFGaussianNISPVSSPSF201707Args>(store);
  std::shared_ptr<CLSF> LSF = LSFFactory.Create(lsfType, args);

  const Float64 arcsecPix = 0.355;
  const Float64 angstromPix = 13.4;
  Float64 instrumentSigma = (lambda * 3.939e-4 + 2.191);
  Float64 sourcesizeSigma = sourcesize * angstromPix / arcsecPix;
  Float64 width_ref = sqrt(instrumentSigma * instrumentSigma +
                           sourcesizeSigma * sourcesizeSigma);

  BOOST_CHECK(LSF->IsValid() == true);
  BOOST_CHECK(LSF->GetWidth(lambda) == width_ref);
  BOOST_CHECK(LSF->GetWidth(lambda, true) == width_ref);
}

BOOST_AUTO_TEST_CASE(LSF_GaussianNISPSIM2016_test) {
  std::string lsfType = "GaussianNISPSIM2016";

  // TEST OK
  Float64 lambda = 7000.;
  std::shared_ptr<TLSFArguments> args = std::make_shared<TLSFArguments>();
  std::shared_ptr<CLSF> LSF = LSFFactory.Create(lsfType, args);

  Float64 width_ref = (lambda * 8.121e-4 + 7.4248) / 2.35;

  BOOST_CHECK(LSF->IsValid() == true);
  BOOST_CHECK(LSF->GetWidth(lambda) == width_ref);
  BOOST_CHECK(LSF->GetWidth(lambda, true) == width_ref);
}

BOOST_AUTO_TEST_CASE(LSF_GaussianVariableWidth_test) {
  std::string lsfType = "GaussianVariableWidth";

  // TEST OK
  Int32 n = 10;
  TFloat64List spcAxis(n);
  TFloat64List widthList(n);
  for (Int32 i = 0; i < n; i++) {
    spcAxis[i] = 12000 + 1000 * i;
    widthList[i] = spcAxis[i] * 3.939e-4 + 2.191;
  }

  Float64 lambda = 14500.;
  std::shared_ptr<TLSFArguments> args =
      std::make_shared<TLSFGaussianVarWidthArgs>(spcAxis, widthList);
  std::shared_ptr<CLSF> LSF = LSFFactory.Create(lsfType, args);

  BOOST_CHECK(LSF->IsValid() == true);
  BOOST_CHECK(LSF->checkAvailability(13000) == true);
  BOOST_CHECK(LSF->getSpectralRange() == TFloat64Range(12000.0, 21000.0));
  for (Int32 i = 0; i < spcAxis.size(); i++) {
    BOOST_CHECK(LSF->GetWidth(spcAxis[i]) == widthList[i]);
  }

  Float64 width_ref = lambda * 3.939e-4 + 2.191;
  BOOST_CHECK_CLOSE(LSF->GetWidth(lambda), width_ref, precision);
  BOOST_CHECK(LSF->GetWidth(lambda, true) == width_ref);
  width_ref = spcAxis[0] * 3.939e-4 + 2.191;
  BOOST_CHECK(LSF->GetWidth(10000.0, true) == width_ref);
  width_ref = spcAxis[n - 1] * 3.939e-4 + 2.191;
  BOOST_CHECK(LSF->GetWidth(25000.0, true) == width_ref);

  // TEST KO
  TFloat64List widthList_2(n + 1, 0.5);
  args = std::make_shared<TLSFGaussianVarWidthArgs>(spcAxis, widthList_2);

  BOOST_CHECK_THROW(LSFFactory.Create(lsfType, args), GlobalException);

  widthList_2.resize(0);
  args = std::make_shared<TLSFGaussianVarWidthArgs>(spcAxis, widthList_2);
  BOOST_CHECK_THROW(LSFFactory.Create(lsfType, args), GlobalException);

  TFloat64List widthList_3(n, 0.0);
  args = std::make_shared<TLSFGaussianVarWidthArgs>(spcAxis, widthList_3);
  std::shared_ptr<CLSF> LSF_2 = LSFFactory.Create(lsfType, args);
  BOOST_CHECK(LSF_2->IsValid() == false);

  bool available = LSF_2->checkAvailability(99.);
  BOOST_CHECK(available == false);

  BOOST_CHECK_THROW(LSF->GetWidth(99.), GlobalException);
}
}
