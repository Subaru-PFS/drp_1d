
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
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/processflow/inputcontext.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"
#include "RedshiftLibrary/spectrum/rebin/rebin.h"
#include "RedshiftLibrary/spectrum/rebin/rebinLinear.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

static TFloat64List mySpectralList;
static TFloat64List mySpectralListNoLog;
static TFloat64List myFluxList;
static TFloat64List myNoiseList;

// Create LSF based on parameters file
class MyLSF {
public:
  static std::shared_ptr<CLSF>
  createLSF(std::shared_ptr<CParameterStore> paramStore) {
    TFloat64List wave = {1213, 1214, 1215, 1216, 1217};
    TFloat64List width = {1, 2, 3, 4, 5};
    std::string lsfType = paramStore->Get<std::string>("LSF.LSFType");
    std::shared_ptr<TLSFArguments> args;
    if (lsfType == "GaussianVariableWidth")
      args = std::make_shared<TLSFGaussianVarWidthArgs>(wave, width);
    else if (lsfType == "GaussianNISPSIM2016")
      args = std::make_shared<TLSFArguments>();
    else if (lsfType == "GaussianNISPVSSPSF201707")
      args = std::make_shared<TLSFGaussianNISPVSSPSF201707Args>(paramStore);
    else if (lsfType == "GaussianConstantResolution")
      args = std::make_shared<TLSFGaussianConstantResolutionArgs>(paramStore);
    else
      args = std::make_shared<TLSFGaussianConstantWidthArgs>(paramStore);
    return LSFFactory.Create(lsfType, args);
  }

private:
};

// create Spectral Axis
class MySpectralAxis {
public:
  static CSpectrumSpectralAxis createSpcAxis(bool isLog = true) {
    mySpectralList = {1213, 1214, 1215, 1216, 1217, 1218};
    mySpectralListNoLog = {1213, 1214, 1216, 1218, 1219, 1220};
    if (isLog)
      return CSpectrumSpectralAxis(mySpectralList);
    return CSpectrumSpectralAxis(mySpectralListNoLog);
  }
};

// create Noise Axis
class MyNoiseAxis {
public:
  static CSpectrumNoiseAxis createNoiseAxis() {
    myNoiseList = {1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4};
    return CSpectrumNoiseAxis(myNoiseList);
  }
};

// create Flux Axis
class MyFluxAxis {
public:
  static CSpectrumFluxAxis createFluxAxis() {
    myFluxList = {1e-2, 2e-2, 3e-2, 4e-2, 5e-2, 6e-2};
    CSpectrumFluxAxis fluxAxis(myFluxList);
    fluxAxis.setError(MyNoiseAxis::createNoiseAxis());
    return fluxAxis;
  }
};

// * createSpectrum : create Spectrum based on 3 previous classes
// * createSpectrum : create Spectrum based on 3 previous classes with LSF based
// on parameters file
class MySpectrum {
public:
  static CSpectrum createSpectrum(bool isLog = true) {
    return CSpectrum(MySpectralAxis::createSpcAxis(isLog),
                     MyFluxAxis::createFluxAxis());
  }

  static CSpectrum createSpectrumWithLSF(std::shared_ptr<CLSF> lsf,
                                         bool isLog = true) {
    return CSpectrum(MySpectralAxis::createSpcAxis(isLog),
                     MyFluxAxis::createFluxAxis(), lsf);
  }

  static std::shared_ptr<CSpectrum> createSharedSpectrum(bool isLog = true) {
    return std::make_shared<CSpectrum>(MySpectralAxis::createSpcAxis(isLog),
                                       MyFluxAxis::createFluxAxis());
  }
};

class MyTemplate {
public:
  static CTemplate createTemplate() {
    return CTemplate("name", "category", MySpectralAxis::createSpcAxis(),
                     MyFluxAxis::createFluxAxis());
  }
  static CTemplate createEmptyTemplate() {
    return CTemplate("name", "category");
  }
};

// create Calzetti correction
class MyCalzettiCorrection {
public:
  static std::shared_ptr<CSpectrumFluxCorrectionCalzetti> createCalzettiCorr() {
    TFloat64List lbda = {1214., 1215., 1216., 1217., 1218.};
    TFloat64List flux = {0.1, 0.2, 0.5, 0.3, 0.8};
    CalzettiCorrection calzettiCorr(lbda, flux);
    return std::make_shared<CSpectrumFluxCorrectionCalzetti>(calzettiCorr, 0.,
                                                             0.1, 10);
  }
};

// create Meiskin correction
class MyMeiskinCorrection {
public:
  static std::shared_ptr<CSpectrumFluxCorrectionMeiksin>
  createMeiskinCorr(std::shared_ptr<CLSF> lsf) {
    std::vector<MeiksinCorrection> meiskinCorr;
    TFloat64List lbda = {1214., 1215., 1216., 1217., 1218.};
    TFloat64List flux = {0.1, 0.2, 0.5, 0.3, 0.8};
    meiskinCorr.push_back(MeiksinCorrection(lbda, {flux, flux, flux}));
    meiskinCorr.push_back(MeiksinCorrection(lbda, {flux, flux, flux}));
    TFloat64List z_bins = {1, 2, 3};
    std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin =
        std::make_shared<CSpectrumFluxCorrectionMeiksin>(meiskinCorr, z_bins);
    TFloat64Range lbdaRange(1214, 1216);
    igmCorrectionMeiksin->convolveByLSF(lsf, lbdaRange);
    return igmCorrectionMeiksin;
  }
};

// Create light context with the following possibilities :
// * create LSF
// * create SpectralAxis
// * create Flux Axis
// * create Noise Axis
// * create Spectrum with or without LSF
// * create template
// * create Calzetti and Meiskin corrections
class MyInputContext : public CInputContext {
public:
  using CInputContext::CInputContext;

  std::shared_ptr<CLSF> createLSF() {
    return MyLSF::createLSF(GetParameterStore());
  }

  CSpectrumSpectralAxis createSpcAxis() {
    return MySpectralAxis::createSpcAxis();
  }

  CSpectrumNoiseAxis createNoiseAxis() {
    return MyNoiseAxis::createNoiseAxis();
  }

  CSpectrumFluxAxis createFluxAxis() { return MyFluxAxis::createFluxAxis(); }

  CSpectrum createSpectrum(bool isLog = true) {
    return MySpectrum::createSpectrum(isLog);
  }

  std::shared_ptr<CSpectrum> createSharedSpectrum(bool isLog = true) {
    return MySpectrum::createSharedSpectrum(isLog);
  }

  CTemplate createTemplate() { return MyTemplate::createTemplate(); }

  CTemplate createEmptyTemplate() { return MyTemplate::createTemplate(); }

  CSpectrum createSpectrumWithLSF(bool isLog = true) {
    return MySpectrum::createSpectrumWithLSF(
        MyLSF::createLSF(GetParameterStore()), isLog);
  }

  std::shared_ptr<CSpectrumFluxCorrectionCalzetti> createCalzettiCorr() {
    return MyCalzettiCorrection::createCalzettiCorr();
  }

  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> createMeiskinCorr() {
    return MyMeiskinCorrection::createMeiskinCorr(
        MyLSF::createLSF(GetParameterStore()));
  }
};
