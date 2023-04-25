
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
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/continuum/irregularsamplingmedian.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "tests/src/tool/inputContextLight.h"
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

// build spectrum
// - spectralAxis
// - rawFluxAxis
// - LSF
// - photoData
// - to estimate :
//  - continuumFluxAxis
//  - withoutContinuumFluxAxis

void print_flux(TAxisSampleList sample) {
  BOOST_TEST_MESSAGE("=======");
  for (Int32 i = 0; i < sample.size(); i++) {
    BOOST_TEST_MESSAGE(sample[i]);
  }
  BOOST_TEST_MESSAGE("=======");
}

const std::string jsonString =
    "{\"smoothWidth\" : 0.5,"
    "\"continuumRemoval\" : { \"medianKernelWidth\" : 74.0, "
    "\"medianEvenReflection\" : false, "
    "\"method\" : \"IrregularSamplingMedian\"}}";

class fixture_SpectrumTest {
public:
  void reinitializeAxis() {
    spcAxis = fixture_SpectralAxis().spcAxis;
    fluxAxis = fixture_FluxAxis().fluxAxis;
  }

  TScopeStack scopeStack; // = fixture_ScopeStore().scopeStack;
  CSpectrum spc = fixture_Spectrum().spc;
  CSpectrum spcWithVariableWidthLSF =
      fixture_SpectrumWithLSF().spcWithVariableWidthLSF;
  CSpectrum spcLight = fixture_SpectrumLight().spcLight;
  CSpectrumFluxAxis fluxAxis = fixture_FluxAxis().fluxAxis;
  CSpectrumSpectralAxis spcAxis = fixture_SpectralAxis().spcAxis;
  CSpectrumNoiseAxis noiseAxis = fixture_NoiseAxis().noiseAxis;
  std::shared_ptr<CParameterStore> paramStore =
      fixture_ParamStore(jsonString, scopeStack).paramStore;
  std::shared_ptr<CLSF> LSF = fixture_LSFGaussianVariableWidth().LSF;

  // list
  TFloat64List fluxList = fixture_FluxAxis().fluxAxisList;
  TFloat64List spectralList = fixture_SpectralAxis().spcAxisList;
  TFloat64List noiseList = fixture_NoiseAxis().noiseList;

  // Size
  Int32 spcAxisSize = fixture_SpectralAxis().spcAxisSize;
};

BOOST_FIXTURE_TEST_SUITE(Spectrum, fixture_SpectrumTest)

BOOST_AUTO_TEST_CASE(constructor_test) {
  // create spectrum
  spc.InitSpectrumContinuum(*paramStore);
  spc.EstimateContinuum();
  TFloat64List maskList(spcAxisSize, 1);

  // constructor
  CSpectrum spc_2;
  BOOST_CHECK(spc_2.m_estimationMethod == "");
  BOOST_CHECK(spc_2.m_Name == "");
  BOOST_CHECK(spc_2.GetSampleCount() == 0);

  CSpectrum spc_3("spectra");
  BOOST_CHECK(spc_3.m_estimationMethod == "");
  BOOST_CHECK(spc_3.m_Name == "spectra");
  BOOST_CHECK(spc_3.GetSampleCount() == 0);

  maskList[0] = 0.;
  CSpectrum spc_4(spc, maskList);
  BOOST_CHECK(spc_4.GetSampleCount() == spcAxisSize - 1);
  BOOST_CHECK(spc_4.GetFluxAxis()[0] == fluxList[1]);

  // CSpectrum spc_5 = ctx.createSpectrumWithLSF();
  std::shared_ptr<const CLSF> LSF_out = spcWithVariableWidthLSF.GetLSF();
  BOOST_CHECK(spcWithVariableWidthLSF.GetSampleCount() == spcAxisSize);
  BOOST_CHECK(spcWithVariableWidthLSF.GetFluxAxis().GetSamplesVector() ==
              fluxList);
  BOOST_CHECK(LSF_out->GetWidth(1214) == 2);

  // CSpectrumSpectralAxis spectralAxis_2 = ctx.createSpcAxis();
  spcAxis.SetSize(2);
  BOOST_CHECK_THROW(CSpectrum spc_5(spcAxis, fluxAxis), GlobalException);

  // copy and copy assignement
  CSpectrum spc_6(spc);
  BOOST_CHECK(spc_6.GetSampleCount() == spcAxisSize);
  BOOST_CHECK(spc_6.GetFluxAxis().GetSamplesVector() == fluxList);

  TAxisSampleList rawFlux =
      spcWithVariableWidthLSF.GetRawFluxAxis_().GetSamplesVector();
  rawFlux.pop_back();
  spcWithVariableWidthLSF.GetRawFluxAxis_().setSamplesVector(rawFlux);
  BOOST_CHECK_THROW(CSpectrum spc_6b(spcWithVariableWidthLSF), GlobalException);

  CSpectrum spc_7;
  spc_7 = spc;
  BOOST_CHECK(spc_7.GetSampleCount() == spcAxisSize);
  BOOST_CHECK(spc_7.GetFluxAxis().GetSamplesVector() == fluxList);

  // move and move assignement
  CSpectrum spc_8(std::move(spc_6));
  BOOST_CHECK(spc_8.GetSampleCount() == spcAxisSize);
  BOOST_CHECK(spc_8.GetFluxAxis().GetSamplesVector() == fluxList);
  BOOST_CHECK(spc_6.GetSampleCount() == 0);

  BOOST_CHECK_THROW(CSpectrum spc_8b(std::move(spcWithVariableWidthLSF)),
                    GlobalException);

  CSpectrum spc_9;
  spc_9 = std::move(spc_8);
  BOOST_CHECK(spc_9.GetSampleCount() == spcAxisSize);
  BOOST_CHECK(spc_9.GetFluxAxis().GetSamplesVector() == fluxList);
  BOOST_CHECK(spc_8.GetSampleCount() == 0);
}

BOOST_AUTO_TEST_CASE(setXXX_test) {
  // add LSF
  spc.SetLSF(fixture_LSFGaussianVariableWidth().LSF);
  std::shared_ptr<const CLSF> LSF_out = spc.GetLSF();
  BOOST_CHECK(LSF_out->GetWidth(1214) == 2);

  // add photometric data
  const TStringList name = {"band1", "band2"};
  const TFloat64List flux = {1e-15, 1e-14};
  const TFloat64List fluxErr = {1e-18, 2e-18};
  std::shared_ptr<CPhotometricData> photoData =
      std::make_shared<CPhotometricData>(name, flux, fluxErr);

  spc.SetPhotData(photoData);
  std::shared_ptr<const CPhotometricData> photoData_out = spc.GetPhotData();
  BOOST_CHECK(photoData_out->GetFlux("band1") == 1e-15);
  BOOST_CHECK(photoData_out->GetFlux("band2") == 1e-14);

  // SetErrorAxis (cp)
  spc.SetErrorAxis(noiseAxis);
  BOOST_CHECK(spc.GetFluxAxis().GetError().GetSamplesVector() ==
              noiseAxis.GetSamplesVector());

  // SetErrorAxis (mv)
  noiseAxis[2] = 0.0;
  spc.SetErrorAxis(std::move(noiseAxis));
  BOOST_CHECK(spc.GetFluxAxis().GetError().GetSamplesVector()[2] == 0.0);

  bool isNoiseEmpty = spc.IsNoiseEmpty();
  BOOST_CHECK(isNoiseEmpty == false);

  // SetSpectralAxis (cp)
  CSpectrumSpectralAxis spectralAxis_2;
  BOOST_CHECK_THROW(spc.SetSpectralAxis(spectralAxis_2), GlobalException);

  spectralAxis_2.SetSize(spcAxisSize);
  spc.SetSpectralAxis(spectralAxis_2);
  BOOST_CHECK(spc.GetSpectralAxis().GetSamplesVector() ==
              TFloat64List(spcAxisSize, 0.0));

  // SetSpectralAxis (mv)
  spectralAxis_2[0] = 1211;
  spc.SetSpectralAxis(std::move(spectralAxis_2));
  BOOST_CHECK(spc.GetSpectralAxis().GetSamplesVector()[0] == 1211);
  BOOST_CHECK(spectralAxis_2.GetSamplesCount() == 0);

  CSpectrumSpectralAxis spectralAxis_3;
  BOOST_CHECK_THROW(spc.SetSpectralAxis(std::move(spectralAxis_3)),
                    GlobalException);

  // SetFluxAxis (cp)
  reinitializeAxis();
  spc.SetSpectralAxis(spcAxis);
  CSpectrumFluxAxis fluxAxis_2;
  BOOST_CHECK_THROW(spc.SetFluxAxis(fluxAxis_2), GlobalException);

  fluxAxis_2.SetSize(spcAxisSize);
  spc.SetFluxAxis(fluxAxis_2);
  BOOST_CHECK(spc.GetFluxAxis().GetSamplesVector() ==
              TFloat64List(spcAxisSize, 0.0));

  // SetFluxAxis (mv)
  fluxAxis_2[0] = 1.2;
  spc.SetFluxAxis(std::move(fluxAxis_2));
  BOOST_CHECK(spc.GetFluxAxis().GetSamplesVector()[0] == 1.2);

  CSpectrumFluxAxis fluxAxis_3;
  BOOST_CHECK_THROW(spc.SetFluxAxis(std::move(fluxAxis_3)), GlobalException);

  // SetSpectralAndFluxAxes
  BOOST_TEST_MESSAGE("=======");
  BOOST_TEST_MESSAGE("SetSpectralAndFluxAxes");
  BOOST_TEST_MESSAGE("=======");
  spectralAxis_2 = spc.GetSpectralAxis();
  fluxAxis_2 = spc.GetFluxAxis();
  BOOST_CHECK_THROW(spc.SetSpectralAndFluxAxes(spectralAxis_3, fluxAxis_3),
                    GlobalException);
  BOOST_CHECK_THROW(spc.SetSpectralAndFluxAxes(spectralAxis_2, fluxAxis_3),
                    GlobalException);

  spectralAxis_2[0] = 1210.;
  fluxAxis_2[0] = 1.5;
  spc.SetSpectralAndFluxAxes(spectralAxis_2, fluxAxis_2);
  BOOST_CHECK(spc.GetSpectralAxis().GetSamplesVector()[0] == 1210.);
  BOOST_CHECK(spc.GetFluxAxis().GetSamplesVector()[0] == 1.5);

  bool isEmpty = spc.IsFluxEmpty();
  BOOST_CHECK(isEmpty == false);

  // SetType
  reinitializeAxis();
  spc.SetSpectralAndFluxAxes(spcAxis, fluxAxis);
  spc.InitSpectrumContinuum(*paramStore);
  spc.SetContinuumEstimationMethod("raw");
  spc.EstimateContinuum();
  CSpectrumFluxAxis flux_out;
  spc.SetType(CSpectrum::EType::nType_continuumOnly);
  BOOST_CHECK(spc.GetType() == 2);
  flux_out = spc.GetFluxAxis();
  BOOST_CHECK(flux_out.GetSamplesVector() ==
              spc.GetContinuumFluxAxis().GetSamplesVector());
  spc.ResetContinuum();
  flux_out = spc.GetFluxAxis_();
  BOOST_CHECK(flux_out.GetSamplesVector() ==
              spc.GetContinuumFluxAxis_().GetSamplesVector());

  spc.SetType(CSpectrum::EType::nType_noContinuum);
  BOOST_CHECK(spc.GetType() == 3);
  spc.ResetContinuum();
  flux_out = spc.GetFluxAxis();
  BOOST_CHECK(flux_out.GetSamplesVector() ==
              spc.GetWithoutContinuumFluxAxis().GetSamplesVector());
  spc.ResetContinuum();
  flux_out = spc.GetFluxAxis_();
  BOOST_CHECK(flux_out.GetSamplesVector() ==
              spc.GetWithoutContinuumFluxAxis_().GetSamplesVector());

  spc.SetType(CSpectrum::EType::nType_raw);
  BOOST_CHECK(spc.GetType() == 1);
  flux_out = spc.GetFluxAxis();
  BOOST_CHECK(flux_out.GetSamplesVector() ==
              spc.GetRawFluxAxis().GetSamplesVector());
  flux_out = spc.GetFluxAxis_();
  BOOST_CHECK(flux_out.GetSamplesVector() ==
              spc.GetRawFluxAxis_().GetSamplesVector());
}

BOOST_AUTO_TEST_CASE(continuum_test) {
  // InitSpectrum
  TFloat64List spectralList_2 = {0, 1, 2};
  TFloat64List fluxList_2 = {1, 2, 1};
  CSpectrum spc_2(spectralList_2, fluxList_2);
  TAxisSampleList rawFlux2 = spc_2.GetRawFluxAxis_().GetSamplesVector();
  rawFlux2.pop_back();
  spc_2.GetRawFluxAxis_().setSamplesVector(rawFlux2);
  BOOST_CHECK_THROW(spc_2.InitSpectrumContinuum(*paramStore), GlobalException);
  spc.InitSpectrumContinuum(*paramStore);
  BOOST_CHECK(spc.GetContinuumEstimationMethod() == "IrregularSamplingMedian");
  BOOST_CHECK(spc.GetMedianWinsize() == 74.);
  BOOST_CHECK(spc.GetMedianEvenReflection() == false);

  // Estimate continuum
  spc.SetMedianEvenReflection(true);
  spc.EstimateContinuum();
  TFloat64List continuum = spc.GetContinuumFluxAxis().GetSamplesVector();
  TFloat64List rawFlux = spc.GetRawFluxAxis().GetSamplesVector();
  TFloat64List fluwWithoutContinuum(rawFlux.size(), 0.0);
  for (Int32 i = 0; i < rawFlux.size(); i++) {
    fluwWithoutContinuum[i] = rawFlux[i] - continuum[i];
  }
  BOOST_CHECK(spc.GetWithoutContinuumFluxAxis().GetSamplesVector() ==
              fluwWithoutContinuum);

  // Invert
  bool result = spc.InvertFlux();
  for (Int32 i = 0; i < rawFlux.size(); i++) {
    BOOST_CHECK(spc.GetRawFluxAxis().GetSamplesVector()[i] == -rawFlux[i]);
    BOOST_CHECK(spc.GetContinuumFluxAxis().GetSamplesVector()[i] ==
                -continuum[i]);
    BOOST_CHECK(spc.GetWithoutContinuumFluxAxis().GetSamplesVector()[i] ==
                -fluwWithoutContinuum[i]);
  }
  result = spc.InvertFlux();

  // ApplyAmplitude
  spc.ApplyAmplitude(2.);
  for (Int32 i = 0; i < rawFlux.size(); i++) {
    BOOST_CHECK(spc.GetRawFluxAxis().GetSamplesVector()[i] == 2 * rawFlux[i]);
    BOOST_CHECK(spc.GetContinuumFluxAxis().GetSamplesVector()[i] ==
                2 * continuum[i]);
    BOOST_CHECK(spc.GetWithoutContinuumFluxAxis().GetSamplesVector()[i] ==
                2 * fluwWithoutContinuum[i]);
  }

  // ValidateSpectrum
  spc.ApplyAmplitude(0.5);
  TFloat64Range lambdaRange(spectralList[0], spectralList[spcAxisSize - 1]);
  BOOST_CHECK_NO_THROW(spc.ValidateSpectrum(lambdaRange, true));
  // not IsValid
  TAxisSampleList rawFlux3 = spc.GetRawFluxAxis_().GetSamplesVector();
  rawFlux3.push_back(5592.);
  spc.GetRawFluxAxis_().setSamplesVector(rawFlux3);
  BOOST_CHECK_THROW(spc.ValidateSpectrum(lambdaRange, true), GlobalException);
  rawFlux3.pop_back();
  spc.GetRawFluxAxis_().setSamplesVector(rawFlux3);
  // correct flux
  rawFlux3[1] = fluxList[1];
  spc.GetRawFluxAxis_().setSamplesVector(rawFlux3);
  spc.ValidateSpectrum(lambdaRange, true);
  BOOST_CHECK(spc.GetRawFluxAxis_().GetSamplesVector()[1] == fluxList[1]);
  // not ValidateFlux
  rawFlux3[1] = std::numeric_limits<double>::infinity();
  spc.GetRawFluxAxis_().setSamplesVector(rawFlux3);
  BOOST_CHECK_THROW(spc.ValidateSpectrum(lambdaRange, false), GlobalException);
  rawFlux3[1] = fluxList[1];
  spc.GetRawFluxAxis_().setSamplesVector(rawFlux3);
  // not ValidateNoise
  TFloat64List error = spc.GetRawFluxAxis_().GetError().GetSamplesVector();
  error[1] = std::numeric_limits<double>::infinity();
  spc.GetRawFluxAxis_().setError(CSpectrumNoiseAxis(error));
  BOOST_CHECK_THROW(spc.ValidateSpectrum(lambdaRange, false), GlobalException);
  error[1] = noiseList[1];
  spc.GetRawFluxAxis_().setError(CSpectrumNoiseAxis(error));
  // LSF spectralAxis don't cover lambdaRange
  spc.SetLSF(LSF);
  TFloat64Range lambdaRange2(4680, 4712);
  BOOST_CHECK_THROW(spc.ValidateSpectrum(lambdaRange2, false), GlobalException);

  // SetContinuumEstimationMethod
  spc.SetContinuumEstimationMethod("raw");
  spc.EstimateContinuum();
  for (Int32 i = 0; i < rawFlux.size(); i++) {
    BOOST_CHECK(spc.GetRawFluxAxis().GetSamplesVector()[i] == rawFlux[i]);
    BOOST_CHECK(spc.GetContinuumFluxAxis().GetSamplesVector()[i] == rawFlux[i]);
    BOOST_CHECK(spc.GetWithoutContinuumFluxAxis().GetSamplesVector()[i] == 0.0);
  }

  spc.SetContinuumEstimationMethod("zero");
  spc.EstimateContinuum();
  for (Int32 i = 0; i < rawFlux.size(); i++) {
    BOOST_CHECK(spc.GetRawFluxAxis().GetSamplesVector()[i] == rawFlux[i]);
    BOOST_CHECK(spc.GetContinuumFluxAxis().GetSamplesVector()[i] == 0.);
    BOOST_CHECK(spc.GetWithoutContinuumFluxAxis().GetSamplesVector()[i] ==
                rawFlux[i]);
  }

  CSpectrum spc_3(spc);
  spc.SetContinuumEstimationMethod(spc_3.GetContinuumFluxAxis());
  spc.EstimateContinuum();
  for (Int32 i = 0; i < rawFlux.size(); i++) {
    BOOST_CHECK(spc.GetRawFluxAxis().GetSamplesVector()[i] == rawFlux[i]);
    BOOST_CHECK(spc.GetContinuumFluxAxis().GetSamplesVector()[i] == 0.);
    BOOST_CHECK(spc.GetWithoutContinuumFluxAxis().GetSamplesVector()[i] ==
                rawFlux[i]);
  }

  TAxisSampleList continuumFlux =
      spc_3.GetContinuumFluxAxis_().GetSamplesVector();
  continuumFlux.push_back(1.);
  spc_3.GetContinuumFluxAxis_().setSamplesVector(continuumFlux);
  spc_3.SetContinuumEstimationMethod("zero");
  BOOST_CHECK_THROW(
      spc.SetContinuumEstimationMethod(spc_3.GetContinuumFluxAxis()),
      GlobalException);

  spc.SetContinuumEstimationMethod("unkown");
  BOOST_CHECK_THROW(spc.EstimateContinuum(), GlobalException);

  // SetName
  const char *test_name = "spectrum_1";
  spc.SetName(test_name);
  BOOST_CHECK(spc.GetName() == test_name);

  // SetFullPath
  const char *test_path = "chemin";
  spc.SetFullPath(test_path);
  BOOST_CHECK(spc.GetFullPath() == test_path);
}

BOOST_AUTO_TEST_CASE(Calcul) {
  // GetMeanAndStdFluxInRange
  TFloat64Range range1(spectralList[0], spectralList[spcAxisSize - 1]);
  TFloat64Range range2(spectralList[0] - 5., spectralList[spcAxisSize - 1] + 5);
  TFloat64Range range3(spectralList[0] + 5., spectralList[spcAxisSize - 1] + 5);
  Float64 mean = 0.0;
  Float64 std = 0.0;

  bool result = spc.GetMeanAndStdFluxInRange(range1, mean, std);
  BOOST_CHECK(result == true);
  Float64 mean_out, std_out;
  CMask mask;
  CSpectrumSpectralAxis spectralAxis = spc.GetSpectralAxis();
  spectralAxis.GetMask(range1, mask);
  CSpectrumFluxAxis fluxAxis = spc.GetFluxAxis();
  result = fluxAxis.ComputeMeanAndSDev(mask, mean_out, std_out);
  BOOST_CHECK_CLOSE(mean, mean_out, 1e-12);
  BOOST_CHECK_CLOSE(std, std_out, 1e-12);

  result = spc.GetMeanAndStdFluxInRange(range2, mean, std);
  BOOST_CHECK(result == false);

  result = spc.GetMeanAndStdFluxInRange(range3, mean, std);
  BOOST_CHECK(result == false);

  // GetLinearRegInRange
  TFloat64Range range4(spectralList[0], spectralList[1]);
  Float64 a, b;
  result = spc.GetLinearRegInRange(range4, a, b);
  BOOST_CHECK(result == true);
  BOOST_CHECK_CLOSE(
      a, (fluxAxis[1] - fluxAxis[0]) / (spectralAxis[1] - spectralAxis[0]),
      1e-12);
  BOOST_CHECK_CLOSE(b, fluxAxis[1] - a * spectralAxis[1], 1e-12);

  result = spc.GetLinearRegInRange(range2, a, b);
  BOOST_CHECK(result == false);

  result = spc.GetLinearRegInRange(range3, a, b);
  BOOST_CHECK(result == false);

  //-----------
  CSpectrum object_CSpectrum;
  BOOST_TEST_MESSAGE("index:" << object_CSpectrum.GetSampleCount());
  BOOST_CHECK(object_CSpectrum.GetSampleCount() == 0);

  int nbmin = 0;
  int nbmax = 12;

  CSpectrumFluxAxis m_FluxAxis(nbmax);
  CSpectrumSpectralAxis m_SpectralAxis(nbmax, false);

  TFloat64List error;
  error.resize(nbmax);
  for (int i = nbmin; i < nbmax; ++i) {
    m_SpectralAxis[i] = i + 1;

    if (i < 5) {
      m_FluxAxis[i] = 0.0;
      error[i] = 0.0;
    } else if (i == 7) {
      m_FluxAxis[i] = std::nan("1");
      error[i] = std::nan("2");
    } else if (i == 9) {
      m_FluxAxis[i] = std::numeric_limits<double>::infinity();
      error[i] = std::numeric_limits<double>::infinity();
    } else {
      m_FluxAxis[i] = i + 2;
      error[i] = 1e-12;
    }

    BOOST_TEST_MESSAGE("m_SpectralAxis[i]:" << as_const(m_SpectralAxis)[i]);
  }
  m_FluxAxis.setError(CSpectrumNoiseAxis(error));

  BOOST_TEST_MESSAGE("index1:" << m_SpectralAxis.GetSamplesCount());
  BOOST_TEST_MESSAGE("index2:" << m_FluxAxis.GetSamplesCount());

  object_CSpectrum.SetSpectralAndFluxAxes(m_SpectralAxis, m_FluxAxis);

  const CSpectrumFluxAxis &const_FluxAxis = m_FluxAxis;
  const CSpectrumSpectralAxis &const_SpectralAxis = m_SpectralAxis;
  const CSpectrumNoiseAxis &const_noiseAxis = m_FluxAxis.GetError();

  //--------------------//

  CSpectrumFluxAxis _FluxAxis2(nbmax);
  CSpectrumFluxAxis _FluxAxis3(nbmax);
  CSpectrumFluxAxis _FluxAxis4(nbmax);
  CSpectrumFluxAxis _FluxAxis5(nbmax);
  CSpectrumFluxAxis _FluxAxis6(nbmax);
  CSpectrumFluxAxis _FluxAxis7(nbmax);

  TFloat64List error2(nbmax, 0.0);
  TFloat64List error3(nbmax, 0.0);
  TFloat64List error6(nbmax, 0.0);

  for (int i = nbmin; i < nbmax; ++i) {
    _FluxAxis2[i] = (i + 2) * 1e+3;
    _FluxAxis3[i] = _FluxAxis2[i];
    _FluxAxis4[i] = 0.0;
    _FluxAxis6[i] = _FluxAxis2[i];
    _FluxAxis7[i] = _FluxAxis2[i];
    error2[i] = 1e-5;
    error3[i] = 0.0;

    if (i < 5) {
      _FluxAxis5[i] = std::nan("5");
      error6[i] = std::nan("6");
    } else if (i == 5) {
      _FluxAxis5[i] = 1e+3;
      error6[i] = 1e-9;
    } else {
      _FluxAxis5[i] = std::numeric_limits<double>::infinity();
      error6[i] = std::numeric_limits<double>::infinity();
    }
  }

  _FluxAxis2.setError(CSpectrumNoiseAxis(error2));
  _FluxAxis3.setError(CSpectrumNoiseAxis(error3));
  _FluxAxis4.setError(CSpectrumNoiseAxis(error2));
  _FluxAxis5.setError(CSpectrumNoiseAxis(error2));
  _FluxAxis6.setError(CSpectrumNoiseAxis(error6));

  CSpectrum object_CSpectrum4(m_SpectralAxis, _FluxAxis2);
  CSpectrum object_CSpectrum5(m_SpectralAxis, _FluxAxis3);
  CSpectrum object_CSpectrum6(m_SpectralAxis, _FluxAxis4);
  CSpectrum object_CSpectrum7(m_SpectralAxis, _FluxAxis5);
  CSpectrum object_CSpectrum8(m_SpectralAxis, _FluxAxis6);
  CSpectrum object_CSpectrum2b(m_SpectralAxis, _FluxAxis7);

  //--------------------//
  // test ValidateSpectralAxis //
  BOOST_CHECK_THROW(object_CSpectrum.ValidateSpectralAxis(11, 14.1),
                    GlobalException); // cas où l'intervalle est à l'extérieur
  BOOST_CHECK_THROW(
      object_CSpectrum.ValidateSpectralAxis(10, 14.1),
      GlobalException); // cas où borne sup de l'intervalle est à l'extérieur
  BOOST_CHECK_THROW(
      object_CSpectrum.ValidateSpectralAxis(0, 10) == false,
      GlobalException); // cas où borne inf de l'intervalle est à l'extérieur

  //--------------------//
  // test ValidateFlux
  BOOST_CHECK_THROW(object_CSpectrum.ValidateFlux(1, 11.1),
                    GlobalException); // cas dans tout l'intervalle
  BOOST_CHECK_THROW(object_CSpectrum.ValidateFlux(1, 5.1),
                    GlobalException); // cas dans l'intervalle 1 à 5 avec 0.0
  BOOST_CHECK_NO_THROW(
      object_CSpectrum.ValidateFlux(1, 6.1)); // cas dans l'intervalle 1 à 6
  BOOST_CHECK_THROW(
      object_CSpectrum.ValidateFlux(6, 11.1),
      GlobalException); // cas dans l'intervalle 6 à 11 avec nan et inf
  BOOST_CHECK_NO_THROW(
      object_CSpectrum.ValidateFlux(6, 7.1)); // cas dans l'intervalle 6 à 7
  BOOST_CHECK_THROW(object_CSpectrum.ValidateFlux(7, 9.1),
                    GlobalException); // cas dans l'intervalle 7 à 9 avec nan
  BOOST_CHECK_NO_THROW(object_CSpectrum.ValidateFlux(
      9, 9.1)); // cas où l'intervalle est un point
  BOOST_CHECK_THROW(object_CSpectrum.ValidateFlux(9, 11.1),
                    GlobalException); // cas dans l'intervalle 9 à 11 avec inf
  BOOST_CHECK_THROW(object_CSpectrum.ValidateFlux(10, 10.1),
                    GlobalException); // cas où l'intervalle est un point inf
  TAxisSampleList fluxAxis2 =
      object_CSpectrum.GetFluxAxis_().GetSamplesVector();
  fluxAxis2.pop_back();
  object_CSpectrum.GetFluxAxis_().setSamplesVector(fluxAxis2);
  BOOST_CHECK_THROW(object_CSpectrum.ValidateFlux(1, 11.1),
                    GlobalException); // cas où les tailles sont différentes

  fluxAxis2.push_back(13.);
  object_CSpectrum.GetFluxAxis_().setSamplesVector(fluxAxis2);
  CSpectrumNoiseAxis err = object_CSpectrum.GetFluxAxis_().GetError();
  TAxisSampleList errAxisList = err.GetSamplesVector();
  errAxisList.push_back(1e-12);
  err.setSamplesVector(errAxisList);
  BOOST_CHECK_THROW(object_CSpectrum.GetFluxAxis_().setError(err),
                    GlobalException);

  //--------------------//
  // test ValidateNoise

  BOOST_CHECK_THROW(object_CSpectrum.ValidateNoise(1, 11.1),
                    GlobalException); // cas dans tout l'intervalle
  BOOST_CHECK_THROW(object_CSpectrum.ValidateNoise(1, 5.1),
                    GlobalException); // cas dans l'intervalle 1 à 5 avec 0.0
  BOOST_CHECK_THROW(object_CSpectrum.ValidateNoise(1, 6.1),
                    GlobalException); // cas dans l'intervalle 1 à 6
  BOOST_CHECK_THROW(
      object_CSpectrum.ValidateNoise(6, 11.1),
      GlobalException); // cas dans l'intervalle 6 à 11 avec nan et inf
  BOOST_CHECK_NO_THROW(
      object_CSpectrum.ValidateNoise(6, 7.1)); // cas dans l'intervalle 6 à 7
  BOOST_CHECK_THROW(object_CSpectrum.ValidateNoise(7, 9.1),
                    GlobalException); // cas dans l'intervalle 7 à 9 avec nan
  BOOST_CHECK_NO_THROW(object_CSpectrum.ValidateNoise(
      9, 9.1)); // cas où l'intervalle est un point
  BOOST_CHECK_THROW(object_CSpectrum.ValidateNoise(9, 11.1),
                    GlobalException); // cas dans l'intervalle 9 à 11 avec inf
  BOOST_CHECK_THROW(object_CSpectrum.ValidateNoise(10, 10.1),
                    GlobalException); // cas où l'intervalle est un point inf

  //--------------------//
  // test correctSpectrum

  // check throw : spectrum is not valid
  fluxAxis2.pop_back();
  object_CSpectrum.GetFluxAxis_().setSamplesVector(fluxAxis2);
  BOOST_CHECK_THROW(object_CSpectrum.correctSpectrum(1, 11.2), GlobalException);
  fluxAxis2.push_back(13.);
  object_CSpectrum.GetFluxAxis_().setSamplesVector(fluxAxis2);

  // cas où toutes les valeurs du flux et de l'erreur sont valides
  BOOST_CHECK(object_CSpectrum4.correctSpectrum(1, 11.2) == false);
  // cas où toutes les valeurs de l'erreur sont nulles, et le flux valide
  BOOST_CHECK_THROW(object_CSpectrum5.correctSpectrum(1, 11.2),
                    GlobalException);
  // cas où toutes les valeurs du flux sont nulles, et l'erreur valide
  BOOST_CHECK(object_CSpectrum6.correctSpectrum(1, 11.2) == false);

  // cas où une seule valeur du flux est valide, et l'erreur valide
  BOOST_CHECK(object_CSpectrum7.correctSpectrum(0.7, 11.3) == true);
  const TFloat64List f7 = object_CSpectrum7.GetFluxAxis().GetSamplesVector();
  const TFloat64List f7c{
      1e+2, 1e+2, 1e+2, 1e+2, 1e+2, 1e+3,
      1e+2, 1e+2, 1e+2, 1e+2, 1e+2, std::numeric_limits<double>::infinity()};
  BOOST_CHECK(f7 == f7c);
  const TFloat64List e7 =
      object_CSpectrum7.GetFluxAxis().GetError().GetSamplesVector();
  const TFloat64List e7c{1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-5,
                         1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-5};
  BOOST_CHECK(e7 == e7c);

  // cas où une seule valeur de l'erreur est valide, et le flux valide
  BOOST_CHECK(object_CSpectrum8.correctSpectrum(0.8, 11.3) == true);
  const TFloat64List f8 = object_CSpectrum8.GetFluxAxis().GetSamplesVector();
  const TFloat64List f8c = {7.0 * 1e+2, 7.0 * 1e+2, 7.0 * 1e+2, 7.0 * 1e+2,
                            7.0 * 1e+2, 7.0 * 1e+3, 7.0 * 1e+2, 7.0 * 1e+2,
                            7.0 * 1e+2, 7.0 * 1e+2, 7.0 * 1e+2, 13.0 * 1e+3};
  BOOST_CHECK(f8 == f8c);
  const TFloat64List e8 =
      object_CSpectrum8.GetFluxAxis().GetError().GetSamplesVector();
  const TFloat64List e8c{
      1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-9,
      1e-8, 1e-8, 1e-8, 1e-8, 1e-8, std::numeric_limits<double>::infinity()};
  BOOST_CHECK(e8 == e8c);

  // cas dans l'intervalle 1 à 5 avec 0.0
  CSpectrum object_CSpectrum2copy = object_CSpectrum;
  BOOST_CHECK_THROW(object_CSpectrum2copy.correctSpectrum(1, 5.4),
                    GlobalException);

  // cas dans l'intervalle 1 à 6
  CSpectrum object_CSpectrum9 = object_CSpectrum;
  BOOST_CHECK(object_CSpectrum9.correctSpectrum(1, 6.4) == true);
  // const access
  const CSpectrumFluxAxis &const_FluxAxis9 = object_CSpectrum9.GetFluxAxis();
  const CSpectrumNoiseAxis &const_noiseAxis9 =
      object_CSpectrum9.GetFluxAxis().GetError();

  for (int i = nbmin; i < nbmax; ++i) {
    if (i < 5) { // valeurs corrigées
      BOOST_CHECK(const_FluxAxis9[i] == 0.7);
      BOOST_CHECK(const_noiseAxis9[i] == 1e-11);
    } else if (i == 7) { // nan en position 8
      BOOST_CHECK(std::isnan(const_FluxAxis9[i]) == true);
      BOOST_CHECK(std::isnan(const_noiseAxis9[i]) == true);
    } else {
      BOOST_CHECK(const_FluxAxis9[i] == const_FluxAxis[i]);
      BOOST_CHECK(const_noiseAxis9[i] == const_noiseAxis[i]);
    }
  }

  // cas dans l'intervalle 6 à 7
  CSpectrum object_CSpectrum10 = object_CSpectrum;
  BOOST_CHECK(object_CSpectrum10.correctSpectrum(6, 7.4) == false);
  // const access
  const CSpectrumFluxAxis &const_FluxAxis10 = object_CSpectrum10.GetFluxAxis();
  const CSpectrumNoiseAxis &const_noiseAxis10 =
      object_CSpectrum10.GetFluxAxis().GetError();

  for (int i = nbmin; i < nbmax; ++i) {
    if (i == 7) { // nan en position 8
      BOOST_CHECK(std::isnan(const_FluxAxis10[i]) == true);
      BOOST_CHECK(std::isnan(const_noiseAxis10[i]) == true);
    } else {
      BOOST_CHECK(const_FluxAxis10[i] == const_FluxAxis[i]);
      BOOST_CHECK(const_noiseAxis10[i] == const_noiseAxis[i]);
    }
  }

  // cas dans l'intervalle 7 à 9 avec nan
  CSpectrum object_CSpectrum11 = object_CSpectrum;
  BOOST_CHECK(object_CSpectrum11.correctSpectrum(7, 9.4) == true);
  // const access
  const CSpectrumFluxAxis &const_FluxAxis11 = object_CSpectrum11.GetFluxAxis();
  const CSpectrumNoiseAxis &const_noiseAxis11 =
      object_CSpectrum11.GetFluxAxis().GetError();

  for (int i = nbmin; i < nbmax; ++i) {
    if (i == 7) { // nan en position 8 : valeur corrigée
      BOOST_CHECK(const_FluxAxis11[i] == 0.8);
      BOOST_CHECK(const_noiseAxis11[i] == 1e-11);
    } else {
      BOOST_CHECK(const_FluxAxis11[i] == const_FluxAxis[i]);
      BOOST_CHECK(const_noiseAxis11[i] == const_noiseAxis[i]);
    }
  }

  // cas où l'intervalle est un point
  CSpectrum object_CSpectrum12 = object_CSpectrum;
  BOOST_CHECK(object_CSpectrum12.correctSpectrum(9, 9.4) == false);
  // const access
  const CSpectrumFluxAxis &const_FluxAxis12 = object_CSpectrum12.GetFluxAxis();
  const CSpectrumNoiseAxis &const_noiseAxis12 =
      object_CSpectrum12.GetFluxAxis().GetError();

  for (int i = nbmin; i < nbmax; ++i) {
    if (i == 7) { // nan en position 8
      BOOST_CHECK(std::isnan(const_FluxAxis12[i]) == true);
      BOOST_CHECK(std::isnan(const_noiseAxis12[i]) == true);
    } else {
      BOOST_CHECK(const_FluxAxis12[i] == const_FluxAxis[i]);
      BOOST_CHECK(const_noiseAxis12[i] == const_noiseAxis[i]);
    }
  }

  // cas dans l'intervalle 9 à 11 avec inf
  CSpectrum object_CSpectrum13 = object_CSpectrum;
  BOOST_CHECK(object_CSpectrum13.correctSpectrum(9, 11.4) == true);
  // const access
  const CSpectrumFluxAxis &const_FluxAxis13 = object_CSpectrum13.GetFluxAxis();
  const CSpectrumNoiseAxis &const_noiseAxis13 =
      object_CSpectrum13.GetFluxAxis().GetError();

  for (int i = nbmin; i < nbmax; ++i) {
    if (i == 9) { // inf en position 10 : valeur corrigée
      BOOST_CHECK(const_FluxAxis13[i] == 1.0);
      BOOST_CHECK(const_noiseAxis13[i] == 1e-11);
    } else if (i == 7) { // nan en position 8
      BOOST_CHECK(std::isnan(const_FluxAxis13[i]) == true);
      BOOST_CHECK(std::isnan(const_noiseAxis13[i]) == true);
    } else {
      BOOST_CHECK(const_FluxAxis13[i] == const_FluxAxis[i]);
      BOOST_CHECK(const_noiseAxis13[i] == const_noiseAxis[i]);
    }
  }

  // cas après correction où l'intervalle est un point inf
  BOOST_CHECK(object_CSpectrum13.correctSpectrum(10, 10.5) == false);

  // cas après correction dans l'intervalle 6 à 9 avec nan
  BOOST_CHECK(object_CSpectrum11.correctSpectrum(6, 9.5) == false);

  // cas après correction dans l'intervalle 1 à 7
  BOOST_CHECK(object_CSpectrum9.correctSpectrum(1, 7.5) == false);

  // cas où l'intervalle est à l'extérieur
  BOOST_CHECK(object_CSpectrum.correctSpectrum(11, 14.5) == false);

  //--------------------//
  // test GetLambdaRange

  Float64 intervalle =
      object_CSpectrum.GetSpectralAxis().GetLambdaRange().GetLength();

  BOOST_CHECK_CLOSE(intervalle, object_CSpectrum.GetLambdaRange().GetLength(),
                    1e-12);

  BOOST_TEST_MESSAGE(
      "result9A:"
      << object_CSpectrum.GetSpectralAxis().GetLambdaRange().GetLength());
  BOOST_TEST_MESSAGE(
      "index55A:"
      << object_CSpectrum.GetSpectralAxis().GetLambdaRange().GetBegin());
  BOOST_TEST_MESSAGE(
      "index55B:"
      << object_CSpectrum.GetSpectralAxis().GetLambdaRange().GetEnd());

  //--------------------//
  // test GetMeanResolution

  Float64 result6 = object_CSpectrum.GetMeanResolution();
  Float64 result7 = object_CSpectrum.GetSpectralAxis().GetMeanResolution();

  BOOST_CHECK_CLOSE(result6, result7, 1e-12);

  BOOST_TEST_MESSAGE("result6:" << result6);
  BOOST_TEST_MESSAGE("result7:" << result7);

  //--------------------//
  // test GetResolution

  Float64 result8 = object_CSpectrum.GetResolution();
  Float64 result9 = object_CSpectrum.GetSpectralAxis().GetResolution();

  BOOST_CHECK_CLOSE(result8, result9, 1e-12);

  BOOST_TEST_MESSAGE("result8:" << result8);
  BOOST_TEST_MESSAGE("result9:" << result9);

  //--------------------//
  // //test removeContinuum

  CContinuumIrregularSamplingMedian remover2;

  for (int i = nbmin; i < nbmax; ++i) {
    m_FluxAxis[i] = 2.0;
  }
  object_CSpectrum.SetSpectralAndFluxAxes(m_SpectralAxis, m_FluxAxis);

  BOOST_CHECK(object_CSpectrum.RemoveContinuum(remover2) == true);
  BOOST_TEST_MESSAGE(
      "test Remove:" << object_CSpectrum.RemoveContinuum(remover2));
}

BOOST_AUTO_TEST_CASE(ExtractTest) {
  const CSpectrumSpectralAxis spcAxis({1., 2., 3., 4., 5.});
  const CSpectrumNoiseAxis noiseAxis({-1., -2., -3., -4., -5.});
  const CSpectrumFluxAxis fluxAxis(CSpectrumAxis({2., 4., 6., 8., 10.}),
                                   noiseAxis);
  const CSpectrum spc(spcAxis, fluxAxis);

  Int32 istart = 1;
  Int32 iend = 3;
  Int32 s = iend - istart + 1;

  const TFloat64List correctSpcAxis{2., 3., 4.};
  const TFloat64List correctFluxAxis{4., 6., 8};
  const TFloat64List correctNoiseAxis{-2., -3., -4.};

  ////////////////////////////
  // read results
  const CSpectrum extractedSpc = spc.extract(istart, iend);
  const TFloat64List &extractedSpcAxis =
      extractedSpc.GetSpectralAxis().GetSamplesVector();
  const TFloat64List &extractedFluxAxis =
      extractedSpc.GetFluxAxis().GetSamplesVector();
  const TFloat64List &extractedNoiseAxis =
      extractedSpc.GetFluxAxis().GetError().GetSamplesVector();

  // check results
  // check size is correct
  BOOST_CHECK(extractedSpcAxis.size() == s);
  BOOST_CHECK(extractedFluxAxis.size() == s);
  BOOST_CHECK(extractedNoiseAxis.size() == s);

  BOOST_CHECK_EQUAL_COLLECTIONS(extractedSpcAxis.begin(),
                                extractedSpcAxis.end(), correctSpcAxis.begin(),
                                correctSpcAxis.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(extractedFluxAxis.begin(),
                                extractedFluxAxis.end(),
                                correctFluxAxis.begin(), correctFluxAxis.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(
      extractedNoiseAxis.begin(), extractedNoiseAxis.end(),
      correctNoiseAxis.begin(), correctNoiseAxis.end());
}

BOOST_AUTO_TEST_CASE(rebin_test) {
  CSpectrum rebinedSpectrum;

  TFloat64Range range1(1213., 1218.);
  CSpectrumSpectralAxis tgtSpectralAxis_1({1212, 1213, 1213.5, 1214, 1214.5,
                                           1215, 1215.5, 1216, 1216.5, 1217,
                                           1217.5, 1218, 1219});
  CMask rebinedMask;

  std::string interp = "lin"; // lin, spline, precomputedfinegrid, ngp
  std::string errorRebinMethod = "rebin";
  spcLight.setRebinInterpMethod("lin");

  // check throw : spectrum is not valid
  TAxisSampleList fluxAxis2 = spcLight.GetFluxAxis_().GetSamplesVector();
  fluxAxis2.pop_back();
  spcLight.GetFluxAxis_().setSamplesVector(fluxAxis2);
  BOOST_CHECK_THROW(
      spcLight.Rebin(range1, tgtSpectralAxis_1, rebinedSpectrum, rebinedMask),
      GlobalException);
  fluxAxis2.push_back(6e-2);
  spcLight.GetFluxAxis_().setSamplesVector(fluxAxis2);

  // interp = "lin" et errorRebinMethod = "rebin"
  spcLight.Rebin(range1, tgtSpectralAxis_1, rebinedSpectrum, rebinedMask,
                 errorRebinMethod);
  TFloat64List rebinedFlux =
      rebinedSpectrum.GetRawFluxAxis().GetSamplesVector();
  Int32 N = rebinedFlux.size();
  BOOST_CHECK(N == 13);
  BOOST_CHECK(rebinedFlux[0] == 0. && rebinedFlux[N - 1] == 0.);
  for (Int32 i = 1; i < N - 1; i++) {
    BOOST_CHECK_CLOSE(rebinedFlux[i], 0.01 + (i - 1) * 5e-3, 1e-8);
  }
  TFloat64List rebinedError = rebinedSpectrum.GetErrorAxis().GetSamplesVector();
  BOOST_CHECK(rebinedError[0] == INFINITY && rebinedError[N - 1] == INFINITY);
  for (Int32 i = 1; i < N - 1; i++) {
    BOOST_CHECK(rebinedError[i] == 1e-4);
  }
  BOOST_CHECK(rebinedMask[0] == 0. && rebinedMask[N - 1] == 0.);
  for (Int32 i = 1; i < N - 1; i++) {
    BOOST_CHECK(rebinedMask[i] == 1);
  }

  // test setRebinType
  BOOST_CHECK_NO_THROW(spcLight.setRebinInterpMethod("lin"));
  BOOST_CHECK_NO_THROW(spcLight.setRebinInterpMethod("precomputedfinegrid"));
  BOOST_CHECK_NO_THROW(spcLight.setRebinInterpMethod("spline"));
  BOOST_CHECK_NO_THROW(spcLight.setRebinInterpMethod("ngp"));
  BOOST_CHECK_THROW(spcLight.setRebinInterpMethod("linn"), GlobalException);
}

BOOST_AUTO_TEST_SUITE_END()
