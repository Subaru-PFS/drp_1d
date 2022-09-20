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
// #include "RedshiftLibrary/continuum/irregularsamplingmedian.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"
#include "RedshiftLibrary/spectrum/rebin/rebin.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/tests/test-tools.h"
#include <boost/test/unit_test.hpp>

using namespace NSEpic;
using namespace CPFTest;

TFloat64List spectralList = {10, 20, 30};
TFloat64List fluxList = {1, 2, 3};
TFloat64List noiseList = {0.1, 0.1, 0.1};
TFloat64List maskList = {0, 1, 1};
TFloat64List widthList = {1, 1, 1};

const std::string jsonString =
    "{\"smoothWidth\" : 0.5,"
    "\"continuumRemoval\" : { \"medianKernelWidth\" : 74.0, "
    "\"medianEvenReflection\" : false, "
    "\"method\" : \"IrregularSamplingMedian\"}}";

class MyInputContext {
public:
  std::shared_ptr<CParameterStore> paramStore;
  std::shared_ptr<CLSF> LSF;

  void InitContext() {
    TScopeStack scopeStack;
    paramStore = std::make_shared<CParameterStore>(scopeStack);
    paramStore->FromString(jsonString);

    std::string lsfType = "GaussianVariableWidth";
    std::shared_ptr<TLSFArguments> args =
        std::make_shared<TLSFGaussianVarWidthArgs>(spectralList, widthList);
    LSF = LSFFactory.Create(lsfType, args);
  }

  std::shared_ptr<CParameterStore> GetParameterStore() { return paramStore; }
  std::shared_ptr<CLSF> GetLSF() { return LSF; }
};

BOOST_AUTO_TEST_SUITE(Rebin)

BOOST_AUTO_TEST_CASE(rebinLinear_test) {
  // Initialize context
  MyInputContext ctx;
  ctx.InitContext();

  // create spectrum
  CSpectrumSpectralAxis spectralAxis(spectralList);
  CSpectrumNoiseAxis noiseAxis(noiseList);
  CSpectrumFluxAxis fluxAxis(fluxList);
  fluxAxis.GetError() = noiseAxis;
  CSpectrum spc(spectralAxis, fluxAxis);

  TFloat64Range range1(10., 30.);
  CSpectrumSpectralAxis tgtSpectralAxis_1({9., 10., 15., 20., 25., 30., 31.});
  CMask rebinedMask;
  CSpectrum rebinedSpectrum;
  std::string interp = "lin"; // lin, spline, precomputedfinegrid, ngp
  std::string errorRebinMethod = "rebin";

  // check throw : range is not included in spectral axis
  spc.setRebinType("lin");
  TFloat64Range range2(9., 11.);
  BOOST_CHECK_THROW(spc.m_rebin->compute(range2, tgtSpectralAxis_1,
                                         rebinedSpectrum, rebinedMask,
                                         errorRebinMethod),
                    GlobalException);

  // interp = "lin" et errorRebinMethod = "rebin
  spc.m_rebin->compute(range1, tgtSpectralAxis_1, rebinedSpectrum, rebinedMask,
                       errorRebinMethod);
  TFloat64List rebinedFlux =
      rebinedSpectrum.GetRawFluxAxis().GetSamplesVector();
  Int32 N = rebinedFlux.size();
  BOOST_CHECK(N == 7);
  BOOST_CHECK(rebinedFlux[0] == 0. && rebinedFlux[N - 1] == 0.);
  for (Int32 i = 1; i < N - 1; i++) {
    BOOST_CHECK(rebinedFlux[i] == 1 + (i - 1) * 0.5);
  }
  TFloat64List rebinedError = rebinedSpectrum.GetErrorAxis().GetSamplesVector();
  BOOST_CHECK(rebinedError[0] == INFINITY && rebinedError[N - 1] == INFINITY);
  for (Int32 i = 1; i < N - 1; i++) {
    BOOST_CHECK(rebinedError[i] == 0.1);
  }
  BOOST_CHECK(rebinedMask[0] == 0. && rebinedMask[N - 1] == 0.);
  for (Int32 i = 1; i < N - 1; i++) {
    BOOST_CHECK(rebinedMask[i] == 1);
  }

  // interp = "lin" et errorRebinMethod = "rebinVariance"
  CSpectrumSpectralAxis tgtSpectralAxis_2({9., 10., 15., 20., 25., 30.});
  errorRebinMethod = "rebinVariance";
  spc.m_rebin->compute(range1, tgtSpectralAxis_2, rebinedSpectrum, rebinedMask,
                       errorRebinMethod);
  rebinedFlux = rebinedSpectrum.GetRawFluxAxis().GetSamplesVector();
  N = rebinedFlux.size();
  BOOST_CHECK(N == 6);
  BOOST_CHECK(rebinedFlux[0] == 0.);
  for (Int32 i = 1; i < N; i++) {
    BOOST_CHECK(rebinedFlux[i] == 1 + (i - 1) * 0.5);
  }
  rebinedError = rebinedSpectrum.GetErrorAxis().GetSamplesVector();
  TFloat64List t_list_ref = {0., 0.5, 0., 0.5, 0.};
  BOOST_CHECK(rebinedError[0] == INFINITY);
  for (Int32 i = 1; i < N; i++) {
    Float64 t = t_list_ref[i - 1];
    Float64 error_ref =
        sqrt(noiseAxis[0] * noiseAxis[0] * (1 - t) * (1 - t) +
             noiseAxis[1] * noiseAxis[1] * t *
                 t); // noise is constant here (no need to loop over)
    error_ref *= sqrt(2.);
    BOOST_CHECK(rebinedError[i] == error_ref);
  }
  BOOST_CHECK(rebinedMask[0] == 0.);
  for (Int32 i = 1; i < N; i++) {
    BOOST_CHECK(rebinedMask[i] == 1);
  }
}

BOOST_AUTO_TEST_CASE(rebinFineGrid_test) {
  // Initialize context
  MyInputContext ctx;
  ctx.InitContext();

  // create spectrum
  CSpectrumSpectralAxis spectralAxis(spectralList);
  CSpectrumNoiseAxis noiseAxis(noiseList);
  CSpectrumFluxAxis fluxAxis(fluxList);
  fluxAxis.GetError() = noiseAxis;
  CSpectrum spc(spectralAxis, fluxAxis);

  TFloat64Range range1(10., 30.);
  CSpectrumSpectralAxis tgtSpectralAxis_1({9., 10., 15., 20., 25., 30., 31.});
  CMask rebinedMask;
  CSpectrum rebinedSpectrum;
  std::string interp = "precomputedfinegrid";
  std::string errorRebinMethod = "no";

  // interp = "precomputedfinegrid" et errorRebinMethod = "no"
  spc.setRebinType("precomputedfinegrid");
  spc.m_rebin->compute(range1, tgtSpectralAxis_1, rebinedSpectrum, rebinedMask,
                       errorRebinMethod);
  TFloat64List rebinedFlux =
      rebinedSpectrum.GetRawFluxAxis().GetSamplesVector();
  Int32 N = rebinedFlux.size();
  BOOST_CHECK(N == 7);
  BOOST_CHECK(rebinedFlux[0] == 0. && rebinedFlux[N - 1] == 0.);
  for (Int32 i = 1; i < N - 1; i++) {
    BOOST_CHECK(rebinedFlux[i] == 1 + (i - 1) * 0.5);
  }
  BOOST_CHECK(rebinedMask[0] == 0. && rebinedMask[N - 1] == 0.);
  for (Int32 i = 1; i < N - 1; i++) {
    BOOST_CHECK(rebinedMask[i] == 1);
  }

  // interp = "precomputedfinegrid" et errorRebinMethod != "no"
  errorRebinMethod = "rebin";
  bool result = spc.Rebin(range1, tgtSpectralAxis_1, rebinedSpectrum,
                          rebinedMask, interp, errorRebinMethod);
  BOOST_CHECK(result == false);

  // check throw : bad RebinFineGrid
  spc.m_rebin->m_pfgFlux = {};
  result = spc.m_rebin->compute(range1, tgtSpectralAxis_1, rebinedSpectrum,
                                rebinedMask, errorRebinMethod);
  BOOST_CHECK(result == false);

  // RebinFinegrid
  spc.m_rebin->rebinFineGrid();
  TFloat64List pgfFlux = spc.m_rebin->m_pfgFlux;
  Int32 n = pgfFlux.size();
  for (Int32 i = 0; i < 201; i++) {
    BOOST_CHECK_CLOSE(pgfFlux[i], 1 + i * 0.01, 1e-12);
  }
  for (Int32 i = 201; i < n; i++) {
    BOOST_CHECK(pgfFlux[i] == 0);
  }

  // clearFineGrid
  spc.m_rebin->clearFineGrid();
  BOOST_CHECK(spc.m_rebin->m_pfgFlux.size() == 0);
  BOOST_CHECK(spc.m_rebin->m_FineGridInterpolated == false);
}

BOOST_AUTO_TEST_CASE(rebinSpline_test) {
  // Initialize context
  MyInputContext ctx;
  ctx.InitContext();

  // create spectrum
  CSpectrumSpectralAxis spectralAxis(spectralList);
  CSpectrumNoiseAxis noiseAxis(noiseList);
  CSpectrumFluxAxis fluxAxis(fluxList);
  fluxAxis.GetError() = noiseAxis;
  CSpectrum spc(spectralAxis, fluxAxis);

  TFloat64Range range1(10., 30.);
  CSpectrumSpectralAxis tgtSpectralAxis_1({9., 10., 15., 20., 25., 30., 31.});
  CMask rebinedMask;
  CSpectrum rebinedSpectrum;
  std::string interp = "spline";
  std::string errorRebinMethod = "no";

  // interp = "spline" et errorRebinMethod = "no"
  spc.setRebinType("spline");
  spc.m_rebin->compute(range1, tgtSpectralAxis_1, rebinedSpectrum, rebinedMask,
                       errorRebinMethod);
  TFloat64List rebinedFlux =
      rebinedSpectrum.GetRawFluxAxis().GetSamplesVector();
  Int32 N = rebinedFlux.size();
  BOOST_CHECK(N == 7);
  BOOST_CHECK(rebinedFlux[0] == 0. && rebinedFlux[N - 1] == 0.);
  for (Int32 i = 1; i < N - 1; i++) {
    BOOST_CHECK(rebinedFlux[i] == 1 + (i - 1) * 0.5);
  }
  BOOST_CHECK(rebinedMask[0] == 0. && rebinedMask[N - 1] == 0.);
  for (Int32 i = 1; i < N - 1; i++) {
    BOOST_CHECK(rebinedMask[i] == 1);
  }

  // interp = "spline" et errorRebinMethod != "no"
  errorRebinMethod = "rebin";
  bool result = spc.m_rebin->compute(range1, tgtSpectralAxis_1, rebinedSpectrum,
                                     rebinedMask, errorRebinMethod);
  BOOST_CHECK(result == false);
}

BOOST_AUTO_TEST_CASE(rebinNgp_test) {
  // Initialize context
  MyInputContext ctx;
  ctx.InitContext();

  // create spectrum
  CSpectrumSpectralAxis spectralAxis(spectralList);
  CSpectrumNoiseAxis noiseAxis(noiseList);
  CSpectrumFluxAxis fluxAxis(fluxList);
  fluxAxis.GetError() = noiseAxis;
  CSpectrum spc(spectralAxis, fluxAxis);

  TFloat64Range range1(10., 30.);
  CSpectrumSpectralAxis tgtSpectralAxis_1({9., 10., 15., 20., 25., 30., 31.});
  CMask rebinedMask;
  CSpectrum rebinedSpectrum;
  std::string interp = "ngp";
  std::string errorRebinMethod = "rebin";

  // interp = "ngp" et errorRebinMethod = "rebin"
  spc.setRebinType("ngp");
  spc.m_rebin->compute(range1, tgtSpectralAxis_1, rebinedSpectrum, rebinedMask,
                       errorRebinMethod);
  TFloat64List rebinedFlux =
      rebinedSpectrum.GetRawFluxAxis().GetSamplesVector();
  Int32 N = rebinedFlux.size();
  BOOST_CHECK(N == 7);
  BOOST_CHECK(rebinedFlux[0] == 0. && rebinedFlux[N - 1] == 0.);
  for (Int32 i = 1; i < N - 1; i++) {
    Int32 i2 = i / 2;
    if (2 * i2 == i)
      BOOST_CHECK(rebinedFlux[i] == rebinedFlux[i + 1]);
    else
      BOOST_CHECK(rebinedFlux[i] == 1 + (i - 1) * 0.5);
  }
  TFloat64List rebinedError = rebinedSpectrum.GetErrorAxis().GetSamplesVector();
  BOOST_CHECK(rebinedError[0] == INFINITY && rebinedError[N - 1] == INFINITY);
  for (Int32 i = 1; i < N - 1; i++) {
    BOOST_CHECK(rebinedError[i] == 0.1);
  }
  BOOST_CHECK(rebinedMask[0] == 0. && rebinedMask[N - 1] == 0.);
  for (Int32 i = 1; i < N - 1; i++) {
    BOOST_CHECK(rebinedMask[i] == 1);
  }

  // interp = "ngp" et errorRebinMethod = "rebinVariance"
  CSpectrumSpectralAxis tgtSpectralAxis_2({9., 10., 15., 20., 25., 30.});
  errorRebinMethod = "rebinVariance";
  spc.Rebin(range1, tgtSpectralAxis_2, rebinedSpectrum, rebinedMask, interp,
            errorRebinMethod);
  rebinedFlux = rebinedSpectrum.GetRawFluxAxis().GetSamplesVector();
  N = rebinedFlux.size();
  BOOST_CHECK(N == 6);
  BOOST_CHECK(rebinedFlux[0] == 0.);
  for (Int32 i = 1; i < N; i++) {
    Int32 i2 = i / 2;
    if (2 * i2 == i)
      BOOST_CHECK(rebinedFlux[i] == rebinedFlux[i + 1]);
    else
      BOOST_CHECK(rebinedFlux[i] == 1 + (i - 1) * 0.5);
  }
  rebinedError = rebinedSpectrum.GetErrorAxis().GetSamplesVector();
  BOOST_CHECK(rebinedError[0] == INFINITY);
  for (Int32 i = 1; i < N; i++) {
    BOOST_CHECK(rebinedError[i] == 0.1 * sqrt(2.));
  }
  BOOST_CHECK(rebinedMask[0] == 0.);
  for (Int32 i = 1; i < N; i++) {
    BOOST_CHECK(rebinedMask[i] == 1);
  }
}

BOOST_AUTO_TEST_SUITE_END()