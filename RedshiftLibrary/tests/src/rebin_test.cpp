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
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"
#include "RedshiftLibrary/spectrum/rebin/rebin.h"
#include "RedshiftLibrary/spectrum/rebin/rebinFineGrid.h"
#include "RedshiftLibrary/spectrum/rebin/rebinLinear.h"
#include "RedshiftLibrary/spectrum/rebin/rebinNgp.h"
#include "RedshiftLibrary/spectrum/rebin/rebinSpline.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_FIXTURE_TEST_SUITE(Rebin, fixture_SpectrumLight)

BOOST_AUTO_TEST_CASE(rebin_test) {
  CSpectrum spc;
  std::unique_ptr<CRebin> rebin =
      std::unique_ptr<CRebin>(new CRebinLinear(spc));

  // test convert
  BOOST_CHECK_NO_THROW(rebin = std::move(*rebin).convert("lin"));
  BOOST_CHECK(rebin->getType() == "lin");
  BOOST_CHECK_NO_THROW(rebin =
                           std::move(*rebin).convert("precomputedfinegrid"));
  BOOST_CHECK(rebin->getType() == "precomputedfinegrid");
  BOOST_CHECK_NO_THROW(rebin = std::move(*rebin).convert("spline"));
  BOOST_CHECK(rebin->getType() == "spline");
  BOOST_CHECK_NO_THROW(rebin = std::move(*rebin).convert("ngp"));
  BOOST_CHECK(rebin->getType() == "ngp");
  BOOST_CHECK_THROW(std::move(*rebin).convert("linn"), GlobalException);

  // test create
  BOOST_CHECK_NO_THROW(CRebin::create("lin", spc));
  BOOST_CHECK_NO_THROW(CRebin::create("precomputedfinegrid", spc));
  BOOST_CHECK_NO_THROW(CRebin::create("spline", spc));
  BOOST_CHECK_NO_THROW(CRebin::create("ngp", spc));
  BOOST_CHECK_THROW(CRebin::create("linn", spc), GlobalException);
}

BOOST_AUTO_TEST_CASE(rebinLinear_test) {
  CSpectrum rebinedSpectrum;

  TFloat64Range range1(10., 30.);
  CSpectrumSpectralAxis tgtSpectralAxis_1({9., 10., 15., 20., 25., 30., 31.});
  CMask rebinedMask;
  std::string interp = "lin"; // lin, spline, precomputedfinegrid, ngp
  std::string errorRebinMethod = "rebin";

  // check throw : range is not included in spectral axis
  std::unique_ptr<CRebin> rebin =
      std::unique_ptr<CRebin>(new CRebinLinear(spcLight));
  TFloat64Range range2(9., 11.);
  BOOST_CHECK_THROW(rebin->compute(range2, tgtSpectralAxis_1, rebinedSpectrum,
                                   rebinedMask, errorRebinMethod),
                    GlobalException);

  // interp = "lin" et errorRebinMethod = "rebin
  rebin->compute(range1, tgtSpectralAxis_1, rebinedSpectrum, rebinedMask,
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
  rebin->compute(range1, tgtSpectralAxis_2, rebinedSpectrum, rebinedMask,
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
  CSpectrumNoiseAxis noiseAxis = spcLight.GetErrorAxis();
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

  // interp = "lin" et errorRebinMethod = "no"
  errorRebinMethod = "no";
  rebin->compute(range1, tgtSpectralAxis_1, rebinedSpectrum, rebinedMask,
                 errorRebinMethod);
  rebinedFlux = rebinedSpectrum.GetRawFluxAxis().GetSamplesVector();
  N = rebinedFlux.size();
  BOOST_CHECK(N == 7);
  BOOST_CHECK(rebinedFlux[0] == 0. && rebinedFlux[N - 1] == 0.);
  for (Int32 i = 1; i < N - 1; i++) {
    BOOST_CHECK(rebinedFlux[i] == 1 + (i - 1) * 0.5);
  }
  rebinedError = rebinedSpectrum.GetErrorAxis().GetSamplesVector();
  TFloat64List errorRef(7, 1.);
  BOOST_CHECK(rebinedError[0] == errorRef[0]);
}

BOOST_AUTO_TEST_CASE(rebinFineGrid_test) {
  CSpectrum rebinedSpectrum;

  TFloat64Range range1(10., 30.);
  CSpectrumSpectralAxis tgtSpectralAxis_1({9., 10., 15., 20., 25., 30., 31.});
  CMask rebinedMask;
  std::string interp = "precomputedfinegrid";
  std::string errorRebinMethod = "no";

  // interp = "precomputedfinegrid" et errorRebinMethod = "no"
  std::unique_ptr<CRebinFineGrid> rebin =
      std::unique_ptr<CRebinFineGrid>(new CRebinFineGrid(spcLight));
  spcLight.setRebinInterpMethod("precomputedfinegrid");
  rebin->compute(range1, tgtSpectralAxis_1, rebinedSpectrum, rebinedMask,
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
  BOOST_CHECK_THROW(rebin->compute(range1, tgtSpectralAxis_1, rebinedSpectrum,
                                   rebinedMask, errorRebinMethod),
                    GlobalException);

  // check throw : bad RebinFineGrid
  rebin->m_pfgFlux = {};
  BOOST_CHECK_THROW(rebin->compute(range1, tgtSpectralAxis_1, rebinedSpectrum,
                                   rebinedMask, errorRebinMethod),
                    GlobalException);

  // RebinFinegrid
  rebin->rebinFineGrid();
  TFloat64List pgfFlux = rebin->m_pfgFlux;
  Int32 n = pgfFlux.size();
  for (Int32 i = 0; i < 201; i++) {
    BOOST_CHECK_CLOSE(pgfFlux[i], 1 + i * 0.01, 1e-12);
  }
  for (Int32 i = 201; i < n; i++) {
    BOOST_CHECK(pgfFlux[i] == 0);
  }

  // clearFineGrid
  rebin->reset();
  BOOST_CHECK(rebin->m_pfgFlux.size() == 0);
  BOOST_CHECK(rebin->m_FineGridInterpolated == false);
}

BOOST_AUTO_TEST_CASE(rebinSpline_test) {
  CSpectrum rebinedSpectrum;

  // create spectrum
  CSpectrumSpectralAxis spectralAxis(spectralList);
  CSpectrumNoiseAxis noiseAxis(noiseList);
  CSpectrumFluxAxis fluxAxis(fluxList);
  fluxAxis.GetError() = noiseAxis;
  CSpectrum rebinedSpectrum;
  CSpectrum spc(spectralAxis, fluxAxis);

  TFloat64Range range1(10., 30.);
  CSpectrumSpectralAxis tgtSpectralAxis_1({9., 10., 15., 20., 25., 30., 31.});
  CMask rebinedMask;
  std::string interp = "spline";
  std::string errorRebinMethod = "no";

  // interp = "spline" et errorRebinMethod = "no"
  std::unique_ptr<CRebin> rebin =
      std::unique_ptr<CRebin>(new CRebinSpline(spcLight));
  rebin->compute(range1, tgtSpectralAxis_1, rebinedSpectrum, rebinedMask,
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
  BOOST_CHECK_THROW(rebin->compute(range1, tgtSpectralAxis_1, rebinedSpectrum,
                                   rebinedMask, errorRebinMethod),
                    GlobalException);
}

BOOST_AUTO_TEST_CASE(rebinNgp_test) {
  CSpectrum rebinedSpectrum;

  TFloat64Range range1(10., 30.);
  CSpectrumSpectralAxis tgtSpectralAxis_1({9., 10., 15., 20., 25., 30., 31.});
  CMask rebinedMask;
  std::string interp = "ngp";
  std::string errorRebinMethod = "rebin";

  // interp = "ngp" et errorRebinMethod = "rebin"
  std::unique_ptr<CRebin> rebin =
      std::unique_ptr<CRebin>(new CRebinNgp(spcLight));
  rebin->compute(range1, tgtSpectralAxis_1, rebinedSpectrum, rebinedMask,
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
  rebin->compute(range1, tgtSpectralAxis_2, rebinedSpectrum, rebinedMask,
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

  // interp = "ngp" et errorRebinMethod = "no"
  errorRebinMethod = "no";
  rebin->compute(range1, tgtSpectralAxis_1, rebinedSpectrum, rebinedMask,
                 errorRebinMethod);
  rebinedFlux = rebinedSpectrum.GetRawFluxAxis().GetSamplesVector();
  N = rebinedFlux.size();
  BOOST_CHECK(N == 7);
  BOOST_CHECK(rebinedFlux[0] == 0. && rebinedFlux[N - 1] == 0.);
  for (Int32 i = 1; i < N - 1; i++) {
    Int32 i2 = i / 2;
    if (2 * i2 == i)
      BOOST_CHECK(rebinedFlux[i] == rebinedFlux[i + 1]);
    else
      BOOST_CHECK(rebinedFlux[i] == 1 + (i - 1) * 0.5);
  }
  rebinedError = rebinedSpectrum.GetErrorAxis().GetSamplesVector();
  TFloat64List errorRef(7, 1.);
  BOOST_CHECK(rebinedError == errorRef);
}

BOOST_AUTO_TEST_SUITE_END()