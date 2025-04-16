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
#include <boost/test/unit_test.hpp>

#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/continuum/irregularsamplingmedian.h"
#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

using namespace NSEpic;
using namespace std;

void print_flux(TAxisSampleList sample) {
  BOOST_TEST_MESSAGE("=======");
  for (Int32 i = 0; i < ssize(sample); i++) {
    BOOST_TEST_MESSAGE(sample[i]);
  }
  BOOST_TEST_MESSAGE("=======");
}

BOOST_AUTO_TEST_SUITE(continuum_test)

BOOST_AUTO_TEST_CASE(median_test) {
  CContinuumIrregularSamplingMedian sample;
  TFloat64List y_in = {1., 3., 5., 7., 9., 11., 13., 15., 17., 19., 21.};
  TFloat64List y_out;
  Int32 n_points;

  BOOST_TEST_MESSAGE("TEST MEDIAN SMOOTH");

  // odd n_points     : rest part is on right side
  // even n_points    : no rest part

  // n_points = 0
  n_points = 0;
  y_out = sample.MedianSmooth(y_in, n_points);
  BOOST_TEST_MESSAGE("n_points = 0");
  BOOST_CHECK(y_out == y_in);
  print_flux(y_out);

  // n_points = 1 > y_out == y_in
  n_points = 1;
  y_out = sample.MedianSmooth(y_in, n_points);
  BOOST_TEST_MESSAGE("n_points = 1");
  BOOST_CHECK(y_out == y_in);
  print_flux(y_out);

  // n_points = 3 (first and last range are truncated) :
  //           y_in       0   1   2   3   4   5   6   7   8   9   10
  //                      |   |   |   |   |   |   |   |   |   |   |
  // points for i=0   |___|___|...|
  // points for i=1       |___|___|...|
  // points for i=8                                   |___|___|...|
  // points for i=9                                       |___|___|...|
  // points for i=10                                          |___|___|...|
  n_points = 3;
  y_out = sample.MedianSmooth(y_in, n_points);
  BOOST_TEST_MESSAGE("n_points = 3");
  // i = 0 : n_points = 2 in CMedian::Find
  BOOST_CHECK(y_out[0] == (y_in[0] + y_in[1]) / (n_points - 1));
  // i in [1,9] : n_points =3 in CMedian::Find
  for (Int32 i = 1; i < ssize(y_out) - 1; i++) {
    BOOST_CHECK(y_out[i] == y_in[i]);
  }
  // i = 10 : n_points = 2 in CMedian::Find
  BOOST_CHECK(y_out.back() ==
              (y_in[y_in.size() - 2] + y_in[y_in.size() - 1]) / (n_points - 1));
  print_flux(y_out);

  // max range : 2*n_points - 1
  //----------------------------

  // n_points = 21 :
  // final n_points = 11 -> odd median : y_out = y_in[5]
  TFloat64List res(11, 11.);
  n_points = 21;
  y_out = sample.MedianSmooth(y_in, n_points);
  BOOST_TEST_MESSAGE("n_points = 21");
  BOOST_CHECK(y_out == res);
  print_flux(y_out);

  // even y_in test
  TFloat64List y_in_2 = {1., 3., 5., 7., 9., 11., 13., 15., 17., 19.};
  TFloat64List res_2(10, 10.);
  TFloat64List y_out_2;

  // n_points = 19 :
  // final n_points = 19 -> even  median : y_out = 0.5 * (y_in[4] + y_in[5])
  n_points = 19;
  y_out_2 = sample.MedianSmooth(y_in_2, n_points);
  BOOST_TEST_MESSAGE("n_points = 19");
  BOOST_CHECK(y_out_2 == res_2);
  print_flux(y_out_2);

  // unordered vector
  //----------------------------
  y_in = {1., 5., 3., 2., 4., 6., 8., 1., 2., 1., 0.};
  n_points = 3;
  y_out = sample.MedianSmooth(y_in, n_points);
  BOOST_TEST_MESSAGE("n_points = 3");
  // i = 0 : n_points = 2 in CMedian::Find
  BOOST_CHECK(y_out[0] == (y_in[0] + y_in[1]) / (n_points - 1));
  // i in [1,9] : n_points = 3 in CMedian::Find
  TFloat64List y_in_sorted(3);
  for (Int32 i = 1; i < ssize(y_out) - 1; i++) {
    std::copy(y_in.begin() + (i - 1), y_in.begin() + (i + 2),
              y_in_sorted.begin());
    std::sort(y_in_sorted.begin(), y_in_sorted.end());
    BOOST_CHECK(y_out[i] == y_in_sorted[1]);
  }
  // i = 10 : n_points = 2 in CMedian::Find
  BOOST_CHECK(y_out.back() ==
              (y_in[y_in.size() - 2] + y_in[y_in.size() - 1]) / (n_points - 1));
  print_flux(y_out);
}

BOOST_AUTO_TEST_CASE(mean_test) {
  CContinuumIrregularSamplingMedian sample;
  TFloat64List y_in = {1., 3., 5., 7., 9., 11., 13., 15., 17., 19., 21.};
  TFloat64List y_out;
  Int32 n_range;

  BOOST_TEST_MESSAGE("TEST MEAN SMOOTH");

  // odd n_range  : rest part is on right side
  // even n_range : no rest part

  // n_range null
  n_range = 0;
  y_out = sample.MeanSmooth(y_in, n_range);
  BOOST_TEST_MESSAGE("n_range = 0");
  BOOST_CHECK(y_out == y_in);
  print_flux(y_out);

  // n_range = 1 > y_out == y_in
  n_range = 1;
  y_out = sample.MeanSmooth(y_in, n_range);
  BOOST_TEST_MESSAGE("n_range = 1");
  BOOST_CHECK(y_out == y_in);
  print_flux(y_out);

  // n_range = 2 (first range is truncated) :
  //           y_in       0   1   2   3   4   5   6   7   8   9   10
  //                      |   |   |   |   |   |   |   |   |   |   |
  // range for i=0    |___|...|
  // range for i=1        |___|...|
  // range for i=9                                        |___|...|
  // range for i=10                                           |___|...|
  n_range = 2;
  y_out = sample.MeanSmooth(y_in, n_range);
  BOOST_TEST_MESSAGE("n_range = 2");
  BOOST_CHECK(y_out.front() == y_in.front());
  for (Int32 i = 1; i < ssize(y_out); i++) {
    BOOST_CHECK(y_out[i] == (y_in[i] + y_in[i - 1]) / n_range);
  }
  print_flux(y_out);

  // max range : 2*n_range -1
  //----------------------------

  // n_range = 21 :
  TFloat64List res(11, 11.);
  n_range = 21;
  y_out = sample.MeanSmooth(y_in, n_range);
  BOOST_TEST_MESSAGE("n_range = 21");
  BOOST_CHECK(y_out == res);
  print_flux(y_out);

  // even y_in test
  TFloat64List y_in_2 = {1., 3., 5., 7., 9., 11., 13., 15., 17., 19.};
  TFloat64List res_2(10, 10.);
  TFloat64List y_out_2;

  // n_range = 19 :
  n_range = 19;
  y_out_2 = sample.MeanSmooth(y_in_2, n_range);
  BOOST_TEST_MESSAGE("n_range = 19");
  BOOST_CHECK(y_out_2 == res_2);
  print_flux(y_out_2);
}

BOOST_AUTO_TEST_CASE(evenMirror_test) {
  CContinuumIrregularSamplingMedian sample;
  TFloat64List y_in;
  TFloat64List y_out;
  Int32 beg;
  Int32 N;
  Int32 Nreflex;

  BOOST_TEST_MESSAGE("TEST EVEN MIRROR");

  // TEST Nreflex < N
  y_in = {1., 2., 3., 2., 1.};
  beg = 0;
  N = y_in.size();
  Nreflex = 2;
  y_out = sample.EvenMirror(y_in.begin(), y_in.end(), Nreflex);
  BOOST_CHECK(y_out.size() == 2 * Nreflex + y_in.size());
  TFloat64List y_rev(Nreflex);
  // first part
  copy(y_in.begin(), y_in.begin() + Nreflex, y_rev.begin());
  reverse(y_rev.begin(), y_rev.end());
  for (Int32 i = 0; i < Nreflex; i++) {
    BOOST_CHECK(y_out[i] == y_rev[i]);
  }
  // second part
  for (Int32 i = 0; i < N; i++) {
    BOOST_CHECK(y_out[i + Nreflex] == y_in[i]);
  }
  // third part
  copy(y_in.end() - Nreflex, y_in.end(), y_rev.begin());
  reverse(y_rev.begin(), y_rev.end());
  for (Int32 i = 0; i < Nreflex; i++) {
    BOOST_CHECK(y_out[i + N + Nreflex] == y_rev[i]);
  }
  print_flux(y_out);

  // TEST Nreflex = N
  y_in = {1., 2., 3., 2., 1.};
  beg = 0;
  N = y_in.size();
  Nreflex = y_in.size();
  y_out = sample.EvenMirror(y_in.begin(), y_in.end(), Nreflex);
  BOOST_CHECK(y_out.size() == 2 * Nreflex + y_in.size());
  y_rev.resize(Nreflex);
  // first part
  copy(y_in.begin(), y_in.begin() + Nreflex, y_rev.begin());
  reverse(y_rev.begin(), y_rev.end());
  for (Int32 i = 0; i < Nreflex; i++) {
    BOOST_CHECK(y_out[i] == y_rev[i]);
  }
  // second part
  for (Int32 i = 0; i < N; i++) {
    BOOST_CHECK(y_out[i + Nreflex] == y_in[i]);
  }
  // third part
  copy(y_in.end() - Nreflex, y_in.end(), y_rev.begin());
  reverse(y_rev.begin(), y_rev.end());
  for (Int32 i = 0; i < Nreflex; i++) {
    BOOST_CHECK(y_out[i + N + Nreflex] == y_rev[i]);
  }
  print_flux(y_out);
}

BOOST_AUTO_TEST_CASE(oddMirror_test) {
  CContinuumIrregularSamplingMedian sample;
  TFloat64List y_in;
  TFloat64List y_out;
  Int32 beg;
  Int32 N;
  Int32 Nreflex;

  BOOST_TEST_MESSAGE("TEST ODD MIRROR");

  // TEST Nreflex < N
  y_in = {1., 2., 3., 2., 1.};
  beg = 0;
  N = y_in.size();
  Nreflex = 2;
  y_out = sample.OddMirror(y_in.begin(), y_in.end(), Nreflex, y_in.front(),
                           y_in.back());
  BOOST_CHECK(y_out.size() == 2 * Nreflex + y_in.size());
  TFloat64List y_rev(Nreflex);
  // first part
  copy(y_in.begin(), y_in.begin() + Nreflex, y_rev.begin());
  reverse(y_rev.begin(), y_rev.end());
  for (Int32 i = 0; i < Nreflex; i++) {
    BOOST_CHECK(y_out[i] == 2 * y_in.front() - y_rev[i]);
  }
  // second part
  for (Int32 i = 0; i < N; i++) {
    BOOST_CHECK(y_out[i + Nreflex] == y_in[i]);
  }
  // third part
  copy(y_in.end() - Nreflex, y_in.end(), y_rev.begin());
  reverse(y_rev.begin(), y_rev.end());
  for (Int32 i = 0; i < Nreflex; i++) {
    BOOST_CHECK(y_out[i + N + Nreflex] == 2 * y_in.back() - y_rev[i]);
  }
  print_flux(y_out);

  // TEST Nreflex = N
  y_in = {1., 2., 3., 2., 1.};
  beg = 0;
  N = y_in.size();
  Nreflex = y_in.size();
  y_out = sample.OddMirror(y_in.begin(), y_in.end(), Nreflex, y_in.front(),
                           y_in.back());
  BOOST_CHECK(y_out.size() == 2 * Nreflex + y_in.size());
  y_rev.resize(Nreflex);
  // first part
  copy(y_in.begin(), y_in.begin() + Nreflex, y_rev.begin());
  reverse(y_rev.begin(), y_rev.end());
  for (Int32 i = 0; i < Nreflex; i++) {
    BOOST_CHECK(y_out[i] == 2 * y_in.front() - y_rev[i]);
  }
  // second part
  for (Int32 i = 0; i < N; i++) {
    BOOST_CHECK(y_out[i + Nreflex] == y_in[i]);
  }
  // third part
  copy(y_in.end() - Nreflex, y_in.end(), y_rev.begin());
  reverse(y_rev.begin(), y_rev.end());
  for (Int32 i = 0; i < Nreflex; i++) {
    BOOST_CHECK(y_out[i + N + Nreflex] == 2 * y_in.back() - y_rev[i]);
  }
  print_flux(y_out);
}

BOOST_AUTO_TEST_CASE(fitBorder_test) {
  CContinuumIrregularSamplingMedian sample;
  TAxisSampleList sAxis = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
  TAxisSampleList fAxis = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
  CSpectrumSpectralAxis spectralAxis(sAxis);
  CSpectrumFluxAxis fluxAxis(fAxis);
  CSpectrum spectra = CSpectrum(spectralAxis, fluxAxis);

  BOOST_TEST_MESSAGE("TEST FIT BORDER");
  TFloat64List flux = fluxAxis.GetSamplesVector();
  Float64 leftBorder = sample.FitBorder(spectra, 0, 10, false);
  BOOST_TEST_MESSAGE("Left : " << leftBorder << "\n");
  BOOST_CHECK_CLOSE(leftBorder, 0., 1e-6);
  Float64 rightBorder = sample.FitBorder(spectra, 0, 10, true);
  BOOST_TEST_MESSAGE("Right : " << rightBorder << "\n");
  BOOST_CHECK_CLOSE(rightBorder, 10., 1e-6);

  fAxis = {0., 1., 2., 3., 4., 5., 4., 3., 2., 1., 0.};
  CSpectrumFluxAxis fluxAxis2(fAxis);
  CSpectrum spectra2 = CSpectrum(spectralAxis, fluxAxis2);
  TFloat64List flux2 = fluxAxis2.GetSamplesVector();
  leftBorder = sample.FitBorder(spectra2, 0, 4, false);
  BOOST_TEST_MESSAGE("Left : " << leftBorder << "\n");
  BOOST_CHECK_CLOSE(leftBorder, 0., 1e-6);
  rightBorder = sample.FitBorder(spectra2, 6, 10, true);
  BOOST_TEST_MESSAGE("Right : " << leftBorder << "\n");
  BOOST_CHECK_CLOSE(rightBorder, 0., 1e-6);
}

//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(removeContinuum_test) {
  // Output Flux
  Int32 n = 10;
  CSpectrumFluxAxis noContinuumFluxAxis(n, 0.);

  // CContinuumIrregularSamplingMedian
  CContinuumIrregularSamplingMedian sample;
  sample.SetMedianEvenReflection(false);
  bool result;
  Float32 width;

  // Cspectrum
  Float64 meanResolution;

  // -------- FUNCTIONAL TESTS ------------

  // TEST 1 : length of fAxis < 11 --> return false
  // (test k<=10)
  TAxisSampleList sAxis = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
  TAxisSampleList fAxis = {-0.5, -0.5, -0.5, -0.5, -5.,
                           -0.5, -0.5, -0.5, -0.5, -0.5};
  CSpectrumSpectralAxis spectralAxis(sAxis);
  CSpectrumFluxAxis fluxAxis(fAxis);
  CSpectrum spectra = CSpectrum(spectralAxis, fluxAxis);
  result = sample.RemoveContinuum(spectra, noContinuumFluxAxis);
  BOOST_CHECK(result == false);
  Int32 k0, k1;
  result = sample.FindEffectiveSpectrumBorder(fluxAxis, k0, k1);
  BOOST_CHECK(result == false);

  // TEST 2 : mean smooth defined to -75.
  // (test meanSmoothAmplitude<=0)
  sAxis = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
  fAxis = {-0.5, -0.5, -0.5, -0.5, -0.5, -5., -0.5, -0.5, -0.5, -0.5};
  spectralAxis = CSpectrumSpectralAxis(sAxis);
  fluxAxis = CSpectrumFluxAxis(fAxis);
  spectra = CSpectrum(spectralAxis, fluxAxis);
  width = -75;
  sample.SetMeanKernelWidth(width);
  result = sample.RemoveContinuum(spectra, noContinuumFluxAxis);
  BOOST_CHECK(result == false);

  // TEST 3 : cycles defines to 0
  // (test m_MedianSmoothCycles<=0)
  width = 75;
  sample.SetMeanKernelWidth(width);
  sample.SetMedianCycleCount(0);
  result = sample.RemoveContinuum(spectra, noContinuumFluxAxis);
  BOOST_CHECK(result == false);

  // TEST 4 : median smooth defined to -75.
  // (test medianSmoothAmplitude<=0)
  width = -75;
  sample.SetMedianKernelWidth(width);
  sample.SetMedianCycleCount(5);
  result = sample.RemoveContinuum(spectra, noContinuumFluxAxis);
  BOOST_CHECK(result == false);

  // TEST 5 : spectra with more than ten 0. values at the beginning of spectra
  // (test k<=10)
  width = 75;
  sample.SetMedianKernelWidth(width);
  sAxis = {1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10.,
           11., 12., 13., 14., 15., 16., 17., 18., 19., 20.};
  fAxis = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
           1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
  spectralAxis = CSpectrumSpectralAxis(sAxis);
  fluxAxis = CSpectrumFluxAxis(fAxis);
  spectra = CSpectrum(spectralAxis, fluxAxis);
  result = sample.RemoveContinuum(spectra, noContinuumFluxAxis);
  BOOST_CHECK(result == false);
  result = sample.FindEffectiveSpectrumBorder(fluxAxis, k0, k1);
  BOOST_CHECK(result == false);

  // -------- SPECTRA WITH TWO PARTS TESTS ------------

  // TEST 6 : even mirror
  width = 75;
  sample.SetMeanKernelWidth(width);
  sample.SetMedianKernelWidth(width);
  sAxis = {0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10.,
           11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21.,
           22., 23., 24., 25., 26., 27., 28., 29., 30.};
  fAxis = {0.,  0.,  0.3, 0.5, 0.2, 0.5, 0.3, 0.2, 0.4, 0.5, 0.3,
           0.2, 0.4, 1.1, 1.2, 1.8, 2.5, 2.2, 1.1, 1.2, 1.1, 1.2,
           1.1, 1.2, 1.1, 1.2, 1.1, 1.2, 0.,  0.,  0.};
  spectralAxis = CSpectrumSpectralAxis(sAxis);
  fluxAxis = CSpectrumFluxAxis(fAxis);
  spectra = CSpectrum(spectralAxis, fluxAxis);
  result = sample.RemoveContinuum(spectra, noContinuumFluxAxis);
  BOOST_CHECK(result == true);
  BOOST_TEST_MESSAGE("TEST REMOVE CONTINUUM WITH ODD REFLECTION");
  print_flux(noContinuumFluxAxis.GetSamplesVector());
  result = sample.FindEffectiveSpectrumBorder(fluxAxis, k0, k1);
  BOOST_CHECK(k0 == 2 && k1 == 27);

  // TEST 7 : odd mirror
  sample.SetMedianEvenReflection(true);
  result = sample.RemoveContinuum(spectra, noContinuumFluxAxis);
  BOOST_CHECK(result == true);
  BOOST_TEST_MESSAGE("TEST REMOVE CONTINUUM WITH EVEN REFLECTION");
  print_flux(noContinuumFluxAxis.GetSamplesVector());
  result = sample.FindEffectiveSpectrumBorder(fluxAxis, k0, k1);
  BOOST_CHECK(k0 == 2 && k1 == 27);
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(ProcessRemoveContinuum) {
  // Output Flux
  Int32 n = 10;
  CSpectrumFluxAxis noContinuumFluxAxis(n, 0.);
  CSpectrumFluxAxis noContinuumFluxAxis2(n, 0.);

  // CContinuumIrregularSamplingMedian
  CContinuumIrregularSamplingMedian sample;
  CContinuumIrregularSamplingMedian sample2;
  bool result;

  // Spectrum
  TAxisSampleList sAxis = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.};
  TAxisSampleList fAxis = {-0.5, -0.5, -0.5, -0.5, -5., -0.5,
                           -0.5, -0.5, -0.5, -0.5, -0.5};
  CSpectrumSpectralAxis spectralAxis(sAxis);
  CSpectrumFluxAxis fluxAxis(fAxis);
  CSpectrum spectra = CSpectrum(spectralAxis, fluxAxis);

  // Test
  result = sample.RemoveContinuum(spectra, noContinuumFluxAxis);
  BOOST_CHECK(result == true);
  result = sample2.ProcessRemoveContinuum(spectra, noContinuumFluxAxis2,
                                          spectra.GetMeanResolution());
  BOOST_CHECK(result == true);
  BOOST_TEST_MESSAGE("TEST ProcessRemoveContinuum");
  for (int i = 0; i < noContinuumFluxAxis.GetSamplesCount(); i++) {
    BOOST_CHECK_CLOSE(abs(noContinuumFluxAxis.GetSamplesVector()[i] -
                          noContinuumFluxAxis2.GetSamplesVector()[i]),
                      0.0, 1e-6);
  }
  print_flux(noContinuumFluxAxis.GetSamplesVector());
  print_flux(noContinuumFluxAxis2.GetSamplesVector());
}

BOOST_AUTO_TEST_SUITE_END()