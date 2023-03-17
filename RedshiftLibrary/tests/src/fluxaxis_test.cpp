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
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/common/mean.h"
#include "RedshiftLibrary/common/median.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/spectrum/fluxaxis.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <numeric>

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(FluxAxis_test)

Float64 precision = 1e-12;

BOOST_AUTO_TEST_CASE(constructor_test) {

  Int32 n = 10;
  TFloat64List sample_ref(n, 0);
  TFloat64List noiseSample_ref(n, 1);

  CSpectrumFluxAxis object_FluxAxis;
  BOOST_CHECK(object_FluxAxis.GetSamplesCount() == 0);

  CSpectrumFluxAxis object_FluxAxis2(n);
  BOOST_CHECK(object_FluxAxis2.GetSamplesCount() == n);
  BOOST_CHECK(object_FluxAxis2.GetSamplesVector() == sample_ref);
  CSpectrumNoiseAxis spectrumNoiseAxis2 = object_FluxAxis2.GetError();
  BOOST_CHECK(spectrumNoiseAxis2.GetSamplesCount() == n);
  BOOST_CHECK(spectrumNoiseAxis2.GetSamplesVector() == noiseSample_ref);

  CSpectrumAxis spectrumAxis(n);
  CSpectrumNoiseAxis spectrumNoiseAxis(n);
  CSpectrumFluxAxis object_FluxAxis2_b(spectrumAxis, spectrumNoiseAxis);
  BOOST_CHECK(object_FluxAxis2_b.GetSamplesCount() == n);
  BOOST_CHECK(object_FluxAxis2_b.GetSamplesVector() == sample_ref);
  CSpectrumNoiseAxis spectrumNoiseAxis2_b = object_FluxAxis2_b.GetError();
  BOOST_CHECK(spectrumNoiseAxis2_b.GetSamplesCount() == n);
  BOOST_CHECK(spectrumNoiseAxis2_b.GetSamplesVector() == noiseSample_ref);

  sample_ref = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  Float64 Array1[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  CSpectrumFluxAxis object_FluxAxis3(Array1, 10);
  BOOST_CHECK(object_FluxAxis3.GetSamplesCount() == 10);
  BOOST_CHECK(object_FluxAxis3.GetSamplesVector() == sample_ref);
  CSpectrumNoiseAxis spectrumNoiseAxis3 = object_FluxAxis3.GetError();
  BOOST_CHECK(spectrumNoiseAxis3.GetSamplesVector() == noiseSample_ref);

  TFloat64List sampleIn1 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  CSpectrumFluxAxis object_FluxAxis3_b(sampleIn1);
  BOOST_CHECK(object_FluxAxis3_b.GetSamplesCount() == 10);
  BOOST_CHECK(object_FluxAxis3_b.GetSamplesVector() == sample_ref);
  CSpectrumNoiseAxis spectrumNoiseAxis3_b = object_FluxAxis3_b.GetError();
  BOOST_CHECK(spectrumNoiseAxis3_b.GetSamplesVector() == noiseSample_ref);

  CSpectrumFluxAxis object_FluxAxis3_c({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
  BOOST_CHECK(object_FluxAxis3_c.GetSamplesCount() == 10);
  BOOST_CHECK(object_FluxAxis3_c.GetSamplesVector() == sample_ref);
  CSpectrumNoiseAxis spectrumNoiseAxis3_c = object_FluxAxis3_c.GetError();
  BOOST_CHECK(spectrumNoiseAxis3_c.GetSamplesVector() == noiseSample_ref);

  noiseSample_ref = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
  Float64 Array2[10] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
  CSpectrumFluxAxis object_FluxAxis4(Array1, 10, Array2, 10);
  Int32 s = object_FluxAxis4.GetSamplesCount();
  BOOST_CHECK(object_FluxAxis4.GetSamplesCount() == 10);
  BOOST_CHECK(object_FluxAxis4.GetSamplesVector() == sample_ref);
  CSpectrumNoiseAxis spectrumNoiseAxis4 = object_FluxAxis4.GetError();
  BOOST_CHECK(spectrumNoiseAxis4.GetSamplesVector() == noiseSample_ref);

  Float64 Array2b[10] = {2, 4, 6, 8, 10, 12, 14, 16, 18};
  BOOST_CHECK_THROW(CSpectrumFluxAxis object_FluxAxis4_b(Array1, 10, Array2, 9),
                    GlobalException);
}

BOOST_AUTO_TEST_CASE(basic_function_test) {
  TFloat64List sample_ref = {1, 2, 3, 4, 5};
  TFloat64List noiseSample_ref = {2, 4, 6, 8, 10};

  Float64 Array1[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  TFloat64List Array2 = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
  CSpectrumFluxAxis object_FluxAxis(Array1, 10);

  //-------------//
  // test setError
  CSpectrumNoiseAxis spectrumNoiseAxis(Array2);
  object_FluxAxis.setError(spectrumNoiseAxis);
  BOOST_CHECK(object_FluxAxis.GetError().GetSamplesVector() == Array2);

  //-------------//
  // test extract
  CSpectrumFluxAxis object_FluxAxis2 = object_FluxAxis.extract(0, 4);
  BOOST_CHECK(object_FluxAxis2.GetSamplesCount() == 5);
  BOOST_CHECK(object_FluxAxis2.GetSamplesVector() == sample_ref);
  CSpectrumNoiseAxis spectrumNoiseAxis2 = object_FluxAxis2.GetError();
  BOOST_CHECK(spectrumNoiseAxis2.GetSamplesVector() == noiseSample_ref);

  //-------------//
  // test SetSize
  object_FluxAxis2.SetSize(10);
  BOOST_CHECK(object_FluxAxis2.GetSamplesCount() == 10);
  spectrumNoiseAxis2 = object_FluxAxis2.GetError();
  BOOST_CHECK(spectrumNoiseAxis2.GetSamplesCount() == 10);

  //-------------//
  // test clear
  object_FluxAxis2.clear();
  BOOST_CHECK(object_FluxAxis2.GetSamplesCount() == 0);
  spectrumNoiseAxis2 = object_FluxAxis2.GetError();
  BOOST_CHECK(spectrumNoiseAxis2.GetSamplesCount() == 0);
}

BOOST_AUTO_TEST_CASE(ApplyMedianSmooth_test) {
  TFloat64List sample_ref = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  CSpectrumFluxAxis object_CSpectrumFluxAxis(sample_ref);

  // kernelHalfWidth = 0 -> false
  bool resultApplyMedianSmoothcas1 =
      object_CSpectrumFluxAxis.ApplyMedianSmooth(0);
  BOOST_CHECK(resultApplyMedianSmoothcas1 == false);
  BOOST_TEST_MESSAGE(
      "resultApplyMedianSmooth cas1:" << resultApplyMedianSmoothcas1);

  // kernelHalfWidth + 1 > sample size -> false
  BOOST_TEST_MESSAGE("object_CSpectrumFluxAxis.GetSamplesCount"
                     << object_CSpectrumFluxAxis.GetSamplesCount());

  bool resultApplyMedianSmoothcas2 =
      object_CSpectrumFluxAxis.ApplyMedianSmooth(15);
  BOOST_CHECK(resultApplyMedianSmoothcas2 == false);
  BOOST_TEST_MESSAGE(
      "resultApplyMedianSmooth cas2:" << resultApplyMedianSmoothcas2);

  // kernelHalfWidth = 1
  // i = 0    -> median between 2 first value = 1.5
  // i = 1..8 -> median between 3 value (i-1, i, i+1) = i
  // i = 9    -> median between 2 last value = 9.5
  bool resultApplyMedianSmoothcas3 =
      object_CSpectrumFluxAxis.ApplyMedianSmooth(1);
  BOOST_CHECK(resultApplyMedianSmoothcas3 == true);
  TFloat64List sample_out = object_CSpectrumFluxAxis.GetSamplesVector();
  // i = 0 : n_points = 2 in CMedian::Find
  BOOST_CHECK(sample_out[0] == (sample_ref[0] + sample_ref[1]) / 2.);
  // i in [1,8] : n_points = 3 in CMedian::Find
  BOOST_CHECK_EQUAL_COLLECTIONS(sample_ref.begin() + 1, sample_ref.end() - 1,
                                sample_out.begin() + 1, sample_out.end() - 1);
  // i = 9 : n_points = 2 in CMedian::Find
  BOOST_CHECK(sample_out[9] == (sample_ref[8] + sample_ref[9]) / 2.);
  BOOST_TEST_MESSAGE(
      "resultApplyMedianSmooth cas3:" << resultApplyMedianSmoothcas3);

  // kernelHalfWidth = 9
  // smooth median computed using all sample, for an even sample :
  //    -> median = (sample[N/2] + sample[N/2 - 1]) / 2
  // i=0..9 -> sample_out[i] = 5.5
  sample_ref = TFloat64List(10, 5.5);
  bool resultApplyMedianSmoothcas4 =
      object_CSpectrumFluxAxis.ApplyMedianSmooth(9);
  BOOST_CHECK(resultApplyMedianSmoothcas4 == true);
  sample_out = object_CSpectrumFluxAxis.GetSamplesVector();
  BOOST_CHECK(sample_out == sample_ref);
  BOOST_TEST_MESSAGE(
      "resultApplyMedianSmooth cas4:" << resultApplyMedianSmoothcas4);
}

BOOST_AUTO_TEST_CASE(ApplyMeanSmooth_test) {
  TFloat64List sample_ref = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  CSpectrumFluxAxis object_CSpectrumFluxAxis(sample_ref);

  // kernelHalfWidth = 0 -> false
  bool resultApplyMeanSmoothcas1 = object_CSpectrumFluxAxis.ApplyMeanSmooth(0);
  BOOST_CHECK(resultApplyMeanSmoothcas1 == false);
  BOOST_TEST_MESSAGE(
      "resultApplyMeanSmooth cas1:" << resultApplyMeanSmoothcas1);

  // kernelHalfWidth + 1 > sample size -> false
  BOOST_TEST_MESSAGE("object_CSpectrumFluxAxis.GetSamplesCount"
                     << object_CSpectrumFluxAxis.GetSamplesCount());

  bool resultApplyMeanSmoothcas2 = object_CSpectrumFluxAxis.ApplyMeanSmooth(15);
  BOOST_CHECK(resultApplyMeanSmoothcas2 == false);
  BOOST_TEST_MESSAGE(
      "resultApplyMeanSmooth cas2:" << resultApplyMeanSmoothcas2);

  // kernelHalfWidth = 1
  // i = 0    -> mean between 2 first value = 1.5
  // i = 1..8 -> mean between 3 value (i-1, i, i+1) = i
  // i = 9    -> mean between 2 last value = 9.5
  sample_ref = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  object_CSpectrumFluxAxis = CSpectrumFluxAxis(sample_ref);
  bool resultApplyMeanSmoothcas3 = object_CSpectrumFluxAxis.ApplyMeanSmooth(1);
  BOOST_CHECK(resultApplyMeanSmoothcas3 == true);
  TFloat64List sample_out = object_CSpectrumFluxAxis.GetSamplesVector();
  // i = 0 :
  BOOST_CHECK(sample_out[0] == (sample_ref[0] + sample_ref[1]) / 2.);
  // i in [1,8] :
  BOOST_CHECK_EQUAL_COLLECTIONS(sample_ref.begin() + 1, sample_ref.end() - 1,
                                sample_out.begin() + 1, sample_out.end() - 1);
  // i = 9 :
  BOOST_CHECK(sample_out[9] == (sample_ref[8] + sample_ref[9]) / 2.);
  BOOST_TEST_MESSAGE(
      "resultApplyMeanSmooth cas3:" << resultApplyMeanSmoothcas3);

  // kernelHalfWidth = 9
  // smooth mean computed using all sample :
  //    -> mean = (sample[0] + ... + sample[N-1]) / N
  // i=0..9 -> sample_out[i] = 5.5
  object_CSpectrumFluxAxis = CSpectrumFluxAxis(sample_ref);
  sample_ref = TFloat64List(10, 5.5);
  bool resultApplyMeanSmoothcas4 = object_CSpectrumFluxAxis.ApplyMeanSmooth(9);
  BOOST_CHECK(resultApplyMeanSmoothcas4 == true);
  sample_out = object_CSpectrumFluxAxis.GetSamplesVector();
  BOOST_CHECK(sample_out == sample_ref);
  BOOST_TEST_MESSAGE(
      "resultApplyMeanSmooth cas4:" << resultApplyMeanSmoothcas4);
}

BOOST_AUTO_TEST_CASE(ComputeMeanAndSDev_test) {

  TFloat64List sample_ref = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  TFloat64List error_ref(10, 0.);
  CSpectrumFluxAxis object_CSpectrumFluxAxis(sample_ref.data(), 10,
                                             error_ref.data(), 10);

  Float64 mean;
  Float64 sdev;
  bool result;

  //--------------------//
  // test ComputeMeanAndSDevWithoutError

  // 1st case : size of mask != size of sample
  CMask mask(5);
  BOOST_CHECK_THROW(
      object_CSpectrumFluxAxis.ComputeMeanAndSDevWithoutError(mask, mean, sdev),
      GlobalException);

  // 2nd case : mask = 0
  mask.SetSize(10);
  result =
      object_CSpectrumFluxAxis.ComputeMeanAndSDevWithoutError(mask, mean, sdev);
  BOOST_CHECK(result == false);
  BOOST_CHECK(mean != mean);
  BOOST_CHECK(sdev != sdev);

  // 3rd case : mask = 1 for i=2 & i=3
  mask[2] = 1;
  mask[3] = 1;
  result =
      object_CSpectrumFluxAxis.ComputeMeanAndSDevWithoutError(mask, mean, sdev);
  Float64 sdev_ref = sqrt((30 - 35) * (30 - 35) + (40 - 35) * (40 - 35));
  BOOST_CHECK_CLOSE(mean, 35, precision);
  BOOST_CHECK_CLOSE(sdev, sdev_ref, precision);

  //--------------------//
  // test ComputeMeanAndSDevWithError

  // 1st case : size of mask != size of sample
  mask.SetSize(5);
  BOOST_CHECK_THROW(
      object_CSpectrumFluxAxis.ComputeMeanAndSDevWithError(mask, mean, sdev),
      GlobalException);

  // 2nd case : mask = 0
  mask.SetSize(10);
  for (Int32 i = 0; i < mask.GetMasksCount(); i++) {
    mask[i] = 0;
  }
  result =
      object_CSpectrumFluxAxis.ComputeMeanAndSDevWithError(mask, mean, sdev);
  BOOST_CHECK(result == false);
  BOOST_CHECK(mean != mean);
  BOOST_CHECK(sdev != sdev);

  // 3rd case : mask = 1 for i=2 & i=3
  mask[2] = 1;
  mask[3] = 1;
  error_ref = TFloat64List(10, 0.5); // weight = 4
  object_CSpectrumFluxAxis.setError(error_ref);
  result =
      object_CSpectrumFluxAxis.ComputeMeanAndSDevWithError(mask, mean, sdev);
  sdev_ref = sqrt((4 * (30 - 35) * (30 - 35) + 4 * (40 - 35) * (40 - 35)) /
                  (8 - 32 / 8));
  BOOST_CHECK_CLOSE(mean, 35, precision);
  BOOST_CHECK_CLOSE(sdev, sdev_ref, precision);

  //--------------------//
  // test ComputeMeanAndSDev

  // With Error -> ComputeMeanAndSDevWithError
  result = object_CSpectrumFluxAxis.ComputeMeanAndSDev(mask, mean, sdev);
  BOOST_CHECK_CLOSE(mean, 35, precision);
  BOOST_CHECK_CLOSE(sdev, sdev_ref, precision);

  // Without Error -> ComputeMeanAndSDevWithoutError
  error_ref = TFloat64List(10, 0.0);
  object_CSpectrumFluxAxis.setError(error_ref);
  result = object_CSpectrumFluxAxis.ComputeMeanAndSDev(mask, mean, sdev);
  sdev_ref = sqrt((30 - 35) * (30 - 35) + (40 - 35) * (40 - 35));
  BOOST_CHECK_CLOSE(mean, 35, precision);
  BOOST_CHECK_CLOSE(sdev, sdev_ref, precision);
}

BOOST_AUTO_TEST_CASE(ComputeRMSDiff_test) {
  Int32 n = 10;

  //--------------------//
  // test ComputeRMSDiff

  // size of sampleA != size of sampleB

  TFloat64List sampleA = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
  CSpectrumFluxAxis object_FluxAxisA(sampleA);

  TFloat64List sampleB = {2., 4., 6., 8., 10., 12., 14., 16., 18.};
  CSpectrumFluxAxis object_FluxAxisB(sampleB);

  BOOST_CHECK_THROW(object_FluxAxisA.ComputeRMSDiff(object_FluxAxisB),
                    GlobalException);

  // size of sampleA = size of sampleB
  sampleB = {2., 4., 6., 8., 10., 12., 14., 16., 18., 20};
  object_FluxAxisB = CSpectrumFluxAxis(sampleB);

  Float64 resultComputeRMSDiff =
      object_FluxAxisA.ComputeRMSDiff(object_FluxAxisB);
  BOOST_TEST_MESSAGE(
      "object_FluxAxisA.GetSamples()[1]=" << object_FluxAxisA.GetSamples()[1]);
  BOOST_TEST_MESSAGE(
      "object_FluxAxisB.GetSamples()[1]=" << object_FluxAxisB.GetSamples()[1]);

  Float64 er2 = 0.f;
  Float64 er = 0.f;

  Float64 weight = (Float64)n;
  for (int j = 0; j < 10; j++) {
    er2 +=
        (object_FluxAxisA.GetSamples()[j] - object_FluxAxisB.GetSamples()[j]) *
        (object_FluxAxisA.GetSamples()[j] - object_FluxAxisB.GetSamples()[j]) /
        weight;
  }
  er = sqrt(er2);

  BOOST_CHECK_CLOSE(resultComputeRMSDiff, er, precision);
}

BOOST_AUTO_TEST_CASE(Subtract_test) {

  //--------------------//
  // test Subtract

  // size of sampleA != size of sampleB

  TFloat64List sampleA = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
  CSpectrumFluxAxis object_FluxAxisA(sampleA);

  TFloat64List sampleB = {2., 4., 6., 8., 10., 12., 14., 16., 18.};
  CSpectrumFluxAxis object_FluxAxisB(sampleB);

  BOOST_CHECK_THROW(object_FluxAxisA.Subtract(object_FluxAxisB),
                    GlobalException);

  // size of sampleA = size of sampleB
  sampleB = {2., 4., 6., 8., 10., 12., 14., 16., 18., 20};
  object_FluxAxisB = CSpectrumFluxAxis(sampleB);

  bool resultSubtract = object_FluxAxisA.Subtract(object_FluxAxisB);

  TFloat64List sample_ref = {-1., -2., -3., -4., -5., -6., -7., -8., -9., -10.};

  BOOST_CHECK(sample_ref == object_FluxAxisA.GetSamplesVector());
}

BOOST_AUTO_TEST_CASE(ComputeMaxAbsValue) {
  //--------------------//
  // test computeMaxAbsValue
  TFloat64List sampleA = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
  CSpectrumFluxAxis object_FluxAxisA(sampleA);

  Float64 maxAbsVal = object_FluxAxisA.computeMaxAbsValue(0, 9);
  BOOST_CHECK_CLOSE(maxAbsVal, 10., precision);

  maxAbsVal = object_FluxAxisA.computeMaxAbsValue(3, 5);
  BOOST_CHECK_CLOSE(maxAbsVal, 6., precision);

  // check with negative values
  TFloat64List sampleB = {1., 2., 3., 4., 5., 6., 7., 8., 9., -10.};
  CSpectrumFluxAxis object_FluxAxisB(sampleB);

  maxAbsVal = object_FluxAxisB.computeMaxAbsValue(0, 9);
  BOOST_CHECK_CLOSE(maxAbsVal, 10., precision);
}

BOOST_AUTO_TEST_CASE(Invert_test) {
  //--------------------//
  // test Invert
  TFloat64List sampleA = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
  CSpectrumFluxAxis object_FluxAxisA(sampleA);
  bool resultInvert = object_FluxAxisA.Invert();

  TFloat64List sample_ref = {-1., -2., -3., -4., -5., -6., -7., -8., -9., -10.};

  BOOST_CHECK(sample_ref == object_FluxAxisA.GetSamplesVector());
}

BOOST_AUTO_TEST_SUITE_END()
