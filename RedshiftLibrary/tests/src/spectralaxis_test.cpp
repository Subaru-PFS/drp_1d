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
#include <algorithm>
#include <cmath>
#include <iterator>

#include <boost/test/unit_test.hpp>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/line/airvacuum.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(spectralaxis_test)

Float64 precision = 1e-8;

BOOST_AUTO_TEST_CASE(Constructor) {
  // Default
  // -------

  // No element
  const CSpectrumSpectralAxis n1Axis;
  BOOST_CHECK(n1Axis.GetSamplesCount() == 0);

  // copy
  CSpectrumAxis spcAxis(TFloat64List{1, 2});
  CSpectrumSpectralAxis n6Axis(spcAxis);
  BOOST_CHECK(n6Axis.GetSamplesCount() == 2);

  // move
  CSpectrumSpectralAxis n6Axisb(CSpectrumAxis(TFloat64List{1, 2}));
  BOOST_CHECK(n6Axisb.GetSamplesCount() == 2);

  // constructor with size arg
  // -----------------------------

  const CSpectrumSpectralAxis n21Axis(1);
  BOOST_CHECK(n21Axis.GetSamplesCount() == 1);
  BOOST_CHECK(n21Axis[0] == 0.0);

  // constructor with size and value args
  // ------------------------------------

  const CSpectrumSpectralAxis n23Axis(2, 0.5);
  BOOST_CHECK(n23Axis.GetSamplesCount() == 2);
  BOOST_CHECK(n23Axis[0] == 0.5);
  BOOST_CHECK(n23Axis[1] == 0.5);

  // constructor with airvaccum args
  // --------------------------------------------

  // without AirVacuum conversion
  TFloat64List n3Array{12500.};
  const CSpectrumSpectralAxis n31Axis(n3Array);
  BOOST_CHECK(n31Axis.GetSamplesCount() == 1);
  BOOST_CHECK(n31Axis[0] == 12500.);

  // with AirVacuum conversion
  const CSpectrumSpectralAxis n32Axis(n3Array, "morton2000");
  BOOST_CHECK(n32Axis.GetSamplesCount() == 1);
  auto converter1 = CAirVacuumConverter::Get("morton2000");
  TFloat64List lambdaAir = converter1->VacToAir(n32Axis.GetSamplesVector());
  BOOST_CHECK_CLOSE(lambdaAir[0], n3Array[0], precision);

  // with move sample_in
  const CSpectrumSpectralAxis n33Axis(TFloat64List{15000.}, "morton2000");
  BOOST_CHECK(n33Axis.GetSamplesCount() == 1);
  lambdaAir = converter1->VacToAir(n33Axis.GetSamplesVector());
  BOOST_CHECK_CLOSE(lambdaAir[0], 15000., precision);

  // with array_in
  Float64 Array1[1] = {16000.};
  const CSpectrumSpectralAxis n34Axis(Array1, 1, "morton2000");
  BOOST_CHECK(n34Axis.GetSamplesCount() == 1);
  lambdaAir = converter1->VacToAir(n34Axis.GetSamplesVector());
  BOOST_CHECK_CLOSE(lambdaAir[0], Array1[0], precision);

  // constructor with ShiftByWaveLength
  // ----------------------------------

  // ShiftByWaveLength linear forward
  const TFloat64List n7Array{2., 3.};
  const CSpectrumSpectralAxis n7Axis(n7Array);
  const CSpectrumSpectralAxis n7ShiftForward =
      n7Axis.ShiftByWaveLength(10.1, CSpectrumSpectralAxis::nShiftForward);
  BOOST_CHECK(n7ShiftForward.GetSamplesCount() == 2);
  BOOST_CHECK_CLOSE(n7ShiftForward[0], 20.2, 1.e-12);
  BOOST_CHECK_CLOSE(n7ShiftForward[1], 30.3, 1.e-12);
}

BOOST_AUTO_TEST_CASE(ShiftByWaveLength_test) {
  const TFloat64List sample_ref{2., 3.};
  const CSpectrumSpectralAxis spcAxisOrigin(sample_ref);
  CSpectrumSpectralAxis spcAxisShifted;

  // test wavelengthOffset < 0.
  BOOST_CHECK_THROW(spcAxisOrigin.ShiftByWaveLength(
                        -2., CSpectrumSpectralAxis::nShiftForward),
                    AmzException);

  // ShiftByWaveLength linear forward
  spcAxisShifted = spcAxisOrigin.ShiftByWaveLength(
      10.1, CSpectrumSpectralAxis::nShiftForward);
  BOOST_CHECK(spcAxisShifted.GetSamplesCount() == 2);
  BOOST_CHECK_CLOSE(spcAxisShifted[0], 20.2, 1.e-12);
  BOOST_CHECK_CLOSE(spcAxisShifted[1], 30.3, 1.e-12);
  spcAxisShifted.clear();

  // ShiftByWaveLength linear backward
  spcAxisShifted =
      spcAxisOrigin.ShiftByWaveLength(2, CSpectrumSpectralAxis::nShiftBackward);
  BOOST_CHECK(spcAxisShifted.GetSamplesCount() == 2);
  BOOST_CHECK_CLOSE(spcAxisShifted[0], 1., 1.e-12);
  BOOST_CHECK_CLOSE(spcAxisShifted[1], 3. / 2, 1.e-12);
  spcAxisShifted.clear();

  // ShiftByWaveLength zero shift
  spcAxisShifted =
      spcAxisOrigin.ShiftByWaveLength(0., CSpectrumSpectralAxis::nShiftForward);
  BOOST_CHECK(spcAxisShifted.GetSamplesCount() == 2);
  BOOST_CHECK_CLOSE(spcAxisShifted[0], 2., 1.e-12);
  BOOST_CHECK_CLOSE(spcAxisShifted[1], 3., 1.e-12);

  // ShiftByWaveLength
  CSpectrumSpectralAxis spcAxisShifted2(spcAxisOrigin);
  spcAxisShifted2 = spcAxisShifted2.ShiftByWaveLength(
      8., CSpectrumSpectralAxis::nShiftForward);
  BOOST_CHECK(spcAxisShifted2.GetSamplesCount() == 2);
  BOOST_CHECK_CLOSE(spcAxisShifted2[0], 16., 1.e-12);
  BOOST_CHECK_CLOSE(spcAxisShifted2[1], 24., 1.e-12);
}

BOOST_AUTO_TEST_CASE(basic_functions_test) {
  TFloat64List sample_in = {1., 2.};
  CSpectrumSpectralAxis spcAxis(sample_in);

  // operator *=
  //------------

  spcAxis *= 2;
  BOOST_CHECK(spcAxis[0] == sample_in[0] * 2);
  BOOST_CHECK(spcAxis[1] == sample_in[1] * 2);

  // resort data when op < 0
  spcAxis.isSorted();
  BOOST_ASSERT(spcAxis.m_isSorted == true);
  spcAxis *= -2;
  BOOST_ASSERT(spcAxis.m_isSorted == false);
  BOOST_CHECK(spcAxis[0] == sample_in[0] * -4);
  BOOST_CHECK(spcAxis[1] == sample_in[1] * -4);

  //
  sample_in.pop_back();
  CSpectrumSpectralAxis spcAxis_2(sample_in);
  BOOST_ASSERT(indeterminate(spcAxis_2.m_isSorted));
  spcAxis_2 *= 0;
  BOOST_ASSERT(spcAxis_2.m_isSorted == true);
  spcAxis_2.SetSize(3);
  spcAxis_2 *= 0;
  BOOST_ASSERT(spcAxis_2.m_isSorted == false);

  // SetSize
  //--------

  // s < 2
  sample_in = {1., 2., 3.};
  CSpectrumSpectralAxis spcAxis_4(sample_in);
  BOOST_ASSERT(indeterminate(spcAxis_4.m_isSorted));
  spcAxis_4.SetSize(1);
  BOOST_ASSERT(spcAxis_4.m_isSorted == true);

  BOOST_ASSERT(indeterminate(spcAxis_4.m_isLogSampled));
  // extract
  // -------

  sample_in = {1., 2., 3., 4.};
  CSpectrumSpectralAxis spcAxis_5(sample_in);
  spcAxis_5.m_isLogSampled = true;
  spcAxis_5.m_isSorted = false;
  spcAxis_5.m_regularLogSamplingStep = 2.;
  CSpectrumSpectralAxis spcAxis_6 = spcAxis_5.extract(0, 2);
  BOOST_ASSERT(spcAxis_6.m_isLogSampled == true);
  BOOST_ASSERT(spcAxis_6.m_isSorted == false);
  BOOST_CHECK(spcAxis_6.m_regularLogSamplingStep == 2.);
  sample_in.pop_back();
  BOOST_CHECK(spcAxis_6.GetSamplesVector() == sample_in);

  // resetAxisProperties
  spcAxis_6.m_isSorted = true;
  spcAxis_6.m_isLogSampled = true;
  spcAxis_6.resetAxisProperties();
  BOOST_ASSERT(indeterminate(spcAxis_6.m_isLogSampled));
  BOOST_ASSERT(indeterminate(spcAxis_6.m_isSorted));
}

BOOST_AUTO_TEST_CASE(MaskAxis_test) {
  TFloat64List sample_in = {1., 2., 3., 4., 5.};
  TFloat64List mask = {0, 0, 1, 0, 0};

  CSpectrumSpectralAxis spcAxis(sample_in);
  CSpectrumSpectralAxis spcAxisMasked = spcAxis.MaskAxis(mask);
  BOOST_CHECK(spcAxisMasked.GetSamplesCount() == 1);
  BOOST_CHECK(spcAxisMasked[0] == sample_in[2]);
}

BOOST_AUTO_TEST_CASE(ApplyOffset) {
  const TFloat64List array{1., 2., 3.};
  CSpectrumSpectralAxis axis(array);
  axis.ApplyOffset(1.);

  const CSpectrumSpectralAxis &const_Axis = axis;
  BOOST_CHECK_CLOSE(const_Axis[0], 2., 1.e-12);
  BOOST_CHECK_CLOSE(const_Axis[1], 3., 1.e-12);
  BOOST_CHECK_CLOSE(const_Axis[2], 4., 1.e-12);

  axis.ApplyOffset(0.);
  BOOST_CHECK_CLOSE(axis[0], 2., 1.e-12);
  BOOST_CHECK_CLOSE(axis[1], 3., 1.e-12);
  BOOST_CHECK_CLOSE(axis[2], 4., 1.e-12);
}

BOOST_AUTO_TEST_CASE(Operator) {
  const TFloat64List n7Array{1.};
  const CSpectrumSpectralAxis n72Axis(n7Array);
  const CSpectrumSpectralAxis n71Axis = n72Axis;
  BOOST_CHECK_CLOSE(n71Axis[0], n72Axis[0], 1.e-12);
}

BOOST_AUTO_TEST_CASE(Resolution) {
  const TFloat64List arr{1., 3., 4., 10., 15., 16.};
  const CSpectrumSpectralAxis axis(arr);
  BOOST_CHECK_CLOSE(axis.GetResolution(), 2.0, 1.e-12);
  BOOST_CHECK_CLOSE(axis.GetMeanResolution(), 3.0, 1.e-12);
  BOOST_CHECK_CLOSE(axis.GetResolution(-1.0), 2.0, 1.e-12);
  BOOST_CHECK_CLOSE(axis.GetResolution(-0.1), 2.0, 1.e-12);
  BOOST_CHECK_CLOSE(axis.GetResolution(0.0), 2.0, 1.e-12);
  BOOST_CHECK_CLOSE(axis.GetResolution(1.0), 2.0, 1.e-12);
  BOOST_CHECK_CLOSE(axis.GetResolution(1.1), 2.0, 1.e-12);
  BOOST_CHECK_CLOSE(axis.GetResolution(2.0), 2.0, 1.e-12);
  BOOST_CHECK_CLOSE(axis.GetResolution(3.0), 2.0, 1.e-12);
  BOOST_CHECK_CLOSE(axis.GetResolution(3.3), 1.0, 1.e-12);
  BOOST_CHECK_CLOSE(axis.GetResolution(4.0), 1.0, 1.e-12);
  BOOST_CHECK_CLOSE(axis.GetResolution(7.7), 6.0, 1.e-12);
  BOOST_CHECK_CLOSE(axis.GetResolution(10.0), 6.0, 1.e-12);
  BOOST_CHECK_CLOSE(axis.GetResolution(11.5), 5.0, 1.e-12);
  BOOST_CHECK_CLOSE(axis.GetResolution(15.0), 5.0, 1.e-12);
  BOOST_CHECK_CLOSE(axis.GetResolution(15.6), 1.0, 1.e-12);
  BOOST_CHECK_CLOSE(axis.GetResolution(17.2), 1.0, 1.e-12);

  // not enough samples
  const TFloat64List array{10.};
  const CSpectrumSpectralAxis axis2(array);
  BOOST_CHECK_CLOSE(axis2.GetResolution(), 0.0, 1.e-12);
  BOOST_CHECK_CLOSE(axis2.GetMeanResolution(), 0.0, 1.e-12);
}

BOOST_AUTO_TEST_CASE(GetMask) {
  const TFloat64List array{1., 3., 4., 10.};
  CMask mask;
  const TFloat64Range range(0.5, 5.);
  Mask result[] = {1, 1, 1, 0};

  const CSpectrumSpectralAxis axis(array);
  axis.GetMask(range, mask);
  BOOST_CHECK(std::equal(result, result + 4, mask.GetMasks()));
}

BOOST_AUTO_TEST_CASE(IntersectMaskAndComputeOverlapFraction) {
  CMask mask(4);
  TFloat64Range range(0.5, 5.);
  const TFloat64List array{1., 3., 4., 10.};

  mask[0] = 1;
  mask[1] = 1;
  mask[2] = 0;
  mask[3] = 0;

  const CSpectrumSpectralAxis axis(array);
  BOOST_CHECK_CLOSE(axis.IntersectMaskAndComputeOverlapFraction(range, mask),
                    2. / 3, 1e-18);

  const TFloat64Range outrange(-5, 0.);
  BOOST_CHECK_CLOSE(axis.IntersectMaskAndComputeOverlapFraction(outrange, mask),
                    0., 1e-18);
}

BOOST_AUTO_TEST_CASE(GetIndexAtWaveLength_and_GetIndexesAtWaveLengthRange) {
  // GetIndexAtWaveLength tests
  const TFloat64List arr{0.0, 2.0, 3.0, 6.0};
  const CSpectrumSpectralAxis axis(arr);
  BOOST_CHECK(axis.GetIndexAtWaveLength(-1.0) == 0);
  BOOST_CHECK(axis.GetIndexAtWaveLength(0.0) == 0);
  BOOST_CHECK(axis.GetIndexAtWaveLength(1.0) == 1);
  BOOST_CHECK(axis.GetIndexAtWaveLength(1.9) == 1);
  BOOST_CHECK(axis.GetIndexAtWaveLength(2.0) == 1);
  BOOST_CHECK(axis.GetIndexAtWaveLength(2.1) == 2);
  BOOST_CHECK(axis.GetIndexAtWaveLength(3.0) == 2);
  BOOST_CHECK(axis.GetIndexAtWaveLength(5.3) == 3);
  BOOST_CHECK(axis.GetIndexAtWaveLength(6.0) == 3);
  BOOST_CHECK(axis.GetIndexAtWaveLength(7.0) == 3);
  BOOST_TEST_MESSAGE("index:" << axis.GetIndexAtWaveLength(-1.0));
  BOOST_TEST_MESSAGE("index:" << axis.GetIndexAtWaveLength(0.0));
  BOOST_TEST_MESSAGE("index:" << axis.GetIndexAtWaveLength(1.0));
  BOOST_TEST_MESSAGE("index:" << axis.GetIndexAtWaveLength(1.9));
  BOOST_TEST_MESSAGE("index:" << axis.GetIndexAtWaveLength(2.0));
  BOOST_TEST_MESSAGE("index:" << axis.GetIndexAtWaveLength(2.1));
  BOOST_TEST_MESSAGE("index:" << axis.GetIndexAtWaveLength(3.0));
  BOOST_TEST_MESSAGE("index:" << axis.GetIndexAtWaveLength(5.3));
  BOOST_TEST_MESSAGE("index:" << axis.GetIndexAtWaveLength(6.0));
  BOOST_TEST_MESSAGE("index:" << axis.GetIndexAtWaveLength(7.0));

  // GetIndexesAtWaveLengthRange tests
  const TFloat64Range range1(1.9, 3.1);
  TInt32Range irange1 = axis.GetIndexesAtWaveLengthRange(range1);
  BOOST_TEST_MESSAGE("index:" << irange1.GetBegin() << "," << irange1.GetEnd());
  BOOST_CHECK(irange1.GetBegin() == 1);
  BOOST_CHECK(irange1.GetEnd() == 3);
  const TFloat64Range range2(-2.0, -1.0);
  TInt32Range irange2 = axis.GetIndexesAtWaveLengthRange(range2);
  BOOST_TEST_MESSAGE("index:" << irange2.GetBegin() << "," << irange2.GetEnd());
  BOOST_CHECK(irange2.GetBegin() == 0);
  BOOST_CHECK(irange2.GetEnd() == 0);
  const TFloat64Range range3(10.0, 20.0);
  TInt32Range irange3 = axis.GetIndexesAtWaveLengthRange(range3);
  BOOST_TEST_MESSAGE("index:" << irange3.GetBegin() << "," << irange3.GetEnd());
  BOOST_CHECK(irange3.GetBegin() == 3);
  BOOST_CHECK(irange3.GetEnd() == 3);
  const TFloat64Range range4(-1.0, 2.2);
  TInt32Range irange4 = axis.GetIndexesAtWaveLengthRange(range4);
  BOOST_TEST_MESSAGE("index:" << irange4.GetBegin() << "," << irange4.GetEnd());
  BOOST_CHECK(irange4.GetBegin() == 0);
  BOOST_CHECK(irange4.GetEnd() == 2);
  const TFloat64Range range5(2.2, 10.0);
  TInt32Range irange5 = axis.GetIndexesAtWaveLengthRange(range5);
  BOOST_TEST_MESSAGE("index:" << irange5.GetBegin() << "," << irange5.GetEnd());
  BOOST_CHECK(irange5.GetBegin() == 2);
  BOOST_CHECK(irange5.GetEnd() == 3);
}

BOOST_AUTO_TEST_CASE(LambdaRange) {
  const TFloat64List arr{0., 1.};
  const CSpectrumSpectralAxis axis(arr);
  const TLambdaRange irange6 = axis.GetLambdaRange();
  BOOST_CHECK_CLOSE(irange6.GetBegin(), arr[0], 1.e-12);
  BOOST_CHECK_CLOSE(irange6.GetEnd(), arr[1], 1.e-12);
  const CSpectrumSpectralAxis n122Axis(1, false);
  const TLambdaRange irange7 = n122Axis.GetLambdaRange();
  BOOST_CHECK_CLOSE(irange7.GetBegin(), 0., 1.e-12);
  BOOST_CHECK_CLOSE(irange7.GetEnd(), 0., 1.e-12);
  const CSpectrumSpectralAxis axis2;
  const TLambdaRange irange8 = axis2.GetLambdaRange();
  BOOST_CHECK_CLOSE(irange8.GetBegin(), 0., 1.e-12);
  BOOST_CHECK_CLOSE(irange8.GetEnd(), 0., 1.e-12);
}

BOOST_AUTO_TEST_CASE(ClampLambdaRange) {
  // Axis from [0. ... 1.]
  TFloat64List n131Array{0., 1.};
  const CSpectrumSpectralAxis axis(n131Array);
  TFloat64Range crange;
  // Range [0.1 ... 0.9] inside axis range
  TFloat64Range irange10(0.1, 0.9);
  BOOST_CHECK_NO_THROW(axis.ClampLambdaRange(irange10, crange));
  BOOST_CHECK_CLOSE(crange.GetBegin(), 0.1, 1.e-12);
  BOOST_CHECK_CLOSE(crange.GetEnd(), 0.9, 1.e-12);
  // Range empty [0. ... 0.]
  TFloat64Range irange11(0., 0.);
  BOOST_CHECK_THROW(axis.ClampLambdaRange(irange11, crange), AmzException);
  // Axis empty
  const CSpectrumSpectralAxis n132Axis(2, 0.0);
  BOOST_CHECK_THROW(n132Axis.ClampLambdaRange(irange10, crange), AmzException);
  // Range [-1.0 ... 0.9] starts before axis [0. ... 1.]
  TFloat64Range irange12(-1.0, 0.9);
  BOOST_CHECK_NO_THROW(axis.ClampLambdaRange(irange12, crange));
  BOOST_CHECK_CLOSE(crange.GetBegin(), 0, 1.e-12);
  BOOST_CHECK_CLOSE(crange.GetEnd(), 0.9, 1.e-12);
  // Range [0.1 ... 1.1] ends after axis [0. ... 1.]
  TFloat64Range irange13(0.1, 1.1);
  BOOST_CHECK_NO_THROW(axis.ClampLambdaRange(irange13, crange));
  BOOST_CHECK_CLOSE(crange.GetBegin(), 0.1, 1.e-12);
  BOOST_CHECK_CLOSE(crange.GetEnd(), 1., 1.e-12);
  // Range [-2. ... -1.] outside and before axis range
  TFloat64Range irange14(-2., -1.);
  BOOST_CHECK_THROW(axis.ClampLambdaRange(irange14, crange), AmzException);
  // Range [2. ... 3.] outside and after axis range
  TFloat64Range irange15(2., 3.);
  BOOST_CHECK_THROW(axis.ClampLambdaRange(irange15, crange), AmzException);
}

BOOST_AUTO_TEST_CASE(isSorted) {
  const TFloat64List arr{0., 1., 4., 3.};
  const CSpectrumSpectralAxis axis(arr);
  BOOST_CHECK(axis.isSorted() == false);
  // m_isSorted known -> return value without checking
  BOOST_CHECK(axis.isSorted() == false);
  // m_isSorted known -> return value without checking
  BOOST_CHECK(axis.isSorted() == false);

  const TFloat64List arr2{0., 0.};
  const CSpectrumSpectralAxis axis2(arr2);
  BOOST_CHECK(axis2.isSorted() == false);

  const TFloat64List arr3;
  const CSpectrumSpectralAxis axis3(arr3);
  BOOST_CHECK(axis3.isSorted() == true);
}

BOOST_AUTO_TEST_CASE(logSampling_test) {
  // tests without checking (CheckLoglambdaSampling)
  //------------------------------------------------

  TFloat64List arr{0., 1., 2.};
  CSpectrumSpectralAxis spcAxis(arr);
  spcAxis.m_isLogSampled = true;
  spcAxis.m_regularLogSamplingStep = 1;

  BOOST_CHECK(spcAxis.IsLogSampled() == true);
  BOOST_CHECK(spcAxis.IsLogSampled(1) == true);
  BOOST_CHECK(spcAxis.IsLogSampled(1.2) == false);
  BOOST_CHECK(spcAxis.GetlogGridStep() == 1);

  // tests with checking
  //--------------------

  CSpectrumSpectralAxis spcAxisLinear({exp(1.), exp(2.), exp(3.)});
  bool result = spcAxisLinear.CheckLoglambdaSampling();
  BOOST_CHECK(result == true);
  BOOST_ASSERT(spcAxisLinear.m_isLogSampled == true);
  BOOST_CHECK_CLOSE(spcAxisLinear.m_regularLogSamplingStep, 1., 1e-12);

  // tests IsLogSampled
  //-------------------

  // step OK
  spcAxisLinear.resetAxisProperties();
  result = spcAxisLinear.IsLogSampled();
  BOOST_CHECK(result == true);
  BOOST_ASSERT(spcAxisLinear.m_isLogSampled == true);
  BOOST_CHECK_CLOSE(spcAxisLinear.m_regularLogSamplingStep, 1., 1e-12);

  // step KO
  spcAxisLinear[2] = exp(3.5);
  result = spcAxisLinear.IsLogSampled();
  BOOST_CHECK(result == false);
  BOOST_ASSERT(spcAxisLinear.m_isLogSampled == false);

  // tests IsLogSampled with logGridStep
  //------------------------------------

  // step OK
  spcAxisLinear[2] = exp(3.);
  result = spcAxisLinear.IsLogSampled(1.);
  BOOST_CHECK(result == true);
  BOOST_ASSERT(spcAxisLinear.m_isLogSampled == true);
  BOOST_CHECK_CLOSE(spcAxisLinear.m_regularLogSamplingStep, 1., 1e-12);

  // logGridStep KO
  result = spcAxisLinear.IsLogSampled(1.8);
  BOOST_CHECK(result == false);

  // IsLogSampled KO
  spcAxisLinear[2] = exp(3.5);
  result = spcAxisLinear.IsLogSampled(1.);
  BOOST_CHECK(result == false);

  // tests GetlogGridStep
  //---------------------

  // IsLogSampled OK
  spcAxisLinear[2] = exp(3.);
  Float64 regLogStep = spcAxisLinear.GetlogGridStep();
  BOOST_CHECK_CLOSE(regLogStep, 1., 1e-12);

  // IsLogSampled KO
  spcAxisLinear[2] = exp(3.5);
  BOOST_CHECK_THROW(spcAxisLinear.GetlogGridStep(), AmzException);

  // GetLogSamplingIntegerRatio
  //---------------------------
  Float64 modulo;

  // IsLogSampled OK
  spcAxisLinear[2] = exp(3.);
  Int32 ratio = spcAxisLinear.GetLogSamplingIntegerRatio(1., modulo);
  BOOST_CHECK(ratio == 1);
  BOOST_CHECK_CLOSE(modulo, 0., precision);
  ratio = spcAxisLinear.GetLogSamplingIntegerRatio(1.3, modulo);
  BOOST_CHECK(ratio == 1);
  BOOST_CHECK_CLOSE(modulo, 0.3, precision);
  ratio = spcAxisLinear.GetLogSamplingIntegerRatio(1.5, modulo);
  BOOST_CHECK(ratio == 2);
  BOOST_CHECK_CLOSE(modulo, -0.5, precision);
  ratio = spcAxisLinear.GetLogSamplingIntegerRatio(1.8, modulo);
  BOOST_CHECK(ratio == 2);
  BOOST_CHECK_CLOSE(modulo, -0.2, precision);

  // IsLogSampled KO
  spcAxisLinear[2] = exp(3.5);
  BOOST_CHECK_THROW(spcAxisLinear.GetLogSamplingIntegerRatio(1.8, modulo);
                    , AmzException);
}

BOOST_AUTO_TEST_CASE(SubSamplingMask_test) {
  Int32 ssratio;
  TInt32Range range(1., 3.);
  CSpectrumSpectralAxis spcAxis({exp(1.), exp(2.), exp(3.), exp(4.), exp(5.)});

  // not LogSampled
  spcAxis[2] = 3.5;
  ssratio = 1;
  BOOST_CHECK_THROW(spcAxis.GetSubSamplingMask(ssratio, range), AmzException);

  // range bound KO
  spcAxis[2] = exp(3.);
  TInt32Range range2(-5., 2.);
  BOOST_CHECK_THROW(spcAxis.GetSubSamplingMask(ssratio, range2), AmzException);
  TInt32Range range3(0., 8.);
  BOOST_CHECK_THROW(spcAxis.GetSubSamplingMask(ssratio, range3), AmzException);

  // ssratio = 1
  TFloat64List mask = spcAxis.GetSubSamplingMask(ssratio, range);
  BOOST_CHECK(mask.size() == 5);
  TFloat64List mask_ref(5, 1.);
  BOOST_CHECK(mask == mask_ref);

  // ssratio = 2
  ssratio = 2;
  mask = spcAxis.GetSubSamplingMask(ssratio, range);
  BOOST_CHECK(mask.size() == 5);
  mask_ref = {0., 1., 0., 1., 0.};
  BOOST_CHECK(mask == mask_ref);

  TFloat64Range lambdarange(exp(1.9), exp(4.2));
  mask = spcAxis.GetSubSamplingMask(ssratio, lambdarange);
  BOOST_CHECK(mask.size() == 5);
  BOOST_CHECK(mask == mask_ref);

  mask = spcAxis.GetSubSamplingMask(ssratio);
  BOOST_CHECK(mask.size() == 5);
  mask_ref = {1., 0., 1., 0., 1.};
  BOOST_CHECK(mask == mask_ref);
}

BOOST_AUTO_TEST_CASE(RecomputePreciseLoglambda_test) {
  // axis not logSampled
  TFloat64List wrong_samples = {exp(1), exp(2.2), exp(3)};
  CSpectrumSpectralAxis wrong_spcAxis(wrong_samples);
  BOOST_CHECK_THROW(wrong_spcAxis.RecomputePreciseLoglambda(), AmzException);

  // not enough points
  wrong_spcAxis[1] = exp(2);
  BOOST_CHECK_THROW(wrong_spcAxis.RecomputePreciseLoglambda(), AmzException);

  //
  TFloat64Range range(10, 10000);
  TFloat64List samples = range.SpreadOverLogEpsilon(0.01);

  TFloat64List samples_truncated(samples.size(), 0.);
  for (Int32 i = 0; i < samples.size(); i++) {
    samples_truncated[i] = (int)(100 * samples[i]) / 100.0;
  }

  Float64 maxAbsRelativeError = 0.0;
  Float64 logGridStep = log(samples_truncated[1] / samples_truncated[0]);
  for (Int32 i = 1; i < samples_truncated.size(); i++) {
    Float64 _logGridStep = log(samples_truncated[i] / samples_truncated[i - 1]);
    Float64 relativeErrAbs =
        std::abs((_logGridStep - logGridStep) / logGridStep);
    maxAbsRelativeError = std::max(relativeErrAbs, maxAbsRelativeError);
  }

  CSpectrumSpectralAxis spcAxis(samples_truncated);
  spcAxis.RecomputePreciseLoglambda();
  TFloat64List samples_out = spcAxis.GetSamplesVector();
  Float64 maxAbsRelativeError_out = 0.0;
  logGridStep = log(samples_out[1] / samples_out[0]);
  for (Int32 i = 1; i < samples_out.size(); i++) {
    Float64 _logGridStep = log(samples_out[i] / samples_out[i - 1]);
    Float64 relativeErrAbs =
        std::abs((_logGridStep - logGridStep) / logGridStep);
    maxAbsRelativeError_out = std::max(relativeErrAbs, maxAbsRelativeError_out);
  }

  BOOST_CHECK(maxAbsRelativeError_out < maxAbsRelativeError);
}

BOOST_AUTO_TEST_SUITE_END()
