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
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iterator>

using namespace NSEpic;

BOOST_AUTO_TEST_CASE(Constructor) {
  // No element
  const CSpectrumSpectralAxis n1Axis;
  BOOST_CHECK(n1Axis.GetSamplesCount() == 0);
  // One element initialized to 0. linear scale
  const CSpectrumSpectralAxis n21Axis(1, false);
  BOOST_CHECK(n21Axis.GetSamplesCount() == 1);
  // One element initialized to 0. log scale
  const CSpectrumSpectralAxis n22Axis(1, true);
  BOOST_CHECK(n22Axis.GetSamplesCount() == 1);
  // One element initialized to one element array
  TFloat64List n3Array{1.};
  const CSpectrumSpectralAxis n31Axis(n3Array, false);
  BOOST_CHECK(n31Axis.GetSamplesCount() == 1);
  const CSpectrumSpectralAxis n32Axis(n3Array, true);
  BOOST_CHECK(n32Axis.GetSamplesCount() == 1);

  const TFloat64List n3ArrayList(1, 1);
  const CSpectrumSpectralAxis n33Axis(n3ArrayList);

  // One element from 2 elements array
  TFloat64List n5Array{1., 2.};
  const CSpectrumSpectralAxis n5Axis(n5Array, false);
  BOOST_CHECK(n5Axis.GetSamplesCount() == 2);

  // ShiftByWaveLength linear forward
  const TFloat64List n7Array{2., 3.};
  const CSpectrumSpectralAxis n7Axis(n7Array, false);
  const CSpectrumSpectralAxis n7ShiftForward(
      n7Axis, 10.1, CSpectrumSpectralAxis::nShiftForward);
  BOOST_CHECK(n7ShiftForward.GetSamplesCount() == 2);
  BOOST_CHECK_CLOSE(n7ShiftForward[0], 20.2, 1.e-12);
  BOOST_CHECK_CLOSE(n7ShiftForward[1], 30.3, 1.e-12);

  // ShiftByWaveLength linear backward
  const CSpectrumSpectralAxis n7ShiftBack(
      n7Axis, 2, CSpectrumSpectralAxis::nShiftBackward);
  BOOST_CHECK(n7ShiftBack.GetSamplesCount() == 2);
  BOOST_CHECK_CLOSE(n7ShiftBack[0], 1., 1.e-12);
  BOOST_CHECK_CLOSE(n7ShiftBack[1], 3. / 2, 1.e-12);

  // ShiftByWaveLength log forward
  const CSpectrumSpectralAxis n7AxisLog(n7Array, true);
  const CSpectrumSpectralAxis n7ShiftLogForward(
      n7AxisLog, exp(10.1), CSpectrumSpectralAxis::nShiftForward);
  BOOST_CHECK(n7ShiftLogForward.GetSamplesCount() == 2);
  BOOST_CHECK_CLOSE(n7ShiftLogForward[0], 12.1, 1.e-12);
  BOOST_CHECK_CLOSE(n7ShiftLogForward[1], 13.1, 1.e-12);

  // ShiftByWaveLength log backward
  const CSpectrumSpectralAxis n7ShiftLogBack(
      n7AxisLog, exp(0.5), CSpectrumSpectralAxis::nShiftBackward);
  BOOST_CHECK(n7ShiftLogBack.GetSamplesCount() == 2);
  BOOST_CHECK_CLOSE(n7ShiftLogBack[0], 1.5, 1.e-12);
  BOOST_CHECK_CLOSE(n7ShiftLogBack[1], 2.5, 1.e-12);

  // ShiftByWaveLength zero shift
  const CSpectrumSpectralAxis n7ZeroShift(n7Axis, 0.,
                                          CSpectrumSpectralAxis::nShiftForward);
  BOOST_CHECK(n7ZeroShift.GetSamplesCount() == 2);
  BOOST_CHECK_CLOSE(n7ZeroShift[0], 2., 1.e-12);
  BOOST_CHECK_CLOSE(n7ZeroShift[1], 3., 1.e-12);
}

BOOST_AUTO_TEST_CASE(ApplyOffset) {
  const TFloat64List array{1., 2., 3.};
  CSpectrumSpectralAxis axis(array, false);
  axis.ApplyOffset(1.);

  const CSpectrumSpectralAxis &const_Axis = axis;
  BOOST_CHECK_CLOSE(const_Axis[0], 2., 1.e-12);
  BOOST_CHECK_CLOSE(const_Axis[1], 3., 1.e-12);
  BOOST_CHECK_CLOSE(const_Axis[2], 4., 1.e-12);
}

BOOST_AUTO_TEST_CASE(Operator) {
  const TFloat64List n7Array{1.};
  const CSpectrumSpectralAxis n72Axis(n7Array, false);
  const CSpectrumSpectralAxis n71Axis = n72Axis;
  BOOST_CHECK_CLOSE(n71Axis[0], n72Axis[0], 1.e-12);
}

BOOST_AUTO_TEST_CASE(Resolution) {
  const TFloat64List arr{1., 3., 4., 10., 15., 16.};
  const CSpectrumSpectralAxis axis(arr, false);
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
  const CSpectrumSpectralAxis axis2(array, false);
  BOOST_CHECK_CLOSE(axis2.GetResolution(), 0.0, 1.e-12);
  BOOST_CHECK_CLOSE(axis2.GetMeanResolution(), 0.0, 1.e-12);
}

BOOST_AUTO_TEST_CASE(GetMask) {
  const TFloat64List array{1., 3., 4., 10.};
  CMask mask;
  const TFloat64Range range(0.5, 5.);
  Mask result[] = {1, 1, 1, 0};

  // linear
  const CSpectrumSpectralAxis axis(array, false);
  axis.GetMask(range, mask);
  BOOST_CHECK(std::equal(result, result + 4, mask.GetMasks()));

  // log
  const CSpectrumSpectralAxis axislog(array, true);
  Mask resultlog[] = {1, 0, 0, 0};
  axislog.GetMask(range, mask);
  BOOST_CHECK(std::equal(resultlog, resultlog + 4, mask.GetMasks()));
}

BOOST_AUTO_TEST_CASE(IntersectMaskAndComputeOverlapRate) {
  CMask mask(4);
  TFloat64Range range(0.5, 5.);
  const TFloat64List array{1., 3., 4., 10.};

  mask[0] = 1;
  mask[1] = 1;
  mask[2] = 0;
  mask[3] = 0;

  // linear
  const CSpectrumSpectralAxis axis(array, false);
  BOOST_CHECK_CLOSE(axis.IntersectMaskAndComputeOverlapRate(range, mask),
                    2. / 3, 1e-18);

  // log
  const CSpectrumSpectralAxis axislog(array, true);
  BOOST_CHECK_CLOSE(axislog.IntersectMaskAndComputeOverlapRate(range, mask), 1.,
                    1e-18);

  const TFloat64Range outrange(-5, 0.);
  BOOST_CHECK_CLOSE(axis.IntersectMaskAndComputeOverlapRate(outrange, mask), 0.,
                    1e-18);
}

BOOST_AUTO_TEST_CASE(GetIndexAtWaveLength_and_GetIndexesAtWaveLengthRange) {
  // GetIndexAtWaveLength tests
  const TFloat64List arr{0.0, 2.0, 3.0, 6.0};
  const CSpectrumSpectralAxis axis(arr, false);
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

BOOST_AUTO_TEST_CASE(Scale) {
  const TFloat64List arr{1., 1.};
  CSpectrumSpectralAxis axis(arr, false);

  BOOST_CHECK(axis.IsInLinearScale() == true);
  BOOST_CHECK(axis.IsInLogScale() == false);
  axis.ConvertToLinearScale();
  BOOST_CHECK(axis.IsInLinearScale() == true);
  BOOST_CHECK(axis.IsInLogScale() == false);
  axis.ConvertToLogScale();
  BOOST_CHECK(axis.IsInLinearScale() == false);
  BOOST_CHECK(axis.IsInLogScale() == true);
  BOOST_CHECK_CLOSE(axis[0], 0., 1e-12);

  CSpectrumSpectralAxis axis2(2, true);

  BOOST_CHECK(axis2.IsInLinearScale() == false);
  BOOST_CHECK(axis2.IsInLogScale() == true);
  axis2.ConvertToLogScale();
  BOOST_CHECK(axis2.IsInLinearScale() == false);
  BOOST_CHECK(axis2.IsInLogScale() == true);
  axis2.ConvertToLinearScale();
  BOOST_CHECK(axis2.IsInLinearScale() == true);
  BOOST_CHECK(axis2.IsInLogScale() == false);
  BOOST_CHECK_CLOSE(axis2[0], 1., 1e-12);
}

BOOST_AUTO_TEST_CASE(LambdaRange) {
  // linear
  const TFloat64List arr{0., 1.};
  const CSpectrumSpectralAxis axis(arr, false);
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

  // log
  const TFloat64List arr3{0., 1., 2.};
  const CSpectrumSpectralAxis axis3(arr3, true);
  BOOST_CHECK(axis3.IsInLogScale() == true);
  const TLambdaRange irange9 = axis3.GetLambdaRange();
  BOOST_CHECK_CLOSE(irange9.GetBegin(), exp(0.), 1.e-12);
  BOOST_CHECK_CLOSE(irange9.GetEnd(), exp(2.), 1.e-12);
}

BOOST_AUTO_TEST_CASE(ClampLambdaRange) {
  // Axis from [0. ... 1.]
  TFloat64List n131Array{0., 1.};
  const CSpectrumSpectralAxis axis(n131Array, false);
  TFloat64Range crange;
  // Range [0.1 ... 0.9] inside axis range
  TFloat64Range irange10(0.1, 0.9);
  BOOST_CHECK(axis.ClampLambdaRange(irange10, crange));
  BOOST_CHECK_CLOSE(crange.GetBegin(), 0.1, 1.e-12);
  BOOST_CHECK_CLOSE(crange.GetEnd(), 0.9, 1.e-12);
  // Range empty [0. ... 0.]
  TFloat64Range irange11(0., 0.);
  BOOST_CHECK(~axis.ClampLambdaRange(irange11, crange));
  // Axis empty
  const CSpectrumSpectralAxis n132Axis(2, false);
  BOOST_CHECK(~n132Axis.ClampLambdaRange(irange10, crange));
  // Range [-1.0 ... 0.9] starts before axis [0. ... 1.]
  TFloat64Range irange12(-1.0, 0.9);
  BOOST_CHECK(axis.ClampLambdaRange(irange12, crange));
  BOOST_CHECK_CLOSE(crange.GetBegin(), 0, 1.e-12);
  BOOST_CHECK_CLOSE(crange.GetEnd(), 0.9, 1.e-12);
  // Range [0.1 ... 1.1] ends after axis [0. ... 1.]
  TFloat64Range irange13(0.1, 1.1);
  BOOST_CHECK(axis.ClampLambdaRange(irange13, crange));
  BOOST_CHECK_CLOSE(crange.GetBegin(), 0.1, 1.e-12);
  BOOST_CHECK_CLOSE(crange.GetEnd(), 1., 1.e-12);
  // Range [-2. ... -1.] outside and before axis range
  TFloat64Range irange14(-2., -1.);
  BOOST_CHECK(axis.ClampLambdaRange(irange14, crange));
  BOOST_CHECK_CLOSE(crange.GetBegin(), 0., 1.e-12);
  BOOST_CHECK_CLOSE(crange.GetEnd(), -1., 1.e-12);
  // Range [-2. ... -1.] outside and after axis range
  TFloat64Range irange15(2., 3.);
  BOOST_CHECK(axis.ClampLambdaRange(irange15, crange));
  BOOST_CHECK_CLOSE(crange.GetBegin(), 2., 1.e-12);
  BOOST_CHECK_CLOSE(crange.GetEnd(), 1., 1.e-12);
}

// BOOST_TEST_MESSAGE("Resolution:"<<xAxis.GetResolution());
// BOOST_TEST_MESSAGE("Mean resolution:"<<xAxis.GetMeanResolution());
// BOOST_TEST_MESSAGE("Interpolate value:"<<yInt[0]);
// BOOST_CHECK( 0 == 1 );
BOOST_AUTO_TEST_CASE(isSorted) {
  const TFloat64List arr{0., 1., 4., 3.};
  const CSpectrumSpectralAxis axis(arr, false);
  BOOST_CHECK(axis.isSorted() == false);

  const TFloat64List arr2{0., 0.};
  const CSpectrumSpectralAxis axis2(arr2, false);
  BOOST_CHECK(axis2.isSorted() == false);

  const TFloat64List arr3;
  const CSpectrumSpectralAxis axis3(arr3, false);
  BOOST_CHECK(axis3.isSorted() == true);
}
