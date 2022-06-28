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
#include "RedshiftLibrary/line/airvacuum.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iterator>

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

  // constructor with logScale arg
  // -----------------------------

  // One element initialized to 0. linear scale
  const CSpectrumSpectralAxis n21Axis(1);
  BOOST_CHECK(n21Axis.GetSamplesCount() == 1);
  BOOST_CHECK(n21Axis[0] == 0.0);
  BOOST_CHECK(n21Axis.IsInLinearScale() == 1);
  // One element initialized to 0. log scale
  const CSpectrumSpectralAxis n22Axis(1, true);
  BOOST_CHECK(n22Axis.GetSamplesCount() == 1);
  BOOST_CHECK(n22Axis[0] == 0.0);
  BOOST_CHECK(n22Axis.IsInLogScale() == 1);

  // constructor with size and value args
  // ------------------------------------

  const CSpectrumSpectralAxis n23Axis(2, 0.5);
  BOOST_CHECK(n23Axis.GetSamplesCount() == 2);
  BOOST_CHECK(n23Axis[0] == 0.5);
  BOOST_CHECK(n23Axis[1] == 0.5);
  BOOST_CHECK(n23Axis.IsInLinearScale() == 1);

  // constructor with logScale ans airvaccum args
  // --------------------------------------------

  // without AirVacuum conversion
  TFloat64List n3Array{12500.};
  const CSpectrumSpectralAxis n31Axis(n3Array, true);
  BOOST_CHECK(n31Axis.GetSamplesCount() == 1);
  BOOST_CHECK(n31Axis[0] == 12500.);

  // with AirVacuum conversion
  const CSpectrumSpectralAxis n32Axis(n3Array, false, "Morton2000");
  BOOST_CHECK(n32Axis.GetSamplesCount() == 1);
  auto converter1 = CAirVacuumConverter::Get("Morton2000");
  TFloat64List lambdaAir = converter1->VacToAir(n32Axis.GetSamplesVector());
  BOOST_CHECK_CLOSE(lambdaAir[0], n3Array[0], precision);

  // with move sample_in
  const CSpectrumSpectralAxis n33Axis(TFloat64List{15000.}, false,
                                      "Morton2000");
  BOOST_CHECK(n33Axis.GetSamplesCount() == 1);
  lambdaAir = converter1->VacToAir(n33Axis.GetSamplesVector());
  BOOST_CHECK_CLOSE(lambdaAir[0], 15000., precision);

  // with array_in
  Float64 Array1[1] = {16000.};
  const CSpectrumSpectralAxis n34Axis(Array1, 1, "Morton2000");
  BOOST_CHECK(n34Axis.GetSamplesCount() == 1);
  lambdaAir = converter1->VacToAir(n34Axis.GetSamplesVector());
  BOOST_CHECK_CLOSE(lambdaAir[0], Array1[0], precision);

  // constructor with ShiftByWaveLength
  // ----------------------------------

  // ShiftByWaveLength linear forward
  const TFloat64List n7Array{2., 3.};
  const CSpectrumSpectralAxis n7Axis(n7Array, false);
  const CSpectrumSpectralAxis n7ShiftForward(
      n7Axis, 10.1, CSpectrumSpectralAxis::nShiftForward);
  BOOST_CHECK(n7ShiftForward.GetSamplesCount() == 2);
  BOOST_CHECK_CLOSE(n7ShiftForward[0], 20.2, 1.e-12);
  BOOST_CHECK_CLOSE(n7ShiftForward[1], 30.3, 1.e-12);
}

BOOST_AUTO_TEST_CASE(ShiftByWaveLength_test) {
  const TFloat64List sample_ref{2., 3.};
  const CSpectrumSpectralAxis spcAxisOrigin(sample_ref, false);
  const CSpectrumSpectralAxis spcAxisOriginLog(sample_ref, true);
  CSpectrumSpectralAxis spcAxisShifted;

  // test wavelengthOffset < 0.
  BOOST_CHECK_THROW(
      spcAxisShifted.ShiftByWaveLength(spcAxisOrigin, -2.,
                                       CSpectrumSpectralAxis::nShiftForward),
      GlobalException);

  // ShiftByWaveLength linear forward
  spcAxisShifted.ShiftByWaveLength(spcAxisOrigin, 10.1,
                                   CSpectrumSpectralAxis::nShiftForward);
  BOOST_CHECK(spcAxisShifted.GetSamplesCount() == 2);
  BOOST_CHECK_CLOSE(spcAxisShifted[0], 20.2, 1.e-12);
  BOOST_CHECK_CLOSE(spcAxisShifted[1], 30.3, 1.e-12);
  spcAxisShifted.clear();

  // ShiftByWaveLength linear backward
  spcAxisShifted.ShiftByWaveLength(spcAxisOrigin, 2,
                                   CSpectrumSpectralAxis::nShiftBackward);
  BOOST_CHECK(spcAxisShifted.GetSamplesCount() == 2);
  BOOST_CHECK_CLOSE(spcAxisShifted[0], 1., 1.e-12);
  BOOST_CHECK_CLOSE(spcAxisShifted[1], 3. / 2, 1.e-12);
  spcAxisShifted.clear();

  // ShiftByWaveLength log forward
  spcAxisShifted.ShiftByWaveLength(spcAxisOriginLog, exp(10.1),
                                   CSpectrumSpectralAxis::nShiftForward);
  BOOST_CHECK(spcAxisShifted.GetSamplesCount() == 2);
  BOOST_CHECK_CLOSE(spcAxisShifted[0], 12.1, 1.e-12);
  BOOST_CHECK_CLOSE(spcAxisShifted[1], 13.1, 1.e-12);
  spcAxisShifted.clear();

  // ShiftByWaveLength log backward
  spcAxisShifted.ShiftByWaveLength(spcAxisOriginLog, exp(0.5),
                                   CSpectrumSpectralAxis::nShiftBackward);
  BOOST_CHECK(spcAxisShifted.GetSamplesCount() == 2);
  BOOST_CHECK_CLOSE(spcAxisShifted[0], 1.5, 1.e-12);
  BOOST_CHECK_CLOSE(spcAxisShifted[1], 2.5, 1.e-12);
  spcAxisShifted.clear();

  // ShiftByWaveLength zero shift
  spcAxisShifted.ShiftByWaveLength(spcAxisOrigin, 0.,
                                   CSpectrumSpectralAxis::nShiftForward);
  BOOST_CHECK(spcAxisShifted.GetSamplesCount() == 2);
  BOOST_CHECK_CLOSE(spcAxisShifted[0], 2., 1.e-12);
  BOOST_CHECK_CLOSE(spcAxisShifted[1], 3., 1.e-12);

  // ShiftByWaveLength
  CSpectrumSpectralAxis spcAxisShifted2(spcAxisOriginLog);
  spcAxisShifted2.ShiftByWaveLength(exp(10.),
                                    CSpectrumSpectralAxis::nShiftForward);
  BOOST_CHECK(spcAxisShifted2.GetSamplesCount() == 2);
  BOOST_CHECK_CLOSE(spcAxisShifted2[0], 12., 1.e-12);
  BOOST_CHECK_CLOSE(spcAxisShifted2[1], 13., 1.e-12);
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

  //
  CSpectrumSpectralAxis spcAxis_3(sample_in, true);
  BOOST_CHECK(spcAxis_3.m_isLogScale == true);
  spcAxis_3 *= 2;
  BOOST_ASSERT(indeterminate(spcAxis_3.m_isLogSampled));

  // SetSize
  //--------

  // s < 2
  sample_in = {1., 2., 3.};
  CSpectrumSpectralAxis spcAxis_4(sample_in, true);
  BOOST_ASSERT(indeterminate(spcAxis_4.m_isSorted));
  spcAxis_4.SetSize(1);
  BOOST_ASSERT(spcAxis_4.m_isSorted == true);

  // s > sample size
  spcAxis_4.SetSize(3);
  BOOST_ASSERT(indeterminate(spcAxis_4.m_isSorted));
  BOOST_ASSERT(indeterminate(spcAxis_4.m_isLogSampled));

  // extract
  // -------

  sample_in = {1., 2., 3., 4.};
  CSpectrumSpectralAxis spcAxis_5(sample_in, true);
  BOOST_CHECK(spcAxis_5.m_isLogScale == true);
  spcAxis_5.m_isLogSampled = true;
  spcAxis_5.m_isSorted = false;
  spcAxis_5.m_regularLogSamplingStep = 2.;
  CSpectrumSpectralAxis spcAxis_6 = spcAxis_5.extract(0, 2);
  BOOST_CHECK(spcAxis_6.m_isLogScale == true);
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
  CSpectrumSpectralAxis spcAxisMasked;
  spcAxis.MaskAxis(mask, spcAxisMasked);
  BOOST_CHECK(spcAxisMasked.GetSamplesCount() == 1);
  BOOST_CHECK(spcAxisMasked[0] == sample_in[2]);
}

BOOST_AUTO_TEST_CASE(ApplyOffset) {
  const TFloat64List array{1., 2., 3.};
  CSpectrumSpectralAxis axis(array, false);
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

  CSpectrumSpectralAxis axis2(2, false);

  BOOST_CHECK(axis2.IsInLinearScale() == true);
  BOOST_CHECK(axis2.IsInLogScale() == false);
  axis2.SetLogScale();
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
  BOOST_CHECK(axis.ClampLambdaRange(irange11, crange) == false);
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

BOOST_AUTO_TEST_CASE(isSorted) {
  const TFloat64List arr{0., 1., 4., 3.};
  const CSpectrumSpectralAxis axis(arr, false);
  BOOST_CHECK(axis.isSorted() == false);
  // m_isSorted known -> return value without checking
  BOOST_CHECK(axis.isSorted() == false);

  const TFloat64List arr2{0., 0.};
  const CSpectrumSpectralAxis axis2(arr2, false);
  BOOST_CHECK(axis2.isSorted() == false);

  const TFloat64List arr3;
  const CSpectrumSpectralAxis axis3(arr3, false);
  BOOST_CHECK(axis3.isSorted() == true);
}

BOOST_AUTO_TEST_CASE(logSampling_test) {
  // tests without checking (CheckLoglambdaSampling)
  //------------------------------------------------

  TFloat64List arr{0., 1., 2.};
  CSpectrumSpectralAxis spcAxis(arr, true);
  spcAxis.m_isLogSampled = true;
  spcAxis.m_regularLogSamplingStep = 1;

  BOOST_CHECK(spcAxis.IsLogSampled() == true);
  BOOST_CHECK(spcAxis.IsLogSampled(1) == true);
  BOOST_CHECK(spcAxis.IsLogSampled(1.2) == false);
  BOOST_CHECK(spcAxis.GetlogGridStep() == 1);

  // tests with checking
  //--------------------

  // linear scale
  CSpectrumSpectralAxis spcAxisLinear({exp(1.), exp(2.), exp(3.)}, false);
  bool result = spcAxisLinear.CheckLoglambdaSampling();
  BOOST_CHECK(result == true);
  BOOST_ASSERT(spcAxisLinear.m_isLogSampled == true);
  BOOST_CHECK_CLOSE(spcAxisLinear.m_regularLogSamplingStep, 1., 1e-12);

  // log scale
  CSpectrumSpectralAxis spcAxisLog({1., 2., 3.}, true);
  result = spcAxisLog.CheckLoglambdaSampling();
  BOOST_CHECK(result == true);
  BOOST_ASSERT(spcAxisLog.m_isLogSampled == true);
  BOOST_CHECK_CLOSE(spcAxisLog.m_regularLogSamplingStep, log(3.) / 2., 1e-12);

  // error > 1e-1
  spcAxisLog[2] = 3.5;
  result = spcAxisLog.CheckLoglambdaSampling();
  BOOST_CHECK(result == false);

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
  BOOST_CHECK_THROW(spcAxisLinear.GetlogGridStep(), GlobalException);

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
                    , GlobalException);
}

BOOST_AUTO_TEST_CASE(SubSamplingMask_test) {
  Int32 ssratio;
  TInt32Range range(1., 3.);
  CSpectrumSpectralAxis spcAxis({1., 2., 3., 4., 5.}, true);

  // not LogSampled
  spcAxis[2] = 3.5;
  ssratio = 1;
  BOOST_CHECK_THROW(spcAxis.GetSubSamplingMask(ssratio, range),
                    GlobalException);

  // range bound KO
  spcAxis[2] = 3.;
  TInt32Range range2(-5., 2.);
  BOOST_CHECK_THROW(spcAxis.GetSubSamplingMask(ssratio, range2),
                    GlobalException);
  TInt32Range range3(0., 8.);
  BOOST_CHECK_THROW(spcAxis.GetSubSamplingMask(ssratio, range3),
                    GlobalException);

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

  TFloat64Range lambdarange(1.9, 4.2);
  mask = spcAxis.GetSubSamplingMask(ssratio, lambdarange);
  BOOST_CHECK(mask.size() == 5);
  BOOST_CHECK(mask == mask_ref);

  mask = spcAxis.GetSubSamplingMask(ssratio);
  BOOST_CHECK(mask.size() == 5);
  mask_ref = {1., 0., 1., 0., 1.};
  BOOST_CHECK(mask == mask_ref);
}

BOOST_AUTO_TEST_CASE(RecomputePreciseLoglambda_test) {
  TFloat64List samples(100, 0);
  for (Int32 i = 0; i < 100; i++) {
    samples[i] = exp(i + 1.);
  }
  CSpectrumSpectralAxis spcAxis(samples, false);
  spcAxis.RecomputePreciseLoglambda();
  TFloat64List samples_out = spcAxis.GetSamplesVector();
  for (Int32 i = 0; i < 100; i++) {
    BOOST_CHECK_CLOSE(samples[i], samples_out[i], precision);
  }

  TFloat64List samples_2(100, 0);
  for (Int32 i = 0; i < 100; i++) {
    samples_2[i] = i + 1.;
  }

  CSpectrumSpectralAxis spcAxis_2(samples_2, true);
  spcAxis_2.RecomputePreciseLoglambda();
}

BOOST_AUTO_TEST_SUITE_END()
