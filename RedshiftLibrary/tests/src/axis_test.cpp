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

#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/spectrum/axis.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(axis_test)

BOOST_AUTO_TEST_CASE(constructor_test) {
  // constructor
  CSpectrumAxis n0Axis;
  BOOST_CHECK(n0Axis.GetSamplesCount() == 0);

  Float64 n1Array[] = {0.5};
  CSpectrumAxis n1Axis = CSpectrumAxis(n1Array, 1);
  BOOST_CHECK(n1Axis.GetSamplesCount() == 1);
  BOOST_CHECK(n1Axis[0] == 0.5);

  CSpectrumAxis n1bAxis = CSpectrumAxis(2, 0.5);
  BOOST_CHECK(n1bAxis.GetSamplesCount() == 2);
  BOOST_CHECK(n1bAxis[0] == 0.5);
  BOOST_CHECK(n1bAxis[1] == 0.5);

  TFloat64List n1Vector = {0.5, 1.5};
  CSpectrumAxis n1cAxis = CSpectrumAxis(n1Vector);
  BOOST_CHECK(n1cAxis.GetSamplesCount() == 2);
  BOOST_CHECK(n1cAxis[0] == 0.5);
  BOOST_CHECK(n1cAxis[1] == 1.5);

  CSpectrumAxis n1dAxis = CSpectrumAxis(TFloat64List{0.5, 1.5});
  BOOST_CHECK(n1dAxis.GetSamplesCount() == 2);
  BOOST_CHECK(n1dAxis[0] == 0.5);
  BOOST_CHECK(n1dAxis[1] == 1.5);

  // copy and copy assignement
  CSpectrumAxis n2Axis(n1Axis);
  BOOST_CHECK(n2Axis.GetSamplesCount() == 1);
  BOOST_CHECK(n2Axis[0] == 0.5);

  CSpectrumAxis n3Axis;
  n3Axis = n2Axis;
  BOOST_CHECK(n3Axis.GetSamplesCount() == 1);
  BOOST_CHECK(n3Axis[0] == 0.5);

  // move and move assignement
  CSpectrumAxis n4Axis(std::move(n1Axis));
  BOOST_CHECK(n1Axis.GetSamplesCount() == 0);
  BOOST_CHECK(n4Axis.GetSamplesCount() == 1);
  BOOST_CHECK(n4Axis[0] == 0.5);

  CSpectrumAxis n1Axis_b = CSpectrumAxis(n1Array, 1);
  CSpectrumAxis n5Axis;
  n5Axis = std::move(n1Axis_b);
  BOOST_CHECK(n1Axis_b.GetSamplesCount() == 0);
  BOOST_CHECK(n5Axis.GetSamplesCount() == 1);
  BOOST_CHECK(n5Axis[0] == 0.5);
}

BOOST_AUTO_TEST_CASE(operator_test) {
  CSpectrumAxis n0Axis;
  BOOST_CHECK(n0Axis.GetSamples() == NULL);
  // One element
  Float64 n1Array[] = {0.5};
  CSpectrumAxis n10Axis(n1Array, 1);
  BOOST_CHECK(n10Axis[0] == n1Array[0]);
  const CSpectrumAxis n11Axis = CSpectrumAxis(n1Array, 1);
  BOOST_CHECK(n11Axis[0] == n1Array[0]);
  // Two elements
  Float64 n2Array[] = {0., 1.};
  CSpectrumAxis n2Axis(n2Array, 2);
  BOOST_CHECK(n2Axis[0] == n2Array[0]);
  BOOST_CHECK(n2Axis[1] == n2Array[1]);

  // operator *
  n2Axis *= 2;
  BOOST_CHECK(n2Axis[0] == n2Array[0] * 2);
  BOOST_CHECK(n2Axis[1] == n2Array[1] * 2);

  // operator /
  n2Axis /= 2;
  BOOST_CHECK(n2Axis[0] == n2Array[0]);
  BOOST_CHECK(n2Axis[1] == n2Array[1]);
}

BOOST_AUTO_TEST_CASE(SampleCount_test) {
  // No element
  CSpectrumAxis n0Axis;
  BOOST_CHECK(n0Axis.GetSamplesCount() == 0);
  // One element
  CSpectrumAxis n10Axis(1);
  BOOST_CHECK(n10Axis.GetSamplesCount() == 1);
  Float64 n1Array[] = {0.5};
  CSpectrumAxis n11Axis(n1Array, 1);
  BOOST_CHECK(n11Axis.GetSamplesCount() == 1);
  // Two elements
  CSpectrumAxis n20Axis(2);
  BOOST_CHECK(n20Axis.GetSamplesCount() == 2);
  Float64 n2Array[] = {0., 1.};
  const CSpectrumAxis n21Axis(n2Array, 2);
  BOOST_CHECK(n21Axis.GetSamplesCount() == 2);
}

BOOST_AUTO_TEST_CASE(GetSample_test) {
  // No element
  CSpectrumAxis n0Axis;
  BOOST_CHECK(n0Axis.GetSamples() == NULL);
  // One element
  Float64 n1Array[] = {0.5};
  CSpectrumAxis n1Axis(n1Array, 1);
  BOOST_CHECK(n1Axis.GetSamples()[0] == n1Array[0]);
  // Two elements
  Float64 n2Array[] = {0., 1.};
  const CSpectrumAxis n2Axis(n2Array, 2);
  BOOST_CHECK(n2Axis.GetSamples()[0] == n2Array[0]);
  BOOST_CHECK(n2Axis.GetSamples()[1] == n2Array[1]);
}

BOOST_AUTO_TEST_CASE(GetSampleVector_test) {
  TFloat64List vectorRes;
  // No element
  CSpectrumAxis n0Axis;
  vectorRes = n0Axis.GetSamplesVector();
  BOOST_CHECK(n0Axis.isEmpty() == 1);
  BOOST_CHECK(vectorRes.size() == 0);
  // One element
  Float64 n1Array[] = {0.5};
  CSpectrumAxis n1Axis(n1Array, 1);
  vectorRes = n1Axis.GetSamplesVector();
  BOOST_CHECK(vectorRes[0] == 0.5);
  // Two elements
  Float64 n2Array[] = {0., 1.};
  const CSpectrumAxis n2Axis(n2Array, 2);
  vectorRes = n2Axis.GetSamplesVector();
  BOOST_CHECK(vectorRes[0] == n2Array[0]);
  BOOST_CHECK(vectorRes[1] == n2Array[1]);
}

BOOST_AUTO_TEST_CASE(SpectrumAxis_test) {
  // test extract
  Float64 n4Array[] = {0, 1, 2, 3};
  CSpectrumAxis n4Axis(n4Array, 4);

  CSpectrumAxis n2Axis = n4Axis.extract(1, 2);
  BOOST_CHECK(n2Axis.GetSamplesCount() == 2);
  BOOST_CHECK(n2Axis[0] == n4Array[1]);
  BOOST_CHECK(n2Axis[1] == n4Array[2]);

  CSpectrumAxis n0Axis;
  n2Axis = n0Axis.extract(1, 2);
  BOOST_CHECK(n2Axis.GetSamplesCount() == 0);

  // test SetSize
  n4Axis.SetSize(2);
  BOOST_CHECK(n4Axis.GetSamplesCount() == 2);
  BOOST_CHECK(n4Axis[0] == n4Array[0]);
  BOOST_CHECK(n4Axis[1] == n4Array[1]);

  n4Axis.SetSize(4);
  BOOST_CHECK(n4Axis.GetSamplesCount() == 4);
  BOOST_CHECK(n4Axis[0] == n4Array[0]);
  BOOST_CHECK(n4Axis[1] == n4Array[1]);
  BOOST_CHECK(n4Axis[2] == 0);
  BOOST_CHECK(n4Axis[3] == 0);

  // test clear
  n4Axis.clear();
  BOOST_CHECK(n4Axis.GetSamplesCount() == 0);

  // test maskVector
  Float64 n6Array[] = {0, 1, 2, 3, 4, 5};
  TFloat64List n6Mask = {0, 0, 1, 1, 0, 0};
  CSpectrumAxis n6Axis(n6Array, 6);
  TFloat64List outputVector =
      n6Axis.maskVector(n6Mask, n6Axis.GetSamplesVector());

  BOOST_CHECK(outputVector.size() == 2);
  BOOST_CHECK(outputVector[0] == n6Array[2]);
  BOOST_CHECK(outputVector[1] == n6Array[3]);

  n6Mask = {1, 1, 0, 0, 0, 1};
  outputVector = n6Axis.maskVector(n6Mask, n6Axis.GetSamplesVector());
  BOOST_CHECK(outputVector.size() == 3);
  BOOST_CHECK(outputVector[0] == n6Array[0]);
  BOOST_CHECK(outputVector[1] == n6Array[1]);
  BOOST_CHECK(outputVector[2] == n6Array[5]);

  n6Mask.resize(5);
  BOOST_CHECK_THROW(n6Axis.maskVector(n6Mask, n6Axis.GetSamplesVector()),
                    GlobalException);

  // test MaskAxis
  n6Mask.resize(6);
  CSpectrumAxis n6AxisMasked(n6Axis.MaskAxis(n6Mask));
  BOOST_CHECK(n6AxisMasked.GetSamplesCount() == 2);
  BOOST_CHECK(n6AxisMasked[0] == n6Array[0]);
  BOOST_CHECK(n6AxisMasked[1] == n6Array[1]);
}

BOOST_AUTO_TEST_SUITE_END()