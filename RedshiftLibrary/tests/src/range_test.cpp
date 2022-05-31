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
#include "RedshiftLibrary/common/range.h"
#include <boost/test/unit_test.hpp>
#include <vector>

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(range_test)

//-----------------------------------------------------------------------------
Float64 precision = 1e-12;

BOOST_AUTO_TEST_CASE(range_test_int1) {
  CRange<Int32> range1;
  CRange<Int32> range2(0, 3);
  CRange<Int32> range3;
  bool result;

  // Test on : CRange(const std::vector<T> & v)
  TInt32List myVector = {1, 2, 3};
  CRange<Int32> range4(myVector);
  BOOST_CHECK((range4.GetBegin() == 1 && range4.GetEnd() == 3));

  // test on : Bool GetIsEmpty() const
  BOOST_CHECK(range2.GetIsEmpty() == false);
  range1.Set(1, 4);

  // test on : const T& GetBegin() const
  BOOST_CHECK((range1.GetBegin() == 1 && range1.GetEnd() == 4));
  range1.SetBegin(0);
  range1.SetEnd(5);
  BOOST_CHECK((range1.GetBegin() == 0 && range1.GetEnd() == 5));
  range1 = range2 - 1;
  BOOST_CHECK((range1.GetBegin() == -1 && range1.GetEnd() == 2));

  // test on : T Interpolate( Float64 t )
  BOOST_CHECK(range2.Interpolate(2.0) == 6);

  // test on : T GetLength() const
  BOOST_CHECK(range2.GetLength() == 3);

  // test on : static bool Intersect(...)
  //  range1 = [-5, -1]
  //  range2 = [0, 3]
  range1.Set(-5, -1);
  result = CRange<Int32>::Intersect(range1, range2, range3);
  BOOST_CHECK(result == false);

  // range1 = [4, 8]
  // range2 = [0, 3]
  range1.Set(4, 8);
  result = CRange<Int32>::Intersect(range1, range2, range3);
  BOOST_CHECK(result == false);

  // range1 = [-1, 2]
  // range2 = [0, 3]
  range1.Set(-1, 2);
  result = CRange<Int32>::Intersect(range1, range2, range3);
  BOOST_CHECK(
      (result == true && range3.GetBegin() == 0 && range3.GetEnd() == 2));

  // range1 = [0, 3]
  // range2 = [0, 3]
  range1.Set(0, 3);
  result = CRange<Int32>::Intersect(range1, range2, range3);
  BOOST_CHECK(
      (result == true && range3.GetBegin() == 0 && range3.GetEnd() == 3));

  // range1 = [1, 2]
  // range2 = [0, 3]
  range1.Set(1, 2);
  result = CRange<Int32>::Intersect(range1, range2, range3);
  BOOST_CHECK(
      (result == true && range3.GetBegin() == 1 && range3.GetEnd() == 2));

  // range1 = [-1, 4]
  // range2 = [0, 3]
  range1.Set(-1, 4);
  result = CRange<Int32>::Intersect(range1, range2, range3);
  BOOST_CHECK(
      (result == true && range3.GetBegin() == 0 && range3.GetEnd() == 3));

  // range1 = [2, 4]
  // range2 = [0, 3]
  range1.Set(2, 4);
  result = CRange<Int32>::Intersect(range1, range2, range3);
  BOOST_CHECK(
      (result == true && range3.GetBegin() == 2 && range3.GetEnd() == 3));

  // test on : Bool IntersectWith(...)
  //  fRange1 = [1, 2]
  //  fRange2 = [0, 3]
  CRange<Float64> fRange1;
  CRange<Float64> fRange2;
  fRange1.Set(1., 2.);
  fRange2.Set(0., 3.);
  result = fRange2.IntersectWith(fRange1);
  BOOST_CHECK(
      (result == true && fRange2.GetBegin() == 1. && fRange2.GetEnd() == 2.));

  // fRange1 = [-1, 2]
  // fRange2 = [0, 3]
  fRange1.Set(-1., 2.);
  fRange2.Set(0., 3.);
  result = fRange2.IntersectWith(fRange1);
  BOOST_CHECK(
      (result == true && fRange2.GetBegin() == 0. && fRange2.GetEnd() == 2.));

  // fRange1 = [1, 4]
  // fRange2 = [0, 3]
  fRange1.Set(1., 4.);
  fRange2.Set(0., 3.);
  result = fRange2.IntersectWith(fRange1);
  BOOST_CHECK(
      (result == true && fRange2.GetBegin() == 1. && fRange2.GetEnd() == 3.));

  // fRange1 = [-1, 4]
  // fRange2 = [0, 3]
  fRange1.Set(-1., 4.);
  fRange2.Set(0., 3.);
  result = fRange2.IntersectWith(fRange1);
  BOOST_CHECK(
      (result == true && fRange2.GetBegin() == 0. && fRange2.GetEnd() == 3.));

  // fRange1 = [-4, -1]
  // fRange2 = [0, 3]
  fRange1.Set(-4., -1.);
  fRange2.Set(0., 3.);
  result = fRange2.IntersectWith(fRange1);
  BOOST_CHECK(
      (result == false && fRange2.GetBegin() == 0. && fRange2.GetEnd() == 3.));

  // fRange1 = [4, 7]
  // fRange2 = [0, 3]
  fRange1.Set(4., 7.);
  fRange2.Set(0., 3.);
  result = fRange2.IntersectWith(fRange1);
  BOOST_CHECK(
      (result == false && fRange2.GetBegin() == 0. && fRange2.GetEnd() == 3.));
}

//-----------------------------------------------------------------------------

// TEST ON OPERATORS
BOOST_AUTO_TEST_CASE(Operator_test) {
  TFloat64Range range1(1, 3);
  TFloat64Range range2;
  std::stringstream in;
  std::stringstream out;

  out << range1;
  BOOST_CHECK(out.str() == "[1,3]");

  in.str("4");
  in >> range1;
  BOOST_CHECK((range1.GetBegin() == 4 && range1.GetEnd() == 4));

  range1.Set(1, 3);
  range2 = range1 + 3;
  BOOST_CHECK((range2.GetBegin() == 4 && range2.GetEnd() == 6));

  range1.Set(1, 3);
  range2 = range1 - 1;
  BOOST_CHECK((range2.GetBegin() == 0 && range2.GetEnd() == 2));

  range1.Set(1, 3);
  range2 = range1 * 3;
  BOOST_CHECK((range2.GetBegin() == 3 && range2.GetEnd() == 9));

  range1.Set(1, 3);
  range2 = range1 / 2;
  BOOST_CHECK((range2.GetBegin() == 0.5 && range2.GetEnd() == 1.5));
}

//-----------------------------------------------------------------------------

// INT32 TEST CASE 1:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       TRUE            FALSE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOver_Int32_test1) {
  Float64 delta = 3.0;

  // Known vector
  TInt32List myVector;
  myVector.resize(1);
  myVector[0] = 1;

  CRange<Int32> myRange(1, 1);
  BOOST_CHECK(myRange.GetIsEmpty() == true);

  TInt32List functionResult = myRange.SpreadOver(delta);
  BOOST_CHECK(myVector == functionResult);
}

//-----------------------------------------------------------------------------

// INT32 TEST CASE 2:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       FALSE            TRUE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOver_test2) {
  Float64 delta = 0.0;

  // Known vector
  TInt32List myVector;
  myVector.resize(1);
  myVector[0] = 1;

  CRange<Int32> myRange(1, 10);
  TInt32List functionResult = myRange.SpreadOver(delta);
  BOOST_CHECK(myVector == functionResult);
}

//-----------------------------------------------------------------------------

// INT32 TEST CASE 3:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       FALSE            FALSE              TRUE
BOOST_AUTO_TEST_CASE(SpreadOver_test3) {
  Float64 delta = 7.0;

  // Known vector
  TInt32List myVector;
  myVector.resize(1);
  myVector[0] = 1;

  CRange<Int32> myRange(1, 3);
  BOOST_CHECK(myRange.GetIsEmpty() == false);

  TInt32List functionResult = myRange.SpreadOver(delta);
  BOOST_CHECK(myVector == functionResult);
}

//-----------------------------------------------------------------------------

// INT32 TEST CASE 4:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       FALSE            FALSE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOver_test4) {
  Float64 delta = 2.0;

  CRange<Int32> myRange(1, 11);
  BOOST_CHECK(myRange.GetIsEmpty() == false);
  BOOST_CHECK(myRange.GetLength() == 10);

  TInt32List functionResult = myRange.SpreadOver(delta);
  BOOST_CHECK(functionResult.front() == myRange.GetBegin());
  for (int i = 0; i < functionResult.size() - 1; i++) {
    BOOST_CHECK(functionResult[i + 1] - functionResult[i] == delta);
  }
  BOOST_CHECK(myRange.GetEnd() - functionResult.back() < delta);
}

//-----------------------------------------------------------------------------

// FLOAT64 TEST CASE 1:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       TRUE            FALSE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOver_float_test1) {
  Float64 delta = 3.0;

  // Known vector
  TFloat64List myVector;
  myVector.resize(1);
  myVector[0] = 1;

  CRange<Float64> myRange(1, 1);
  BOOST_CHECK(myRange.GetIsEmpty() == true);

  TFloat64List functionResult = myRange.SpreadOver(delta);
  BOOST_CHECK(myVector == functionResult);
}

//-----------------------------------------------------------------------------

// FLOAT64 TEST CASE 2:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       FALSE            TRUE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOver_float_test2) {
  Float64 delta = 0.0;

  // Known vector
  TFloat64List myVector;
  myVector.resize(1);
  myVector[0] = 1;

  CRange<Float64> myRange(1, 10);
  TFloat64List functionResult = myRange.SpreadOver(delta);
  BOOST_CHECK(myVector == functionResult);
}

//-----------------------------------------------------------------------------

// FLOAT64 TEST CASE 3:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       FALSE            FALSE              TRUE
BOOST_AUTO_TEST_CASE(SpreadOver_float_test3) {
  Float64 delta = 7.0;

  // Known vector
  TFloat64List myVector;
  myVector.resize(1);
  myVector[0] = 1;

  CRange<Float64> myRange(1, 3);
  BOOST_CHECK(myRange.GetIsEmpty() == false);

  TFloat64List functionResult = myRange.SpreadOver(delta);
  BOOST_CHECK(myVector == functionResult);
}

//-----------------------------------------------------------------------------

// FLOAT64 TEST CASE 4:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       FALSE            FALSE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOver_float_test4) {
  Float64 delta = 2.0;

  CRange<Float64> myRange(1, 11);
  BOOST_CHECK(myRange.GetIsEmpty() == false);
  BOOST_CHECK(myRange.GetLength() == 10);

  TFloat64List functionResult = myRange.SpreadOver(delta);
  BOOST_CHECK_CLOSE(functionResult.front(), myRange.GetBegin(), precision);
  for (int i = 0; i < functionResult.size() - 1; i++) {
    BOOST_CHECK_CLOSE(functionResult[i + 1] - functionResult[i], delta,
                      precision);
  }
  BOOST_CHECK(myRange.GetEnd() - functionResult.back() < delta);
}
//-----------------------------------------------------------------------------

// TEST 1
// if(GetIsEmpty() || delta == 0.0  || GetLength() < (GetBegin() +
// offset)*exp(delta) - (GetBegin() + offset) + epsilon)
//       TRUE            FALSE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOverLog_Float_test1) {
  Float64 delta = 0.5;
  Float64 offset = 0.1;

  // reference
  TFloat64List myVector;
  myVector.resize(1);
  myVector[0] = 1;

  // test
  TFloat64Range myRange(1, 1);
  BOOST_CHECK(myRange.GetIsEmpty() == true);

  TFloat64List functionResult = myRange.SpreadOverLog(delta, offset);
  BOOST_CHECK(myVector == functionResult);
}
//-----------------------------------------------------------------------------

// TEST 2
// if(GetIsEmpty() || delta == 0.0  || GetLength() < (GetBegin() +
// offset)*exp(delta) - (GetBegin() + offset) + epsilon)
//       FALSE            TRUE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOverLog_Float_test2) {
  Float64 delta = 0.0;
  Float64 offset = 0.1;

  // Known vector
  TFloat64List myVector;
  myVector.resize(1);
  myVector[0] = 1;

  TFloat64Range myRange(1, 3);
  TFloat64List functionResult = myRange.SpreadOverLog(delta, offset);
  BOOST_CHECK(myVector == functionResult);
}
//-----------------------------------------------------------------------------

// TEST 3
// if(GetIsEmpty() || delta == 0.0  || GetLength() < (GetBegin() +
// offset)*exp(delta) - (GetBegin() + offset) + epsilon)
//       FALSE            FALSE              TRUE
BOOST_AUTO_TEST_CASE(SpreadOverLog_Float_test3) {
  Float64 delta = 5.0;
  Float64 offset = 0.1;

  // Known vector
  TFloat64List myVector;
  myVector.resize(1);
  myVector[0] = 1;

  TFloat64Range myRange(1, 3);
  TFloat64List functionResult = myRange.SpreadOverLog(delta, offset);
  BOOST_CHECK(myVector == functionResult);
}
//-----------------------------------------------------------------------------

// TEST 4
// if(GetIsEmpty() || delta == 0.0  || GetLength() < (GetBegin() +
// offset)*exp(delta) - (GetBegin() + offset) + epsilon)
//       FALSE            FALSE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOverLog_Float_test4) {
  Float64 delta = 0.5;
  Float64 offset = 0.1;

  // test
  TFloat64Range myRange(1, 11);
  BOOST_CHECK(myRange.GetIsEmpty() == false);
  BOOST_CHECK(myRange.GetLength() == 10);

  TFloat64List functionResult = myRange.SpreadOverLog(delta, offset);
  BOOST_CHECK_CLOSE(functionResult.front(), myRange.GetBegin(), precision);
  for (int i = 0; i < functionResult.size() - 1; i++) {
    BOOST_CHECK_CLOSE((functionResult[i + 1] + offset) /
                          (functionResult[i] + offset),
                      exp(delta), precision);
  }
  BOOST_CHECK((myRange.GetEnd() + offset) / (functionResult.back() + offset) <
              exp(delta));
}
//-----------------------------------------------------------------------------

// TEST SpreadOverLogZplusOne --> offset is locked to one
// if(GetIsEmpty() || delta == 0.0  || GetLength() < (GetBegin() +
// offset)*exp(delta) - (GetBegin() + offset) + epsilon)
//       FALSE            FALSE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOverLogZPlusOne_test) {
  Float64 delta = 0.5;

  // test
  TFloat64Range myRange(1, 11);
  BOOST_CHECK(myRange.GetIsEmpty() == false);
  BOOST_CHECK(myRange.GetLength() == 10);

  TFloat64List functionResult = myRange.SpreadOverLogZplusOne(delta);
  TFloat64List functionResult2 = myRange.SpreadOverLog(delta, 1.);
  BOOST_CHECK(functionResult2 == functionResult);
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Enclosing_interval) {
  Float64 delta = 1.0;

  // Known vector
  TFloat64List myVector(15);
  //    myVector.resize(15);

  for (Int32 i = 0; i < myVector.size(); i++) {
    myVector[i] = 1 + delta * i;
  }

  const Float64 target = 8.;
  Int32 i_min = -1;
  Int32 i_max = -1;

  TFloat64Range range = TFloat64Range(6.5, 10.3);

  range.getEnclosingIntervalIndices(myVector, target, i_min, i_max);
  BOOST_CHECK(myVector[i_min] <= range.GetBegin());
  BOOST_CHECK(myVector[i_max] >= range.GetEnd());
  BOOST_CHECK(myVector[i_min + 1] > range.GetBegin());
  BOOST_CHECK(myVector[i_max - 1] < range.GetEnd());

  range = TFloat64Range(6, 10.3);
  range.getEnclosingIntervalIndices(myVector, target, i_min, i_max);

  BOOST_CHECK(myVector[i_min] <= range.GetBegin());
  BOOST_CHECK(myVector[i_max] >= range.GetEnd());
  BOOST_CHECK(myVector[i_min + 1] > range.GetBegin());
  BOOST_CHECK(myVector[i_max - 1] < range.GetEnd());

  range = TFloat64Range(6, 10);
  range.getEnclosingIntervalIndices(myVector, target, i_min, i_max);

  BOOST_CHECK(myVector[i_min] <= range.GetBegin());
  BOOST_CHECK(myVector[i_max] >= range.GetEnd());
  BOOST_CHECK(myVector[i_min + 1] > range.GetBegin());
  BOOST_CHECK(myVector[i_max - 1] < range.GetEnd());

  // Check warning
  range = TFloat64Range(9, 10.3);
  bool result =
      range.getEnclosingIntervalIndices(myVector, target, i_min, i_max);
  BOOST_CHECK(result == 0);

  range = TFloat64Range(6, 7.5);
  result = range.getEnclosingIntervalIndices(myVector, target, i_min, i_max);
  BOOST_CHECK(result == 0);

  range = TFloat64Range(0, 10.3);
  result = range.getEnclosingIntervalIndices(myVector, target, i_min, i_max);
  BOOST_CHECK(result == 0);

  range = TFloat64Range(6, 16);
  result = range.getEnclosingIntervalIndices(myVector, target, i_min, i_max);
  BOOST_CHECK(result == 0);

  // -- TEST WITHOUT TARGET --

  // Check warning
  range = TFloat64Range(0, 10.3);
  result = range.getEnclosingIntervalIndices(myVector, i_min, i_max);
  BOOST_CHECK(result == 0);

  range = TFloat64Range(6, 16);
  result = range.getEnclosingIntervalIndices(myVector, i_min, i_max);
  BOOST_CHECK(result == 0);

  // TEST OK
  range = TFloat64Range(6.5, 10.3);
  range.getEnclosingIntervalIndices(myVector, i_min, i_max);
  BOOST_CHECK(myVector[i_min] <= range.GetBegin());
  BOOST_CHECK(myVector[i_max] >= range.GetEnd());
  BOOST_CHECK(myVector[i_min + 1] > range.GetBegin());
  BOOST_CHECK(myVector[i_max - 1] < range.GetEnd());

  range = TFloat64Range(6, 10.3);
  range.getEnclosingIntervalIndices(myVector, i_min, i_max);

  BOOST_CHECK(myVector[i_min] <= range.GetBegin());
  BOOST_CHECK(myVector[i_max] >= range.GetEnd());
  BOOST_CHECK(myVector[i_min + 1] > range.GetBegin());
  BOOST_CHECK(myVector[i_max - 1] < range.GetEnd());

  range = TFloat64Range(6.5, 10);
  range.getEnclosingIntervalIndices(myVector, i_min, i_max);

  BOOST_CHECK(myVector[i_min] <= range.GetBegin());
  BOOST_CHECK(myVector[i_max] >= range.GetEnd());
  BOOST_CHECK(myVector[i_min + 1] > range.GetBegin());
  BOOST_CHECK(myVector[i_max - 1] < range.GetEnd());

  range = TFloat64Range(6, 10);
  range.getEnclosingIntervalIndices(myVector, i_min, i_max);

  BOOST_CHECK(myVector[i_min] <= range.GetBegin());
  BOOST_CHECK(myVector[i_max] >= range.GetEnd());
  BOOST_CHECK(myVector[i_min + 1] > range.GetBegin());
  BOOST_CHECK(myVector[i_max - 1] < range.GetEnd());
}

BOOST_AUTO_TEST_CASE(Closed_interval) {
  Float64 delta = 1.0;
  TFloat64List myVector(15);
  for (Int32 i = 0; i < myVector.size(); i++) {
    myVector[i] = 1 + delta * i;
  }

  const Float64 target = 8.;
  Int32 i_min = -1, i_max = -1;

  // Check warning
  TFloat64Range range = TFloat64Range(20, 25);
  bool result = range.getClosedIntervalIndices(myVector, i_min, i_max);
  BOOST_CHECK(result == 0);

  range = TFloat64Range(-10, -5);
  result = range.getClosedIntervalIndices(myVector, i_min, i_max);
  BOOST_CHECK(result == 0);

  // range borders belong to orderded values
  range = TFloat64Range(6.5, 10.3);
  range.getClosedIntervalIndices(myVector, i_min, i_max);
  BOOST_CHECK(myVector[i_min] >= range.GetBegin());
  BOOST_CHECK(myVector[i_max] <= range.GetEnd());
  BOOST_CHECK(i_min == 6);
  BOOST_CHECK(i_max == 9);

  i_min = -1, i_max = -1;
  // range borders belong to orderded values
  range = TFloat64Range(6., 10.);
  range.getClosedIntervalIndices(myVector, i_min, i_max);
  BOOST_CHECK(myVector[i_min] <= range.GetBegin());
  BOOST_CHECK(myVector[i_max] >= range.GetEnd());
  BOOST_CHECK(i_min == 5);
  BOOST_CHECK(i_max == 9);

  // range borders correspond to min/max orderded values
  i_min = -1, i_max = -1;
  range = TFloat64Range(1., 15.);
  range.getClosedIntervalIndices(myVector, i_min, i_max);
  BOOST_CHECK(myVector[i_min] <= range.GetBegin());
  BOOST_CHECK(myVector[i_max] >= range.GetEnd());
  BOOST_CHECK(i_min == 0);
  BOOST_CHECK(i_max == 14);

  i_min = -1, i_max = -1;
  range = TFloat64Range(-2., 17.);
  range.getClosedIntervalIndices(myVector, i_min, i_max);
  BOOST_CHECK(i_min == 0);
  BOOST_CHECK(i_max == 14);

  i_min = -1, i_max = -1;
  range = TFloat64Range(0, 1);
  range.getClosedIntervalIndices(myVector, i_min, i_max);
  BOOST_CHECK(i_min == 0);
  BOOST_CHECK(i_max == 0);
}
BOOST_AUTO_TEST_CASE(maskedRange) {
  TInt32Range range(3, 5);
  TFloat64List otherVector = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  TFloat64List mask = {1, 1, 0, 0, 0, 0,
                       0, 0, 1, 1}; // mask deactivating the range from 2 to 7

  TFloat64List ssVector = {0, 1, 8, 9}; // subsampled vector following mask

  TFloat64Range otherRange(otherVector[range.GetBegin()],
                           otherVector[range.GetEnd()]);
  Int32 kstart = -1, kend = -1;
  bool ret = otherRange.getClosedIntervalIndices(ssVector, kstart, kend);
  BOOST_CHECK(ret == false);
}
BOOST_AUTO_TEST_CASE(maskedRange_oneCommon) {
  TInt32Range range(1, 5);
  TFloat64List otherVector = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  TFloat64List mask = {1, 1, 0, 0, 0, 0,
                       0, 0, 1, 1}; // mask deactivating the range from 2 to 7

  TFloat64List ssVector = {0, 1, 8, 9}; // subsampled vector following mask

  TFloat64Range otherRange(otherVector[range.GetBegin()],
                           otherVector[range.GetEnd()]);
  Int32 kstart = -1, kend = -1;
  Int32 ret = otherRange.getClosedIntervalIndices(ssVector, kstart, kend);
  BOOST_CHECK(ret == true);
  BOOST_CHECK(kstart == 1);
  BOOST_CHECK(kend == 1);
}
BOOST_AUTO_TEST_CASE(overlappingRanges) {
  TInt32Range range(3, 5);

  TInt32Range overlappingRange(4, 5);
  bool ret = CRange<Int32>::isoverlapping(range, overlappingRange);
  BOOST_CHECK(ret == true);

  TInt32Range nonOverlappingRange(6, 7);
  ret = CRange<Int32>::isoverlapping(range, nonOverlappingRange);
  BOOST_CHECK(ret == false);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
