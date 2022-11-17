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
  TInt32Range range1;
  TInt32Range range2(0, 3);
  TInt32Range range3;
  bool result;

  // Test on : CRange(const std::vector<T> & v)
  TInt32List myVector = {1, 2, 3};
  TInt32Range range4(myVector);
  BOOST_CHECK(range4 == TInt32Range(1, 3));

  // test on : Bool GetIsEmpty() const
  BOOST_CHECK(range2.GetIsEmpty() == false);
  range1.Set(1, 4);

  // test on : const T& GetBegin() const
  BOOST_CHECK(range1 == TInt32Range(1, 4));
  range1.SetBegin(0);
  range1.SetEnd(5);
  BOOST_CHECK(range1 == TInt32Range(0, 5));
  range1 = range2 - 1;
  BOOST_CHECK(range1 == TInt32Range(-1, 2));

  // test on : T Interpolate( Float64 t )
  BOOST_CHECK(range2.Interpolate(2.0) == 6);

  // test on : T GetLength() const
  BOOST_CHECK(range2.GetLength() == 3);

  // test on : static bool Intersect(...)
  //  range1 = [-5, -1]
  //  range2 = [0, 3]
  range1.Set(-5, -1);
  result = TInt32Range::Intersect(range1, range2, range3);
  BOOST_CHECK(result == false);

  // range1 = [4, 8]
  // range2 = [0, 3]
  range1.Set(4, 8);
  result = TInt32Range::Intersect(range1, range2, range3);
  BOOST_CHECK(result == false);

  // range1 = [-1, 2]
  // range2 = [0, 3]
  range1.Set(-1, 2);
  result = TInt32Range::Intersect(range1, range2, range3);
  BOOST_CHECK((result == true && range3 == TInt32Range(0, 2)));

  // range1 = [0, 3]
  // range2 = [0, 3]
  range1.Set(0, 3);
  result = TInt32Range::Intersect(range1, range2, range3);
  BOOST_CHECK((result == true && range3 == TInt32Range(0, 3)));

  // range1 = [1, 2]
  // range2 = [0, 3]
  range1.Set(1, 2);
  result = TInt32Range::Intersect(range1, range2, range3);
  BOOST_CHECK((result == true && range3 == TInt32Range(1, 2)));

  // range1 = [-1, 4]
  // range2 = [0, 3]
  range1.Set(-1, 4);
  result = TInt32Range::Intersect(range1, range2, range3);
  BOOST_CHECK((result == true && range3 == TInt32Range(0, 3)));

  // range1 = [2, 4]
  // range2 = [0, 3]
  range1.Set(2, 4);
  result = TInt32Range::Intersect(range1, range2, range3);
  BOOST_CHECK((result == true && range3 == TInt32Range(2, 3)));

  // test on : Bool IntersectWith(...)
  //  fRange1 = [1, 2]
  //  fRange2 = [0, 3]
  TFloat64Range fRange1;
  TFloat64Range fRange2;
  fRange1.Set(1., 2.);
  fRange2.Set(0., 3.);
  result = fRange2.IntersectWith(fRange1);
  BOOST_CHECK((result == true && fRange2 == TFloat64Range(1., 2.)));

  // fRange1 = [-1, 2]
  // fRange2 = [0, 3]
  fRange1.Set(-1., 2.);
  fRange2.Set(0., 3.);
  result = fRange2.IntersectWith(fRange1);
  BOOST_CHECK((result == true && fRange2 == TFloat64Range(0., 2.)));

  // fRange1 = [1, 4]
  // fRange2 = [0, 3]
  fRange1.Set(1., 4.);
  fRange2.Set(0., 3.);
  result = fRange2.IntersectWith(fRange1);
  BOOST_CHECK((result == true && fRange2 == TFloat64Range(1., 3.)));

  // fRange1 = [-1, 4]
  // fRange2 = [0, 3]
  fRange1.Set(-1., 4.);
  fRange2.Set(0., 3.);
  result = fRange2.IntersectWith(fRange1);
  BOOST_CHECK((result == true && fRange2 == TFloat64Range(0., 3.)));

  // fRange1 = [-4, -1]
  // fRange2 = [0, 3]
  fRange1.Set(-4., -1.);
  fRange2.Set(0., 3.);
  result = fRange2.IntersectWith(fRange1);
  BOOST_CHECK((result == false && fRange2 == TFloat64Range(0., 3.)));

  // fRange1 = [4, 7]
  // fRange2 = [0, 3]
  fRange1.Set(4., 7.);
  fRange2.Set(0., 3.);
  result = fRange2.IntersectWith(fRange1);
  BOOST_CHECK((result == false && fRange2 == TFloat64Range(0., 3.)));
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
  BOOST_CHECK((range1 == TFloat64Range(4., 4.)));

  range1.Set(1, 3);
  range2 = range1 + 3;
  BOOST_CHECK((range2 == TFloat64Range(4., 6.)));

  range1.Set(1, 3);
  range2 = range1 - 1;
  BOOST_CHECK((range2 == TFloat64Range(0., 2.)));

  range1.Set(1, 3);
  range2 = range1 * 3;
  BOOST_CHECK((range2 == TFloat64Range(3., 9.)));

  range1.Set(1, 3);
  range2 = range1 / 2;
  BOOST_CHECK((range2 == TFloat64Range(0.5, 1.5)));
}

//-----------------------------------------------------------------------------

// INT32 TEST CASE 1:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       TRUE            FALSE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOver_Int32_test1) {
  Float64 delta = 3.0;

  // Known vector
  TInt32List myVector = {1};

  TInt32Range myRange(1, 1);
  BOOST_CHECK(myRange.GetIsEmpty() == true);

  TInt32List functionResult = myRange.SpreadOver(delta);
  BOOST_CHECK(myVector == functionResult);

  functionResult = myRange.SpreadOver_backward(delta);
  BOOST_CHECK(myVector == functionResult);
}

//-----------------------------------------------------------------------------

// INT32 TEST CASE 2:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       FALSE            TRUE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOver_test2) {
  Float64 delta = 0.0;

  // Known vector
  TInt32List myVector = {1};
  TInt32List myVector2 = {10};

  TInt32Range myRange(1, 10);
  TInt32List functionResult = myRange.SpreadOver(delta);
  BOOST_CHECK(myVector == functionResult);
  functionResult = myRange.SpreadOver_backward(delta);
  BOOST_CHECK(myVector2 == functionResult);
}

//-----------------------------------------------------------------------------

// INT32 TEST CASE 3:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       FALSE            FALSE              TRUE
BOOST_AUTO_TEST_CASE(SpreadOver_test3) {
  Float64 delta = 7.0;

  // Known vector
  TInt32List myVector = {1};
  TInt32List myVector2 = {3};

  TInt32Range myRange(1, 3);
  BOOST_CHECK(myRange.GetIsEmpty() == false);

  TInt32List functionResult = myRange.SpreadOver(delta);
  BOOST_CHECK(myVector == functionResult);
  functionResult = myRange.SpreadOver_backward(delta);
  BOOST_CHECK(myVector2 == functionResult);
}

//-----------------------------------------------------------------------------

// INT32 TEST CASE 4:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       FALSE            FALSE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOver_test4) {
  Float64 delta = 2.0;

  TInt32Range myRange(1, 11);
  BOOST_CHECK(myRange.GetIsEmpty() == false);
  BOOST_CHECK(myRange.GetLength() == 10);

  TInt32List functionResult = myRange.SpreadOver(delta);
  BOOST_CHECK(functionResult.front() == myRange.GetBegin());
  for (int i = 1; i < functionResult.size(); i++) {
    BOOST_CHECK(functionResult[i] - functionResult[i - 1] == delta);
  }
  BOOST_CHECK(myRange.GetEnd() - functionResult.back() < delta);

  functionResult = myRange.SpreadOver_backward(delta);
  BOOST_CHECK(functionResult.back() == myRange.GetEnd());
  for (int i = 0; i < functionResult.size() - 1; i++) {
    BOOST_CHECK(functionResult[i + 1] - functionResult[i] == delta);
  }
  BOOST_CHECK(functionResult.front() - myRange.GetBegin() < delta);
}

//-----------------------------------------------------------------------------

// FLOAT64 TEST CASE 1:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       TRUE            FALSE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOver_float_test1) {
  Float64 delta = 3.0;

  // Known vector
  TFloat64List myVector = {1};

  TFloat64Range myRange(1, 1);
  BOOST_CHECK(myRange.GetIsEmpty() == true);

  TFloat64List functionResult = myRange.SpreadOver(delta);
  BOOST_CHECK(myVector == functionResult);
  functionResult = myRange.SpreadOver_backward(delta);
  BOOST_CHECK(myVector == functionResult);
}

//-----------------------------------------------------------------------------

// FLOAT64 TEST CASE 2:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       FALSE            TRUE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOver_float_test2) {
  Float64 delta = 0.0;

  // Known vector
  TFloat64List myVector = {1};
  TFloat64List myVector2 = {10};

  TFloat64Range myRange(1, 10);
  TFloat64List functionResult = myRange.SpreadOver(delta);
  BOOST_CHECK(myVector == functionResult);
  functionResult = myRange.SpreadOver_backward(delta);
  BOOST_CHECK(myVector2 == functionResult);
}

//-----------------------------------------------------------------------------

// FLOAT64 TEST CASE 3:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       FALSE            FALSE              TRUE
BOOST_AUTO_TEST_CASE(SpreadOver_float_test3) {
  Float64 delta = 7.0;

  // Known vector
  TFloat64List myVector = {1};
  TFloat64List myVector2 = {3};

  TFloat64Range myRange(1, 3);
  BOOST_CHECK(myRange.GetIsEmpty() == false);

  TFloat64List functionResult = myRange.SpreadOver(delta);
  BOOST_CHECK(myVector == functionResult);
  functionResult = myRange.SpreadOver_backward(delta);
  BOOST_CHECK(myVector2 == functionResult);
}

//-----------------------------------------------------------------------------

// FLOAT64 TEST CASE 4:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       FALSE            FALSE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOver_float_test4) {
  Float64 delta = 2.0;

  TFloat64Range myRange(1, 11);
  BOOST_CHECK(myRange.GetIsEmpty() == false);
  BOOST_CHECK(myRange.GetLength() == 10);

  TFloat64List functionResult = myRange.SpreadOver(delta);
  BOOST_CHECK_CLOSE(functionResult.front(), myRange.GetBegin(), precision);
  for (int i = 1; i < functionResult.size(); i++) {
    BOOST_CHECK_CLOSE(functionResult[i] - functionResult[i - 1], delta,
                      precision);
  }
  BOOST_CHECK(myRange.GetEnd() - functionResult.back() < delta);

  functionResult = myRange.SpreadOver_backward(delta);
  BOOST_CHECK_CLOSE(functionResult.back(), myRange.GetEnd(), precision);
  for (int i = 1; i < functionResult.size(); i++) {
    BOOST_CHECK_CLOSE(functionResult[i] - functionResult[i - 1], delta,
                      precision);
  }
  BOOST_CHECK(functionResult.front() - myRange.GetBegin() < delta);
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
  TFloat64List myVector = {1};

  // test
  TFloat64Range myRange(1, 1);
  BOOST_CHECK(myRange.GetIsEmpty() == true);

  TFloat64List functionResult = myRange.SpreadOverLog(delta, offset);
  BOOST_CHECK(myVector == functionResult);

  functionResult = myRange.SpreadOverLog_backward(delta, offset);
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
  TFloat64List myVector = {1};
  TFloat64List myVector2 = {3};

  TFloat64Range myRange(1, 3);
  TFloat64List functionResult = myRange.SpreadOverLog(delta, offset);
  BOOST_CHECK(myVector == functionResult);

  functionResult = myRange.SpreadOverLog_backward(delta, offset);
  BOOST_CHECK(myVector2 == functionResult);
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
  TFloat64List myVector = {1};
  TFloat64List myVector2 = {3};

  TFloat64Range myRange(1, 3);
  TFloat64List functionResult = myRange.SpreadOverLog(delta, offset);
  BOOST_CHECK(myVector == functionResult);

  functionResult = myRange.SpreadOverLog_backward(delta, offset);
  BOOST_CHECK(myVector2 == functionResult);
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
  for (int i = 1; i < functionResult.size(); i++) {
    BOOST_CHECK_CLOSE((functionResult[i] + offset) /
                          (functionResult[i - 1] + offset),
                      exp(delta), precision);
  }
  BOOST_CHECK((myRange.GetEnd() + offset) / (functionResult.back() + offset) <
              exp(delta));

  functionResult = myRange.SpreadOverLog_backward(delta, offset);
  BOOST_CHECK_CLOSE(functionResult.back(), myRange.GetEnd(), precision);
  for (int i = 1; i < functionResult.size(); i++) {
    BOOST_CHECK_CLOSE((functionResult[i] + offset) /
                          (functionResult[i - 1] + offset),
                      exp(delta), precision);
  }
  BOOST_CHECK((functionResult.front() + offset) /
                  (myRange.GetBegin() + offset) <
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
  bool ret = TInt32Range::HasIntersection(range, overlappingRange);
  BOOST_CHECK(ret == true);

  TInt32Range nonOverlappingRange(6, 7);
  ret = TInt32Range::HasIntersection(range, nonOverlappingRange);
  BOOST_CHECK(ret == false);
}

BOOST_AUTO_TEST_CASE(removeDuplicates) {
  // no overlapping
  std::vector<TInt32Range> indexRanges{{4, 5}, {1, 2}};
  std::vector<TInt32Range> resRanges =
      TInt32Range::joinIntersections(indexRanges);
  BOOST_CHECK(resRanges.size() == 2);
  BOOST_CHECK(resRanges[1] == TInt32Range(4, 5));
  BOOST_CHECK(resRanges[0] == TInt32Range(1, 2));

  // partial overlapping
  std::vector<TInt32Range> indexRanges1{{4, 5}, {3, 4}, {1, 2}};
  resRanges = TInt32Range::joinIntersections(indexRanges1);
  BOOST_CHECK(resRanges.size() == 2);
  BOOST_CHECK(resRanges[0] == TInt32Range(1, 2));
  BOOST_CHECK(resRanges[1] == TInt32Range(3, 5));

  // complete overlapping
  std::vector<TInt32Range> indexRanges2{{1, 6}, {3, 4}, {1, 2}};
  resRanges = TInt32Range::joinIntersections(indexRanges2);
  BOOST_CHECK(resRanges.size() == 1);
  BOOST_CHECK(resRanges[0] == TInt32Range(1, 6));

  std::vector<TInt32Range> indexRanges3{{1, 2}, {2, 3}, {3, 4}};
  resRanges = TInt32Range::joinIntersections(indexRanges3);
  BOOST_CHECK(resRanges.size() == 1);
  BOOST_CHECK(resRanges[0] == TInt32Range(1, 4));

  std::vector<TInt32Range> indexRanges4{{1, 4}, {3, 6}, {2, 5}};
  resRanges = TInt32Range::joinIntersections(indexRanges4);
  BOOST_CHECK(resRanges.size() == 1);
  BOOST_CHECK(resRanges[0] == TInt32Range(1, 6));
}

BOOST_AUTO_TEST_CASE(doubleequal) {
  BOOST_CHECK((TInt32Range(3, 5) == TInt32Range(3, 5)) == true);
  BOOST_CHECK((TInt32Range(3, 5) == TInt32Range(3, 6)) == false);
  BOOST_CHECK((TInt32Range(3, 5) == TInt32Range(1, 2)) == false);
  BOOST_CHECK((TInt32Range(3, 5) == TInt32Range()) == false);
}

BOOST_AUTO_TEST_CASE(notequal) {
  BOOST_CHECK((TInt32Range(3, 5) != TInt32Range(3, 5)) == false);
  BOOST_CHECK((TInt32Range(3, 5) != TInt32Range(3, 6)) == true);
  BOOST_CHECK((TInt32Range(3, 5) != TInt32Range(1, 2)) == true);
  BOOST_CHECK((TInt32Range(3, 5) != TInt32Range()) == true);
}

BOOST_AUTO_TEST_CASE(Union) {
  TInt32Range a{3, 5};
  TInt32Range b{1, 2};
  TInt32Range c{2, 4};
  TInt32Range resultUnion;

  bool ret = TInt32Range::getUnion(a, b, resultUnion);
  BOOST_CHECK(ret == false);

  ret = TInt32Range::getUnion(a, c, resultUnion);
  BOOST_CHECK(ret == true);
  BOOST_CHECK(resultUnion == TInt32Range(2, 5));

  ret = TInt32Range::getUnion(b, c, resultUnion);
  BOOST_CHECK(ret == true);
  BOOST_CHECK(resultUnion == TInt32Range(1, 4));
}

BOOST_AUTO_TEST_CASE(signCheck) {
  TFloat64Range a(-3, 5);
  TFloat64Range b(3, 5);
  Float64 offset = 1;
  Float64 delta = 1;
  BOOST_CHECK(a.isSameSign(offset) == false);
  BOOST_CHECK_THROW(a.SpreadOverLog(delta), GlobalException);
  BOOST_CHECK_THROW(a.SpreadOverLog_backward(delta), GlobalException);

  BOOST_CHECK(b.isSameSign(offset) == true);
}

BOOST_AUTO_TEST_CASE(Clamp) {
  TFloat64Range a(1.0, 2.0);

  BOOST_CHECK(a.Clamp(1.5) == 1.5);
  BOOST_CHECK(a.Clamp(1.0) == 1.0);
  BOOST_CHECK(a.Clamp(0.1) == 1.0);
  BOOST_CHECK(a.Clamp(2.0) == 2.0);
  BOOST_CHECK(a.Clamp(3.0) == 2.0);
}

BOOST_AUTO_TEST_CASE(spanCenteredWindow_test) {
  TFloat64Range a(1.0, 2.0);
  TFloat64Range b(exp(1.0), exp(2.0));
  Float64 precision = 1e-12;

  Float64 delta = 0.05;

  // linear
  //-------
  Float64 center = 1.5;
  bool logsampling = false;
  TFloat64List vect = a.spanCenteredWindow(center, logsampling, delta);
  BOOST_CHECK(vect.size() == ((a.GetEnd() - a.GetBegin()) / delta) + 1);
  for (std::size_t i = 1; i < vect.size(); i++)
    BOOST_CHECK_CLOSE(vect[i] - vect[i - 1], delta, precision);

  center = 1.21;
  vect = a.spanCenteredWindow(center, logsampling, delta);
  for (std::size_t i = 1; i < vect.size(); i++)
    BOOST_CHECK_CLOSE(vect[i] - vect[i - 1], delta, precision);

  // log
  //-------
  center = exp(1.5);
  logsampling = true;
  vect = b.spanCenteredWindow(center, logsampling, delta);
  BOOST_CHECK(vect.front() >= b.GetBegin() && vect.back() <= b.GetEnd());
  BOOST_CHECK(std::find(vect.begin(), vect.end(), center) != vect.end());

  center = exp(1.21);
  vect = b.spanCenteredWindow(center, true, delta);
  BOOST_CHECK(vect.front() >= b.GetBegin() && vect.back() <= b.GetEnd());
  BOOST_CHECK(std::find(vect.begin(), vect.end(), center) != vect.end());
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
