#include <boost/test/unit_test.hpp>
#include <RedshiftLibrary/common/range.h>
#include <vector>

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(range_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(range_test_int1)
{
  CRange<Int32> range1;
  CRange<Int32> range2(0, 3);
  CRange<Int32> range3;
  bool result;

  //test on : Bool GetIsEmpty() const
  BOOST_CHECK(range2.GetIsEmpty() == false);
  range1.Set(1, 4);

  //test on : const T& GetBegin() const
  BOOST_CHECK((range1.GetBegin() == 1 && range1.GetEnd() == 4));
  range1.SetBegin(0);
  range1.SetEnd(5);
  BOOST_CHECK((range1.GetBegin() == 0 && range1.GetEnd() == 5));
  range1 = range2 - 1;
  BOOST_CHECK((range1.GetBegin() == -1 && range1.GetEnd() == 2));

  //test on : T Interpolate( Float64 t )
  BOOST_CHECK(range2.Interpolate(2.0) == 6);

  //test on : T GetLength() const
  BOOST_CHECK(range2.GetLength() == 3);

  //test on : static Bool Intersect(...)
  // range1 = [-5, -1]
  // range2 = [0, 3]
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
  BOOST_CHECK((result == true && range3.GetBegin() == 0 && range3.GetEnd() == 2));

  // range1 = [0, 3]
  // range2 = [0, 3]
  range1.Set(0, 3);
  result = CRange<Int32>::Intersect(range1, range2, range3);
  BOOST_CHECK((result == true && range3.GetBegin() == 0 && range3.GetEnd() == 3));

  // range1 = [1, 2]
  // range2 = [0, 3]
  range1.Set(1, 2);
  result = CRange<Int32>::Intersect(range1, range2, range3);
  BOOST_CHECK((result == true && range3.GetBegin() == 1 && range3.GetEnd() == 2));

  // range1 = [-1, 4]
  // range2 = [0, 3]
  range1.Set(-1, 4);
  result = CRange<Int32>::Intersect(range1, range2, range3);
  BOOST_CHECK((result == true && range3.GetBegin() == 0 && range3.GetEnd() == 3));

  // range1 = [2, 4]
  // range2 = [0, 3]
  range1.Set(2, 4);
  result = CRange<Int32>::Intersect(range1, range2, range3);
  BOOST_CHECK((result == true && range3.GetBegin() == 2 && range3.GetEnd() == 3));
}

//-----------------------------------------------------------------------------

// INT32 TEST CASE 1:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       TRUE            FALSE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOver_Int32_test1)
{
    Float64 delta = 3.0;

    // Known vector
    std::vector<Int32> myVector;
    myVector.resize(1);
    myVector[0] = 1;

    CRange<Int32> myRange(1, 1);
    BOOST_CHECK(myRange.GetIsEmpty() == true);

    std::vector<Int32> functionResult = myRange.SpreadOver(delta);
    BOOST_CHECK(functionResult.size() == 1);

    for (int i = 0; i < myVector.size(); ++i)
    {
        BOOST_CHECK(myVector[i] == functionResult[i]);
    }
}

//-----------------------------------------------------------------------------

// INT32 TEST CASE 2:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       FALSE            TRUE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOver_test2)
{
    Float64 delta = 0.0;

    // Known vector
    std::vector<Int32> myVector;
    myVector.resize(1);
    myVector[0] = 1;

    CRange<Int32> myRange(1, 10);
    BOOST_CHECK(myRange.GetIsEmpty() == false);

    std::vector<Int32> functionResult = myRange.SpreadOver(delta);
    BOOST_CHECK(functionResult.size() == 1);

    for (int i = 0; i < myVector.size(); ++i)
    {
        BOOST_CHECK(myVector[i] == functionResult[i]);
    }
}

//-----------------------------------------------------------------------------

// INT32 TEST CASE 3:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       FALSE            FALSE              TRUE
BOOST_AUTO_TEST_CASE(SpreadOver_test3)
{
    Float64 delta = 7.0;

    // Known vector
    std::vector<Int32> myVector;
    myVector.resize(1);
    myVector[0] = 1;

    CRange<Int32> myRange(1, 3);
    BOOST_CHECK(myRange.GetIsEmpty() == false);

    std::vector<Int32> functionResult = myRange.SpreadOver(delta);
    BOOST_CHECK(functionResult.size() == 1);

    for (int i = 0; i < myVector.size(); ++i)
    {
        BOOST_CHECK(myVector[i] == functionResult[i]);
    }
}

//-----------------------------------------------------------------------------

// INT32 TEST CASE 4:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       FALSE            FALSE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOver_test4)
{
    Float64 delta = 2.0;

    // Known vector
    std::vector<Int32> myVector;
    myVector.resize(6);
    for(UInt32 i = 0; i < myVector.size(); i++)
    {
        myVector[i] = 1 + delta * i;
    }

    CRange<Int32> myRange(1, 11);
    BOOST_CHECK(myRange.GetIsEmpty() == false);
    BOOST_CHECK(myRange.GetLength() == 10);

    std::vector<Int32> functionResult = myRange.SpreadOver(delta);
    // BOOST_CHECK(functionResult.size() == 1);
    for (int i = 0; i < myVector.size(); ++i)
    {
        BOOST_CHECK(myVector[i] == functionResult[i]);
    }
}

//-----------------------------------------------------------------------------

// FLOAT64 TEST CASE 1:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       TRUE            FALSE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOver_float_test1)
{
    Float64 delta = 3.0;

    // Known vector
    std::vector<Float64> myVector;
    myVector.resize(1);
    myVector[0] = 1;

    CRange<Float64> myRange(1, 1);
    BOOST_CHECK(myRange.GetIsEmpty() == true);

    std::vector<Float64> functionResult = myRange.SpreadOver(delta);
    BOOST_CHECK(functionResult.size() == 1);

    for (int i = 0; i < myVector.size(); ++i)
    {
        BOOST_CHECK(myVector[i] == functionResult[i]);
    }
}

//-----------------------------------------------------------------------------

// FLOAT64 TEST CASE 2:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       FALSE            TRUE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOver_float_test2)
{
    Float64 delta = 0.0;

    // Known vector
    std::vector<Float64> myVector;
    myVector.resize(1);
    myVector[0] = 1;

    CRange<Float64> myRange(1, 10);
    BOOST_CHECK(myRange.GetIsEmpty() == false);

    std::vector<Float64> functionResult = myRange.SpreadOver(delta);
    BOOST_CHECK(functionResult.size() == 1);

    for (int i = 0; i < myVector.size(); ++i)
    {
        BOOST_CHECK(myVector[i] == functionResult[i]);
    }
}

//-----------------------------------------------------------------------------

// FLOAT64 TEST CASE 3:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       FALSE            FALSE              TRUE
BOOST_AUTO_TEST_CASE(SpreadOver_float_test3)
{
    Float64 delta = 7.0;

    // Known vector
    std::vector<Float64> myVector;
    myVector.resize(1);
    myVector[0] = 1;

    CRange<Float64> myRange(1, 3);
    BOOST_CHECK(myRange.GetIsEmpty() == false);

    std::vector<Float64> functionResult = myRange.SpreadOver(delta);
    BOOST_CHECK(functionResult.size() == 1);

    for (int i = 0; i < myVector.size(); ++i)
    {
        BOOST_CHECK(myVector[i] == functionResult[i]);
    }
}

//-----------------------------------------------------------------------------

// FLOAT64 TEST CASE 4:
// if(GetIsEmpty() || delta == 0.0  || GetLength() < delta)
//       FALSE            FALSE              FALSE
BOOST_AUTO_TEST_CASE(SpreadOver_float_test4)
{
    Float64 delta = 2.0;

    // Known vector
    std::vector<Float64> myVector;
    myVector.resize(6);
    for(UInt32 i = 0; i < myVector.size(); i++)
    {
        myVector[i] = 1 + delta * i;
    }

    CRange<Float64> myRange(1, 11);
    BOOST_CHECK(myRange.GetIsEmpty() == false);
    BOOST_CHECK(myRange.GetLength() == 10);

    std::vector<Float64> functionResult = myRange.SpreadOver(delta);
    // BOOST_CHECK(functionResult.size() == 1);
    for (int i = 0; i < myVector.size(); ++i)
    {
        BOOST_CHECK(myVector[i] == functionResult[i]);
    }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
