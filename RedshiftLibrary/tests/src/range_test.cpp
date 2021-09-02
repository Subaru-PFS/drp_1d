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
#include "RedshiftLibrary/common/range.h"
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


BOOST_AUTO_TEST_CASE(Enclosing_interval)
{
    Float64 delta = 1.0;

    // Known vector
    std::vector<Float64> myVector(15);
    //    myVector.resize(15);
    

    for(UInt32 i = 0; i < myVector.size(); i++)
    {
        myVector[i] = 1 + delta * i;
    }

    const Float64 target = 8.;
    Int32 i_min = -1;
    Int32 i_max = -1;

    TFloat64Range range= TFloat64Range(6.5,10.3);
    
    range.getEnclosingIntervalIndices(myVector,target,i_min,i_max);
    BOOST_CHECK( myVector[i_min] <= range.GetBegin());
    BOOST_CHECK( myVector[i_max] >= range.GetEnd());
    BOOST_CHECK( myVector[i_min+1] > range.GetBegin());
    BOOST_CHECK( myVector[i_max-1] < range.GetEnd());

    range= TFloat64Range(6,10.3);
    range.getEnclosingIntervalIndices(myVector,target,i_min,i_max);

    
    BOOST_CHECK( myVector[i_min] <= range.GetBegin());
    BOOST_CHECK( myVector[i_max] >= range.GetEnd());
    BOOST_CHECK( myVector[i_min+1] > range.GetBegin());
    BOOST_CHECK( myVector[i_max-1] < range.GetEnd());

    range= TFloat64Range(6,10);
    range.getEnclosingIntervalIndices(myVector,target,i_min,i_max);

    BOOST_CHECK( myVector[i_min] <= range.GetBegin());
    BOOST_CHECK( myVector[i_max] >= range.GetEnd());
    BOOST_CHECK( myVector[i_min+1] > range.GetBegin());
    BOOST_CHECK( myVector[i_max-1] < range.GetEnd());

}
BOOST_AUTO_TEST_CASE(Closed_interval)
{
    Float64 delta = 1.0;
    std::vector<Float64> myVector(15);
    for(UInt32 i = 0; i < myVector.size(); i++)
    {
        myVector[i] = 1 + delta * i;
    }

    const Float64 target = 8.;
    Int32 i_min = -1, i_max = -1;
    
    //range borders belong to orderded values
    TFloat64Range range = TFloat64Range(6.5,10.3);    
    range.getClosedIntervalIndices(myVector,i_min,i_max);
    BOOST_CHECK( myVector[i_min] >= range.GetBegin());
    BOOST_CHECK( myVector[i_max] <= range.GetEnd());
    BOOST_CHECK( i_min == 6);
    BOOST_CHECK( i_max == 9);
    
    i_min = -1, i_max = -1;
    //range borders belong to orderded values
    range= TFloat64Range(6.,10.);
    range.getClosedIntervalIndices(myVector,i_min,i_max);     
    BOOST_CHECK( myVector[i_min] <= range.GetBegin());
    BOOST_CHECK( myVector[i_max] >= range.GetEnd());
    BOOST_CHECK( i_min == 5);
    BOOST_CHECK( i_max == 9);

    //range borders correspond to min/max orderded values
    i_min = -1, i_max = -1;
    range= TFloat64Range(1.,15.);
    range.getClosedIntervalIndices(myVector,i_min,i_max); 
    BOOST_CHECK( myVector[i_min] <= range.GetBegin());
    BOOST_CHECK( myVector[i_max] >= range.GetEnd());
    BOOST_CHECK( i_min == 0);
    BOOST_CHECK( i_max == 14);

    i_min = -1, i_max = -1;
    range = TFloat64Range(-2.,17.);    
    range.getClosedIntervalIndices(myVector,i_min,i_max);
    BOOST_CHECK( i_min == 0);
    BOOST_CHECK( i_max == 14);

    i_min = -1, i_max = -1;
    range = TFloat64Range(0,1);    
    range.getClosedIntervalIndices(myVector,i_min,i_max);
    BOOST_CHECK( i_min == 0);
    BOOST_CHECK( i_max == 0);
}
BOOST_AUTO_TEST_CASE(maskedRange)
{
    TInt32Range range(3,5);
    TFloat64List otherVector= {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    TFloat64List mask       = {1, 1, 0, 0, 0, 0, 0, 0, 1, 1}; //mask deactivating the range from 2 to 7
    
    TFloat64List ssVector={0, 1, 8, 9};//subsampled vector following mask
    
    TFloat64Range otherRange(otherVector[range.GetBegin()], otherVector[range.GetEnd()]);
    Int32 kstart = -1, kend = -1;
    bool ret = otherRange.getClosedIntervalIndices(ssVector, kstart, kend);
    BOOST_CHECK( ret == false);
}
BOOST_AUTO_TEST_CASE(maskedRange_oneCommon)
{
    TInt32Range range(1,5);
    TFloat64List otherVector= {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    TFloat64List mask       = {1, 1, 0, 0, 0, 0, 0, 0, 1, 1}; //mask deactivating the range from 2 to 7
    
    TFloat64List ssVector={0, 1, 8, 9};//subsampled vector following mask
    
    TFloat64Range otherRange(otherVector[range.GetBegin()], otherVector[range.GetEnd()]);
    Int32 kstart = -1, kend = -1;
    Int32 ret = otherRange.getClosedIntervalIndices(ssVector, kstart, kend);
    BOOST_CHECK( ret == true);
    BOOST_CHECK( kstart == 1);
    BOOST_CHECK( kend == 1);
}
//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
