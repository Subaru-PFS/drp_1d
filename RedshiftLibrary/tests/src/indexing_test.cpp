#include <RedshiftLibrary/common/indexing.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/execution_monitor.hpp>  
#include <vector>

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(indexing_test)

BOOST_AUTO_TEST_CASE(indexing_test_int)
{
    TInt32List myVector = {0, 2, 2, 3, 4, 4, 5, 6, 6, 7};
    Int32 target = 2, idx;  
    idx = CIndexing<Int32>::getIndex(myVector,target);
    BOOST_CHECK( myVector[idx] == target);
}

BOOST_AUTO_TEST_CASE(indexing_test_float)
{
    TFloat64List myVector = {0.0, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5};
    Float64 target = 2.0, idx;   
    idx = CIndexing<Float64>::getIndex(myVector,target);
    BOOST_CHECK( myVector[idx] == target);
}

bool correctMessage(const std::runtime_error& ex)
{
    BOOST_CHECK_EQUAL(ex.what(), std::string("Could not find index for 2.000000"));
    return true;
}
BOOST_AUTO_TEST_CASE(indexing_test_float_erro)
{
    TFloat64List myVector = {0.0, 2.2, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5};
    const Float64 target = 2.0;
    BOOST_CHECK_EXCEPTION(CIndexing<Float64>::getIndex(myVector,target), std::runtime_error, correctMessage);
}


BOOST_AUTO_TEST_CASE(LowestIndex)
{
    TFloat64List myVector = {0.0, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5};
    const Float64 target = 2.2;
    Int32 i_min = -1;
   
    CIndexing<Float64>::getClosestLowerIndex(myVector,target,i_min);
    BOOST_CHECK( myVector[i_min] <= target);
}

BOOST_AUTO_TEST_CASE(LowerIndex)
{
    TFloat64List myVector = {0.0, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5};
    const Float64 target = 2.019999999;
    Int32 i_min = -1;
   
    i_min = CIndexing<Float64>::getCloserIndex(myVector,target);

    BOOST_CHECK( i_min == 1);
}
BOOST_AUTO_TEST_CASE(LowerIndex2)
{
    TFloat64List myVector = {0.0, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5};
    const Float64 target = 2.419999999;
    Int32 i_min = -1;
   
    i_min = CIndexing<Float64>::getCloserIndex(myVector,target);

    BOOST_CHECK( i_min == 2);
}
BOOST_AUTO_TEST_CASE(LowerIndex_outsideBorders)
{
    TFloat64List myVector = {0.0, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5};
    Int32 i_min = -1;

    Float64 target = 6.519999999;   
    i_min = CIndexing<Float64>::getCloserIndex(myVector,target);
    BOOST_CHECK( i_min == -1);

    target = -0.019999999;
    i_min = CIndexing<Float64>::getCloserIndex(myVector,target);
    BOOST_CHECK( i_min == -1);

    target = 6.5+ 1E-9;   
    i_min = CIndexing<Float64>::getCloserIndex(myVector,target);
    BOOST_CHECK( i_min == 10);

    target = 0-1E-9;
    i_min = CIndexing<Float64>::getCloserIndex(myVector,target);
    BOOST_CHECK( i_min == 0);
}

/////
}