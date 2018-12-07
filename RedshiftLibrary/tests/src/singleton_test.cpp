#include <boost/test/unit_test.hpp>
#include <RedshiftLibrary/common/singleton.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>

using namespace NSEpic;
using namespace std;

template <typename T>
class SingletonTest: public CSingleton<SingletonTest<T>>
{
public:
    void run();
};

template <typename T> void SingletonTest<T>::run(){
    CSingleton<T> singleton = CSingleton<T>();
    BOOST_CHECK(singleton.IsCreated() == true);
    BOOST_TEST_MESSAGE("TODO singleton.GetInstance()");
    //BOOST_CHECK(singleton.GetInstance() == singleton);
}


BOOST_AUTO_TEST_SUITE(singleton_test)

//-----------------------------------------------------------------------------

//test on : static Bool IsCreated()
BOOST_AUTO_TEST_CASE(IsCreated_test)
{
    SingletonTest<Int32> testInt;
    SingletonTest<Float64> testFloat;
    SingletonTest<std::string> testString;
    SingletonTest<TFloat64Range> testRange;

    testInt.run();
    testFloat.run();
    testString.run();
    testRange.run();

}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
