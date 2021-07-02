#include <boost/test/unit_test.hpp>
#include "RedshiftLibrary/common/singleton.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"

using namespace NSEpic;
using namespace std;

template <typename T>
class SingletonType: public CSingleton<SingletonType<T>> 
{ 
    public:
        T m;

    private:
        friend class CSingleton<SingletonType<T>>;
        SingletonType() = default;
        ~SingletonType() = default;
};

template <typename T>
class SingletonTest
{
public:
    void creation();
    void run(const T & a, const T & b );
};

template <typename T> void SingletonTest<T>::creation(){
    BOOST_CHECK(SingletonType<T>::GetInstance().m == T());    
}

template <typename T> void SingletonTest<T>::run(const T & a, const T & b  ){
    T &s1 = SingletonType<T>::GetInstance().m;
    T &s2 = SingletonType<T>::GetInstance().m;
    s1 = a;
    BOOST_CHECK(SingletonType<T>::GetInstance().m == a );
    BOOST_CHECK(s2 == a );
    s2 = b;
    BOOST_CHECK(SingletonType<T>::GetInstance().m == b );
    BOOST_CHECK(s1 == b );
}

//typedef SingletonTest<TFloat64Range> SRange;
template <>
void SingletonTest<TFloat64Range>::creation(){
    BOOST_CHECK(SingletonType<TFloat64Range>::GetInstance().m.GetBegin() == TFloat64Range().GetBegin());
    BOOST_CHECK(SingletonType<TFloat64Range>::GetInstance().m.GetBegin() == TFloat64Range().GetEnd());    
}

template<>
void SingletonTest<TFloat64Range>::run(const TFloat64Range& a, const TFloat64Range& b  ){
    TFloat64Range &s1 = SingletonType<TFloat64Range>::GetInstance().m;
    TFloat64Range &s2 = SingletonType<TFloat64Range>::GetInstance().m;
    s1 = a;
    BOOST_CHECK(SingletonType<TFloat64Range>::GetInstance().m.GetBegin() == a.GetBegin() );
    BOOST_CHECK(SingletonType<TFloat64Range>::GetInstance().m.GetEnd() == a.GetEnd() );

    BOOST_CHECK(s2.GetBegin() == a.GetBegin() );
    BOOST_CHECK(s2.GetEnd() == a.GetEnd() );

    s2 = b;
    BOOST_CHECK(SingletonType<TFloat64Range>::GetInstance().m.GetBegin() == b.GetBegin() );
    BOOST_CHECK(SingletonType<TFloat64Range>::GetInstance().m.GetEnd() == b.GetEnd() );

    BOOST_CHECK(s1.GetBegin() == b.GetBegin() );
    BOOST_CHECK(s1.GetEnd() == b.GetEnd() );

}

BOOST_AUTO_TEST_SUITE(singleton_test)

//-----------------------------------------------------------------------------
SingletonTest<Int32> testInt;
SingletonTest<Float64> testFloat;
SingletonTest<std::string> testString;
SingletonTest<TFloat64Range> testRange;
 
//test on instance creation
BOOST_AUTO_TEST_CASE(IsCreated_test)
{
    testInt.creation();
    testFloat.creation();
    testString.creation();
    testRange.creation();
}


BOOST_AUTO_TEST_CASE(run_test)
{
    testInt.run(1,2);
    testFloat.run(1.,2.);
    testString.run("1", "2");
    testRange.run(TFloat64Range(1.,2.), TFloat64Range(3.,4.));

}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
