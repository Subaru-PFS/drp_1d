
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/log/consolehandler.h>


using namespace NSEpic;

struct MyConfig
{

    CLog                                m_Log;
    CLogConsoleHandler                  m_ConsoleHandler;
    MyConfig():
        m_ConsoleHandler( m_Log )
    {

    }
    ~MyConfig()
    {

    }
};

BOOST_GLOBAL_FIXTURE( MyConfig )

BOOST_AUTO_TEST_SUITE(Common)

BOOST_AUTO_TEST_CASE(Median3)
{
  BOOST_CHECK(true);
}
BOOST_AUTO_TEST_SUITE_END()
