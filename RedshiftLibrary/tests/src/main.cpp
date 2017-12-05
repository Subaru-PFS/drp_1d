#define BOOST_TEST_MODULE Redshift Determination Library test suite
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

BOOST_GLOBAL_FIXTURE( MyConfig );
