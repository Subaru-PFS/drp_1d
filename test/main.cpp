
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <epic/core/log/log.h>
#include <epic/core/log/consolehandler.h>


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
