#ifndef _CORE_TEST_HELPER_
#define _CORE_TEST_HELPER_

#include <epic/core/common/datatypes.h>

#include <epic/core/log/log.h>
#include <epic/core/log/consolehandler.h>
#include <epic/core/common/typeinforegistry.h>

#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/ui/text/TestRunner.h>

#include <string>
#include <vector>

namespace NSEpic
{

/**
 * \ingroup Core
 * Wrapper class for creating unit test applications
 *
 *
 */
class CTestHelper
{

public:

    CTestHelper( int argc, char **argv );
    ~CTestHelper();

    Bool    Run( const char* outputFile );
    Bool    Configure( int argc, char **argv );

private:

    Bool                ConfigureCmdLineOptions( int argc, char **argv );
    Void                ConfigureTestOptions();

    Bool                TestMatchFilter( CppUnit::Test* t );

    Void                RecursiveFilterTest( CppUnit::Test* t );

    CLog                                m_Log;
    CLogConsoleHandler                  m_ConsoleHandler;
    CTypeInfoRegistry                   m_TypeInfoRegistry;


    CppUnit::TestResult                 m_TestResult;
    CppUnit::TestResultCollector        m_TestResultCollector;
    CppUnit::BriefTestProgressListener  m_TestProgressListener;
    CppUnit::TextUi::TestRunner         m_TestRunner;

    std::string                         m_RegexFilter;
    std::vector<CppUnit::Test*>         m_FilteredTest;
};

}

#endif
