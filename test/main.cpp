#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/XmlOutputter.h>

#include <epic/core/log/log.h>
#include <epic/core/log/consolehandler.h>
#include <epic/core/common/typeinforegistry.h>
#include <epic/core/socket/context.h>

int main( int argc, char **argv)
{

    // Create the event manager and test controller
    CppUnit::TestResult controller;

    // Add a listener that collects test result
    CppUnit::TestResultCollector result;
    controller.addListener( &result );

    // Add a listener that print dots as test run.
    CppUnit::BriefTestProgressListener progress;
    controller.addListener( &progress );

    // Create logger and handler
    __NS__::CLog log;
    __NS__::CLogConsoleHandler consoleHandler( log );

    __NS__::CTypeInfoRegistry typeInfoRegistry;

    __NS__::CSocketContext socketContext;

    // Run tests
    CppUnit::TextUi::TestRunner runner;
    runner.addTest( CppUnit::TestFactoryRegistry::getRegistry().makeTest()  );
    runner.run( controller );

    // Print test in a compiler compatible format.
    CppUnit::CompilerOutputter outputter( &result, CppUnit::stdCOut() );
    outputter.write();

    // Uncomment this for XML output
    std::ofstream file( "../build/cppunit-reports/cppunit-report.xml" );
    CppUnit::XmlOutputter xml( &result, file );
    xml.write();
    file.close();

    return result.wasSuccessful() ? 0 : 1;
}

