#include <epic/core/test/helper.h>
#include <epic/core/log/log.h>
#include <epic/core/debug/debug.h>

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/XmlOutputter.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/regex.hpp>

#include <sysexits.h>

using namespace NSEpic;

/**
 * Assignment constructor.
 */
CTestHelper::CTestHelper( int argc, char **argv ) :
    m_ConsoleHandler( m_Log )
{

}

/**
 * Empty destructor.
 */
CTestHelper::~CTestHelper()
{

}

/**
 * Calls signal handler installer, test options configurator and returns a call to ConfigureCmdLineOptions. 
 */
Bool CTestHelper::Configure( int argc, char **argv )
{
    InstallSignalHandler();
    ConfigureTestOptions();

    // Create log handler
    //CLog log;
    //CRef<CLogConsoleHandler> logConsoleHandler;
    //logConsoleHandler = new CLogConsoleHandler ( CLog::GetInstance() );
    //logConsoleHandler->SetLevelMask ( CLog::nLevel_Info );
    //logConsoleHandler->SetLevelMask ( CLog::nLevel_Debug );
    
    return ConfigureCmdLineOptions( argc, argv );
}

/**
 * Configure the test suite according to commandline options passed in.
 */
Bool CTestHelper::ConfigureCmdLineOptions( int argc, char **argv )
{
    // Declare the supported options.
    boost::program_options::options_description optionsDesc( "Apply various operation to an input spectrum" );
    optionsDesc.add_options()
      ( "help,h", "Print help" )
      ( "RegexFilter", boost::program_options::value< std::string>(), "Regular expression to run only a subset of the sepcified tests" );

    boost::program_options::positional_options_description positionalOptionsDesc;
    positionalOptionsDesc.add( "RegexFilter", 1 );

    boost::program_options::variables_map vm;
    try
    {
        boost::program_options::store( boost::program_options::command_line_parser( argc, argv ).options( optionsDesc ).positional( positionalOptionsDesc ).run(), vm );
        boost::program_options::notify( vm );
    }
    catch ( std::exception& e )
    {
        Log.LogError( "Error: %s", e.what() );
        return false;
    }

    if( vm.count( "RegexFilter" ) )
      {
        m_RegexFilter = vm["RegexFilter"].as< std::string >();
      }

    return true;
}

/**
 * 
 */
Void CTestHelper::ConfigureTestOptions()
{
    // Add a listener that collects test result
    m_TestResult.addListener( &m_TestResultCollector );
    // Add a listener that print dots as test run.
    m_TestResult.addListener( &m_TestProgressListener );
}

/**
 * 
 */
Bool CTestHelper::TestMatchFilter( CppUnit::Test* t )
{
    if( m_RegexFilter.size() == 0 )
        return false;

    boost::regex e( m_RegexFilter );
    return boost::regex_match( t->getName().c_str(), e );
}

/**
 * 
 */
Void CTestHelper::RecursiveFilterTest( CppUnit::Test* t )
{
    for( Int32 i=0; i<t->getChildTestCount(); i++ )
    {
        CppUnit::Test* child = t->getChildTestAt( i );

        if( child->getChildTestCount() )
        {
            RecursiveFilterTest( child );
        }
        else
        {
            if( TestMatchFilter( child ) )
            {
                m_TestRunner.addTest( child );
            }
        }
    }
}

/**
 * 
 */
Bool CTestHelper::Run( const char* outputFile )
{
    // Run tests
    if( m_RegexFilter.size() )
    {
        CppUnit::TestSuite* suite = (CppUnit::TestSuite*) CppUnit::TestFactoryRegistry::getRegistry().makeTest();
        RecursiveFilterTest( suite );
    }
    else
    {
        m_TestRunner.addTest( CppUnit::TestFactoryRegistry::getRegistry().makeTest() );
    }

    m_TestRunner.run( m_TestResult );

    // Print test in a compiler compatible format.
    CppUnit::CompilerOutputter outputter( &m_TestResultCollector, CppUnit::stdCOut() );
    outputter.write();

    // Uncomment this for XML output
    if( outputFile )
    {
        std::ofstream file( outputFile );
        CppUnit::XmlOutputter xml( &m_TestResultCollector, file );
        xml.write();
        file.close();
     }

    return m_TestResultCollector.wasSuccessful() ? 0 : 1;
}
