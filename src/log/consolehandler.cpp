#include <epic/core/log/consolehandler.h>

#include <stdio.h>

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT_NOT_INSTANCIABLE( CLogConsoleHandler );

CLogConsoleHandler::CLogConsoleHandler( CLog& logger ) :
    CLogHandler( logger )
{
}

CLogConsoleHandler::~CLogConsoleHandler()
{

}

/**
 * Prints log message in stdout.
 */
Void CLogConsoleHandler::LogEntry( UInt32 lvl, const char* header, const char* msg )
{
    if( header )
        printf("%s", header );

    printf("%s\n", msg );
}
