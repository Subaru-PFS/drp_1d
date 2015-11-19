#include <epic/core/log/consolehandler.h>

#include <stdio.h>

using namespace NSEpic;


CLogConsoleHandler::CLogConsoleHandler( CLog& logger ) :
    CLogHandler( logger )
{
}

CLogConsoleHandler::~CLogConsoleHandler()
{

}

Void CLogConsoleHandler::LogEntry( UInt32 lvl, const char* header, const char* msg )
{
    if( header )
        printf("%s", header );

    printf("%s\n", msg );
}
