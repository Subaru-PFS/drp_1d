#include "RedshiftLibrary/log/consolehandler.h"

#include <stdio.h>

using namespace NSEpic;



/**
 * Prints log message in stdout.
 */
void CLogConsoleHandler::LogEntry( UInt32 lvl, const char* header, const char* msg )
{
    if( header )
        printf("%s", header );

    printf("%s\n", msg );
}
