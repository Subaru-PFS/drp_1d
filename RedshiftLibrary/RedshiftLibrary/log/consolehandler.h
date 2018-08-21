#ifndef _CORE_LOG_CONSOLEHANDLER_
#define _CORE_LOG_CONSOLEHANDLER_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/log/handler.h>

namespace NSEpic
{

/**
 * \ingroup Core
 * Handler that outputs log to stdout.
 */
class CLogConsoleHandler : public CLogHandler
{

public:

    CLogConsoleHandler( CLog& logger );
    ~CLogConsoleHandler();

    void LogEntry( UInt32 lvl, const char* header, const char* msg );

private:


};

}

#endif
