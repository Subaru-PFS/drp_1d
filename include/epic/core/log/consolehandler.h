#ifndef _CORE_LOG_CONSOLEHANDLER_
#define _CORE_LOG_CONSOLEHANDLER_

#include <epic/core/common/datatypes.h>
#include <epic/core/log/handler.h>

namespace NSEpic
{

/**
 * \ingroup Core
 * Handler that output log to stdout
 */
class CLogConsoleHandler : public CLogHandler
{

public:

    CLogConsoleHandler( CLog& logger );
    ~CLogConsoleHandler();

    Void LogEntry( UInt32 lvl, const char* header, const char* msg );

private:


};

}

#endif
