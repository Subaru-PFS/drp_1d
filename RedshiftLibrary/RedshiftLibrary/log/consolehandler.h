#ifndef _REDSHIFT_LOG_CONSOLEHANDLER_
#define _REDSHIFT_LOG_CONSOLEHANDLER_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/log/handler.h"

namespace NSEpic
{

/**
 * \ingroup Redshift
 * Handler that outputs log to stdout.
 */
class CLogConsoleHandler : public CLogHandler
{

public:

    void LogEntry( UInt32 lvl, const char* header, const char* msg );

private:


};

}

#endif
