#ifndef _CORE_LOG_FILEHANDLER_
#define _CORE_LOG_FILEHANDLER_

#include <epic/core/common/datatypes.h>
#include <epic/core/log/handler.h>

#include <fstream>

namespace NSEpic
{

/**
 * \ingroup Core
 * Handler that output log to agiven file
 */
class CLogFileHandler : public CLogHandler
{

    DEFINE_MANAGED_OBJECT( CLogFileHandler )

public:

    CLogFileHandler( CLog& logger, const char* filePath );
    ~CLogFileHandler();

    Void LogEntry( UInt32 lvl, const char* header, const char* msg );

private:

    std::fstream    m_OutputStream;
};

}

#endif
