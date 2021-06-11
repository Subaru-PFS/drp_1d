#ifndef _REDSHIFT_LOG_FILEHANDLER_
#define _REDSHIFT_LOG_FILEHANDLER_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/log/handler.h>

#include <fstream>

namespace NSEpic
{

/**
 * \ingroup Redshift
 * Handler that output log to agiven file
 */
class CLogFileHandler : public CLogHandler
{

public:

    CLogFileHandler( const char* filePath );
    ~CLogFileHandler();

    void LogEntry( UInt32 lvl, const char* header, const char* msg );

private:

    std::fstream    m_OutputStream;
};

}

#endif
