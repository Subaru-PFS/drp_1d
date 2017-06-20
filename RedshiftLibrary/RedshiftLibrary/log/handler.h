#ifndef _CORE_LOG_HANDLER_
#define _CORE_LOG_HANDLER_

#include <RedshiftLibrary/common/datatypes.h>

namespace NSEpic
{

class CLog;

/**
 * \ingroup Core
 * Interface for log implementing custom log handler
 */
class CLogHandler
{

public:

    CLogHandler( CLog& log );
    virtual ~CLogHandler();

    Void    SetLevelMask( UInt32 mask );
    UInt32  GetLevelMask() const;

    virtual Void LogEntry( UInt32 lvl, const char* header, const char* msg ) = 0;

private:

    UInt32      m_LevelMask;
    CLog*       m_Logger;

};

}

#endif
