#ifndef _CORE_LOG_HANDLER_
#define _CORE_LOG_HANDLER_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/managedobject.h>

namespace NSEpic
{

class CLog;

/**
 * \ingroup Core
 * Interface for log implementing custom log handler
 */
class CLogHandler : public CManagedObject
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
