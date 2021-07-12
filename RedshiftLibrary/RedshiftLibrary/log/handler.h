#ifndef _REDSHIFT_LOG_HANDLER_
#define _REDSHIFT_LOG_HANDLER_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/log/log.h"
namespace NSEpic
{


/**
 * \ingroup Redshift
 * Interface for log implementing custom log handler
 */
class CLogHandler
{

public:

    CLogHandler();
    virtual ~CLogHandler();

    void    SetLevelMask( UInt32 mask );
    UInt32  GetLevelMask() const;

    virtual void LogEntry( UInt32 lvl, const char* header, const char* msg ) = 0;

private:

    UInt32      m_LevelMask;
    CLog     &  m_Logger;

};

}

#endif
