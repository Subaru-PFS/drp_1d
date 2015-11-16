#ifndef _CORE_LOG_LOG_
#define _CORE_LOG_LOG_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/singleton.h>
#include <epic/core/common/mutex.h>

#include <stdarg.h>

#define LOG_HANDLER_TABLE_SIZE 8
#define LOG_HANDLER_HEADER_LENGTH 64

#define Log (CLog::GetInstance())

namespace NSEpic
{

class CLogHandler;

/**
 * \ingroup Core
 * Log Class
 */
class CLog : public CSingleton< CLog >
{

public:

    enum ELevel
    {
        nLevel_Error = 1,
        nLevel_Warning,
        nLevel_Info,
        nLevel_All = 0x0fffffff,
        nLevel_Count,
        nLevel_None
    };

    CLog( );
    virtual ~CLog();

    Void LogError( const char* format, ... );
    Void LogWarning( const char* format, ... );
    Void LogInfo( const char* format, ... );

    CMutex& GetSynchMutex();

    Void Indent();
    Void UnIndent();

private:

    friend class CLogHandler;

    Void            LogEntry( ELevel lvl, const char*  format, va_list& args );
    Void            AddHandler( CLogHandler& handler );
    Void            RemoveHandler( CLogHandler& handler );
    const char*     GetHeader( CLog::ELevel lvl );

    //Attributes
    CLogHandler*     m_HandlerTable[LOG_HANDLER_TABLE_SIZE];
    Char            m_CurrentHeader[LOG_HANDLER_HEADER_LENGTH];
    Char*           m_WorkingBuffer;
    Char            m_IndentBuffer[128];
    Int32           m_IndentCount;

    CMutex          m_Mutex;


};


}

#endif
