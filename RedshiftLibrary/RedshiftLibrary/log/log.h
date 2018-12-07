#ifndef _CORE_LOG_LOG_
#define _CORE_LOG_LOG_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/singleton.h>
#include <RedshiftLibrary/common/mutex.h>

#include <stdarg.h>

#define LOG_HANDLER_TABLE_SIZE 8
#define LOG_HANDLER_HEADER_LENGTH 64

#define Log (CLog::GetInstance())

namespace NSEpic
{

class CLogHandler;

/**
 * \ingroup Core
 * Responsible for logging features.
 * This class allows prioritized logging, and selective output of these log messages according to handler's configuration.
 */
class CLog : public CSingleton< CLog >
{

public:

  enum ELevel
  {
    nLevel_Critical = 100,
    nLevel_Error = 90,
    nLevel_Warning = 80,
    nLevel_Info = 70,
    nLevel_Detail = 65,
    nLevel_Debug = 60,
    nLevel_None = 0
  };

    CLog( );
    virtual ~CLog();

    void LogError( const char* format, ... );
    void LogWarning( const char* format, ... );
    void LogInfo( const char* format, ... );
    void LogDetail( const char* format, ... );
    void LogDebug( const char* format, ... );

    CMutex& GetSynchMutex();

    void Indent();
    void UnIndent();

private:

    friend class CLogHandler;

    void            LogEntry( ELevel lvl, const char*  format, va_list& args );
    void            AddHandler( CLogHandler& handler );
    void            RemoveHandler( CLogHandler& handler );
    const char*     GetHeader( CLog::ELevel lvl );

    //Attributes
    CLogHandler*    m_HandlerTable[LOG_HANDLER_TABLE_SIZE];
    Char            m_CurrentHeader[LOG_HANDLER_HEADER_LENGTH];
    Char*           m_WorkingBuffer;
    Char            m_IndentBuffer[128];
    Int32           m_IndentCount;

    CMutex          m_Mutex;


};


}

#endif
