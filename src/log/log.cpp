#include <epic/core/log/log.h>

#include <epic/core/log/handler.h>

#include <cstdio>

#define LOG_WORKING_BUFFER_SIZE 4096

using namespace NSEpic;

/**
 * Creates a singleton that holds a buffer and a table of the log handlers.
 */
CLog::CLog( )
{
    Int32 i;
    m_IndentCount = 0;

    m_WorkingBuffer = new Char[LOG_WORKING_BUFFER_SIZE];
    for( i=0;i<LOG_HANDLER_TABLE_SIZE;i++ )
    {
        m_HandlerTable[i] = NULL;
    }

    m_IndentBuffer[0] = 0;
}

CLog::~CLog()
{
    delete [] m_WorkingBuffer;
}

/**
 * Calls LogEntry with an Error message.
 */
Void CLog::LogError( const char* format, ... )
{
    m_Mutex.Lock();

    va_list args;
    va_start( args, format );
    LogEntry( nLevel_Error, format, args );
    va_end( args );

    m_Mutex.Unlock();
}

/**
 * Calls LogEntry with a Warning message.
 */
Void CLog::LogWarning( const char* format, ... )
{
    m_Mutex.Lock();

    va_list args;
    va_start( args, format );
    LogEntry( nLevel_Warning, format, args );
    va_end( args );

    m_Mutex.Unlock();
}

/**
 * Calls LogEntry with an Information message.
 */
Void CLog::LogInfo( const char* format, ... )
{
    m_Mutex.Lock();

    va_list args;
    va_start( args, format );
    LogEntry( nLevel_Info, format, args );
    va_end( args );

    m_Mutex.Unlock();
}

/**
 * Calls LogEntry with a Debug message.
 */
Void CLog::LogDebug( const char* format, ... )
{
    m_Mutex.Lock();
    va_list args;
    va_start( args, format );
    LogEntry( nLevel_Debug, format, args );
    va_end( args );

    m_Mutex.Unlock();
}

/**
 * Calls LogEntry in every handler registered on the table with the input message, if the handler has a level mask <= the message level.
 */
Void CLog::LogEntry( ELevel lvl, const char* format, va_list& args )
{
  Int32 i;

  vsnprintf( m_WorkingBuffer, LOG_WORKING_BUFFER_SIZE, format, args );

  for( i=0;i<LOG_HANDLER_TABLE_SIZE;i++ )
    {
      if( m_HandlerTable[i] )
        {
	  if ( m_HandlerTable[i]->GetLevelMask() <= lvl )
	    {
	      m_HandlerTable[i]->LogEntry( lvl, GetHeader( lvl ), m_WorkingBuffer );
	    }
        }
    }
}

CMutex& CLog::GetSynchMutex()
{
    return m_Mutex;
}

Void CLog::Indent()
{

    m_IndentBuffer[m_IndentCount] = '\t';
    m_IndentCount++;
    m_IndentBuffer[m_IndentCount] = 0;
}
Void CLog::UnIndent()
{
    if( m_IndentCount == 0 )
        return;

    m_IndentBuffer[m_IndentCount] = 0;
    m_IndentCount--;
    m_IndentBuffer[m_IndentCount] = 0;
}

/**
 * Remove the reference for the input handler if it exists in the handler table.
 */
Void CLog::RemoveHandler( CLogHandler& handler )
{
    Int32 i;

    for( i=0; i<LOG_HANDLER_TABLE_SIZE; i++ )
    {
        if( m_HandlerTable[i] == &handler )
        {
            m_HandlerTable[i] = NULL;
            break;
        }
    }
}

/**
 * Insert a reference to the input handler in the first vacant position in the handler table.
 * Will _not_ prevent a handler to be included multiple times.
 */
Void CLog::AddHandler( CLogHandler& handler )
{
    Int32 i;

    for( i=0; i<LOG_HANDLER_TABLE_SIZE; i++ )
    {
        if( m_HandlerTable[i] == NULL )
        {
            m_HandlerTable[i] = & handler;
            break;
        }
    }
}

/**
 * Given a message priority level, returns an appropriate prefix to the message.
 */
const char* CLog::GetHeader( CLog::ELevel lvl )
{
    switch( lvl )
    {
    case nLevel_Error:
        sprintf( m_CurrentHeader, "Error: %s ", m_IndentBuffer );
        break;
    case nLevel_Warning:
        sprintf( m_CurrentHeader, "Warning: %s ", m_IndentBuffer  );
        break;
    case nLevel_Info:
        sprintf( m_CurrentHeader, "Info: %s ", m_IndentBuffer  );
        break;
    case nLevel_None:
        sprintf( m_CurrentHeader,  "%s", m_IndentBuffer );
        break;
    default:
        return NULL;
    }

    return m_CurrentHeader;
}
