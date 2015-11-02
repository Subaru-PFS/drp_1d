#include <epic/core/log/log.h>

#include <epic/core/log/handler.h>

#include <cstdio>

#define LOG_WORKING_BUFFER_SIZE 4096

using namespace NSEpic;

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

Void CLog::LogError( const char* format, ... )
{
    m_Mutex.Lock();

    va_list args;
    va_start( args, format );
    LogEntry( nLevel_Error, format, args );
    va_end( args );

    m_Mutex.Unlock();
}

Void CLog::LogWarning( const char* format, ... )
{
    m_Mutex.Lock();

    va_list args;
    va_start( args, format );
    LogEntry( nLevel_Warning, format, args );
    va_end( args );

    m_Mutex.Unlock();
}

Void CLog::LogInfo( const char* format, ... )
{
    m_Mutex.Lock();

    va_list args;
    va_start( args, format );
    LogEntry( nLevel_Info, format, args );
    va_end( args );

    m_Mutex.Unlock();
}

Void CLog::LogEntry( ELevel lvl, const char* format, va_list& args )
{
    Int32 i;

    vsnprintf( m_WorkingBuffer, LOG_WORKING_BUFFER_SIZE, format, args );

    for( i=0;i<LOG_HANDLER_TABLE_SIZE;i++ )
    {
        if( m_HandlerTable[i] )
        {
            if( m_HandlerTable[i]->GetLevelMask() & lvl )
                m_HandlerTable[i]->LogEntry( lvl, GetHeader( lvl ), m_WorkingBuffer );
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
