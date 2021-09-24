// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include "RedshiftLibrary/log/log.h"

#include "RedshiftLibrary/log/handler.h"

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
void CLog::LogError( const char* format, ... )
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
void CLog::LogWarning( const char* format, ... )
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
void CLog::LogInfo( const char* format, ... )
{
    m_Mutex.Lock();

    va_list args;
    va_start( args, format );
    LogEntry( nLevel_Info, format, args );
    va_end( args );

    m_Mutex.Unlock();
}

/**
 * Calls LogEntry with an Information message.
 */
void CLog::LogDetail( const char* format, ... )
{
    m_Mutex.Lock();

    va_list args;
    va_start( args, format );
    LogEntry( nLevel_Detail, format, args );
    va_end( args );

    m_Mutex.Unlock();
}

/**
 * Calls LogEntry with a Debug message.
 */
void CLog::LogDebug( const char* format, ... )
{
    m_Mutex.Lock();

    va_list args;
    va_start( args, format );
    LogEntry( nLevel_Debug, format, args );
    va_end( args );

    m_Mutex.Unlock();
}

void CLog::log(const std::string& msg,CLog::ELevel lvl)
{
   m_Mutex.Lock();
   for( int i=0;i<LOG_HANDLER_TABLE_SIZE;i++ )
     {
       if( m_HandlerTable[i] )
         {
           if ( m_HandlerTable[i]->GetLevelMask() <= lvl )
             {
               m_HandlerTable[i]->LogEntry( lvl, GetHeader( lvl ), msg.c_str() );
             }
         }
     }
   m_Mutex.Unlock();
}

void CLog::LogError( const std::string& msg){log(msg,nLevel_Error);}
void CLog::LogDebug( const std::string& msg){log(msg,nLevel_Debug);}
void CLog::LogInfo( const std::string& msg){log(msg,nLevel_Info);}
void CLog::LogWarning( const std::string& msg){log(msg,nLevel_Warning);}
void CLog::LogDetail( const std::string& msg){log(msg,nLevel_Detail);}


/**
 * Calls LogEntry in every handler registered on the table with the input message, if the handler has a level mask <= the message level.
 */
void CLog::LogEntry( ELevel lvl, const char* format, va_list& args )
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

void CLog::Indent()
{

    m_IndentBuffer[m_IndentCount] = '\t';
    m_IndentCount++;
    m_IndentBuffer[m_IndentCount] = 0;
}
void CLog::UnIndent()
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
void CLog::RemoveHandler( CLogHandler& handler )
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
void CLog::AddHandler( CLogHandler& handler )
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
    case nLevel_Detail:
        sprintf( m_CurrentHeader, "Detail: %s ", m_IndentBuffer  );
        break;
    case nLevel_Debug:
        sprintf( m_CurrentHeader, "Debug: %s", m_IndentBuffer );
        break;
    case nLevel_None:
        sprintf( m_CurrentHeader, "%s", m_IndentBuffer );
        break;
    default:
        return NULL;
    }

    return m_CurrentHeader;
}
