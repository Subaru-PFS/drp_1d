#ifndef _CORE_DEBUG_DEBUG_
#define _CORE_DEBUG_DEBUG_

#include <RedshiftLibrary/common/datatypes.h>

#include <stdio.h>

namespace NSEpic
{

/**
 * \ingroup Core
 * @{
 */


// Internal usage should not be used directly
void DebugLogCodeInformation( String file, UInt32 line, String func );
void DebugLogInfo( String format, ... );
void DebugLogWarning( String format, ... );
void DebugLogError( String format, ... );

void DebugBreakExecution();
void DebugRaiseException();

void DebugCreateDump( void* dumpData = NULL ) ;

String DebugGetDumpDirectory();

void InstallSignalHandler();
void SignalHandler(int sig);

/**@}*/
}


#endif
