#ifndef _CORE_DEBUG_DEBUG_
#define _CORE_DEBUG_DEBUG_

#include <epic/core/common/datatypes.h>

#include <stdio.h>

namespace NSEpic
{

/**
 * \ingroup Core
 * @{
 */


// Internal usage should not be used directly
Void DebugLogCodeInformation( String file, UInt32 line, String func );
Void DebugLogInfo( String format, ... );
Void DebugLogWarning( String format, ... );
Void DebugLogError( String format, ... );

Void DebugBreakExecution();
Void DebugRaiseException();

Void DebugCreateDump( Void* dumpData = NULL ) ;

String DebugGetDumpDirectory();

Void InstallSignalHandler();
Void SignalHandler(int sig);

/**@}*/
}


#endif
