#include <RedshiftLibrary/debug/debug.h>
#include <RedshiftLibrary/log/log.h>

#include <assert.h>
#include <signal.h>
#include <stdio.h>
#include <execinfo.h>
#include <stdlib.h>
#include <unistd.h>

using namespace NSEpic;

/**
 * Wrapper around Log.LogInfo.
 */
void NSEpic::DebugLogInfo( String format, ... )
{
    va_list args;
    va_start( args, format );
    Log.LogInfo( format, args );
    va_end( args );
}

/**
 * Wrapper around Log.LogWarning.
 */
void NSEpic::DebugLogWarning( String format, ... )
{
    va_list args;
    va_start( args, format );
    Log.LogWarning( format, args );
    va_end( args );
}

/**
 * Wrapper around Log.LogError.
 */
void NSEpic::DebugLogError( String format, ... )
{
    va_list args;
    va_start( args, format );
    DebugCreateDump( );
    Log.LogError( format, args );
    va_end( args );
}

/**
 * Traceback print.
 */
void NSEpic::DebugLogCodeInformation( String file, UInt32 line, String func )
{
    Log.LogInfo( "%s()\n%s:%d", func, file, line );
}

/**
 * Raises SIGABRT.
 */
void NSEpic::DebugRaiseException()
{
    raise( SIGABRT );
}

/**
 * Raises SIGABRT.
 */
void NSEpic::DebugBreakExecution()
{
    raise( SIGABRT );
}

/**
 * Returns a string containing "./CrashDump/".
 */
String NSEpic::DebugGetDumpDirectory()
{
    static String dir = "./CrashDump/";
    return dir;

}

/**
 * Empty.
 */
void NSEpic::DebugCreateDump( void* dumpData )
{

}

/**
 * Calls signal with a SignalHandler parameter.
 */
void NSEpic::InstallSignalHandler()
{
    signal( SIGSEGV, SignalHandler );
}

/**
 * Helper structure to pass arguments to signal.
 */
void NSEpic::SignalHandler(int sig)
{
  void *array[16];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace( array, 10 );

  // print out all the frames to stderr
  fprintf( stderr, "Error: signal %d:\n", sig );
  backtrace_symbols_fd( array, size, STDERR_FILENO );
  exit( 1 );
}
