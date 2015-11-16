#include <epic/core/debug/debug.h>

#include <epic/core/log/log.h>

#include <assert.h>
#include <signal.h>
#include <stdio.h>
#include <execinfo.h>
#include <stdlib.h>
#include <unistd.h>

using namespace NSEpic;

Void NSEpic::DebugLogInfo( String format, ... )
{
    va_list args;
    va_start( args, format );
    Log.LogInfo(  format, args );
    va_end( args );
}

Void NSEpic::DebugLogWarning( String format, ... )
{
    va_list args;
    va_start( args, format );
    Log.LogWarning( format, args );
    va_end( args );
}

Void NSEpic::DebugLogError( String format, ... )
{
    va_list args;
    va_start( args, format );
    DebugCreateDump( );
    Log.LogError( format, args );
    va_end( args );
}

Void NSEpic::DebugLogCodeInformation( String file, UInt32 line, String func )
{
    Log.LogInfo( "%s()\n%s:%d", func, file, line );
}

Void NSEpic::DebugRaiseException()
{
    raise( SIGABRT );
}

Void NSEpic::DebugBreakExecution()
{
    raise( SIGABRT );
}

String NSEpic::DebugGetDumpDirectory()
{
    static String dir = "./CrashDump/";
    return dir;

}

Void NSEpic::DebugCreateDump( Void* dumpData )
{

}

Void NSEpic::InstallSignalHandler()
{
    signal( SIGSEGV, SignalHandler );
}

Void NSEpic::SignalHandler(int sig)
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
