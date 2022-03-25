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
#include "RedshiftLibrary/debug/debug.h"
#include "RedshiftLibrary/log/log.h"

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
void NSEpic::DebugLogCodeInformation( String file, Int32 line, String func )
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
