#ifndef _REDSHIFT_DEBUG_ASSERT_
#define _REDSHIFT_DEBUG_ASSERT_

#include <RedshiftLibrary/debug/debug.h>

#ifdef DEBUG_BUILD

    #define DebugAssertMsg(Exp, Msg)    \
        if( !(Exp) ) \
        { \
            NSEpic::DebugLogCodeInformation( __FILE__,__LINE__, __FUNCTION__ );\
            NSEpic::DebugLogInfo( "Assert: %s\n", #Exp ); \
            NSEpic::DebugLogInfo( Msg ); \
            NSEpic::DebugCreateDump( ); \
            NSEpic::DebugBreakExecution(); \
        } \

    #define DebugAssert( Exp )  \
        DebugAssertMsg( Exp, "Runtime assertion" ) \

    #define DebugError(Msg) \
        { \
            NSEpic::DebugLogCodeInformation( __FILE__,__LINE__, __FUNCTION__ );\
            NSEpic::DebugLogInfo( "Error:\n" ); \
            NSEpic::DebugLogInfo( #Msg ); \
            NSEpic::DebugCreateDump( ); \
            NSEpic::DebugBreakExecution(); \
        } \

#else

    #define DebugAssertMsg(Exp, Msg)    \
        if( !(Exp) ) \
        { \
            NSEpic::DebugCreateDump( ); \
            NSEpic::DebugRaiseException(); \
        } \

    #define DebugAssert( Exp )  \
        DebugAssertMsg( Exp, "Runtime assertion" ) \

    #define DebugError(Msg) \
        NSEpic::DebugCreateDump( ); \
        NSEpic::DebugRaiseException(); \

#endif


#endif
