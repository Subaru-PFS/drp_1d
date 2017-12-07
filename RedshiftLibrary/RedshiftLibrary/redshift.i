%module redshift
%include <std_string.i>
%include <std_shared_ptr.i>
%shared_ptr(CParameterStore)
%shared_ptr(CLogConsoleHandler)
%shared_ptr(CClassifierStore)

%{
	#include "RedshiftLibrary/log/log.h"
	#include "RedshiftLibrary/log/consolehandler.h"
	#include "RedshiftLibrary/processflow/parameterstore.h"
	#include "RedshiftLibrary/reliability/zclassifierstore.h"
	#include "RedshiftLibrary/processflow/context.h"
	#include "RedshiftLibrary/processflow/processflow.h"
	using namespace NSEpic;
%}

typedef int Int32;
typedef short Int16;
typedef signed char Int8 ;
typedef long long Int64 ;
typedef unsigned long long UInt64 ;
typedef unsigned int UInt32 ;
typedef unsigned short UInt16 ;
typedef unsigned char UInt8 ;
typedef float	Float32 ;
typedef double	Float64 ;
typedef char Char;
typedef unsigned char Byte;
typedef void Void;
typedef unsigned int Bool;
typedef const char* String;

class CLog {
public:
  enum ELevel   {
     nLevel_Critical = 100,
     nLevel_Error = 90,
     nLevel_Warning = 80,
     nLevel_Info = 70,
     nLevel_Debug = 60,
     nLevel_None = 0
   };
  CLog( );
};

class CLogConsoleHandler {
public:
  CLogConsoleHandler( CLog& logger );
  Void SetLevelMask( UInt32 mask );
};

class CParameterStore {
public:
  CParameterStore();
  Bool Load( const std::string& path );
  Bool Save( const std::string& path ) const;
};

class CClassifierStore {
public:
  CClassifierStore();
  Bool Load ( const char* dirPath );
};

class CProcessFlowContext {
public:
  CProcessFlowContext();
  bool Init( const char* spectrumPath,
	const char* noisePath,
	std::string processingID,
	const char* tempalteCatalogPath,
	const char* rayCatalogPath,
	std::shared_ptr<CParameterStore> paramStore,
        std::shared_ptr<CClassifierStore> zqualStore  );
  CDataStore& GetDataStore();
};

class CProcessFlow {
public:
  CProcessFlow();
  Bool Process( CProcessFlowContext& ctx );
};

class CDataStore
{
public:
  CDataStore( COperatorResultStore& resultStore, CParameterStore& parameStore );
  Void SaveRedshiftResult( const std::string& dir );
  Void SaveReliabilityResult( const std::string& dir );
  Void SaveAllResults( const std::string& dir ) const;

};

class COperatorResultStore
{
public:
  COperatorResultStore();
};



