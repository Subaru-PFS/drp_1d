%module redshift
%include <std_string.i>
%include <std_shared_ptr.i>
%include typemaps.i
%shared_ptr(CParameterStore)
%shared_ptr(CLogConsoleHandler)
%shared_ptr(CClassifierStore)
%shared_ptr(CTemplateCatalog)
%shared_ptr(CRayCatalog)

%apply std::string &OUTPUT { std::string& out_str };
%apply Int64 &OUTPUT { Int64& out_int };
%apply Float64 &OUTPUT { Float64& out_float };


%{
	#include "RedshiftLibrary/log/log.h"
	#include "RedshiftLibrary/log/consolehandler.h"
	#include "RedshiftLibrary/processflow/parameterstore.h"
	#include "RedshiftLibrary/reliability/zclassifierstore.h"
	#include "RedshiftLibrary/processflow/context.h"
	#include "RedshiftLibrary/processflow/processflow.h"
	#include "RedshiftLibrary/ray/catalog.h"
  	#include "RedshiftLibrary/spectrum/template/catalog.h"
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
  %rename(Get_String) Get( const std::string& name, std::string& out_str, std::string = "");
  %rename(Get_Int64) Get( const std::string& name, Int64& out_int, Int64 defaultValue = 0);
  %rename(Get_Float64) Get( const std::string& name, Int64& out_int, Int64 defaultValue = 0);
public:
  CParameterStore();
  Bool Load( const std::string& path );
  Bool Save( const std::string& path ) const;
  Bool Get( const std::string& name, std::string& out_str, std::string defaultValue = "" );
  Bool Get( const std::string& name, Int64& out_int, Int64 defaultValue = 0 );
  Bool Get( const std::string& name, Float64& out_float, Float64 defaultValue  = 0 );

};

class CClassifierStore {
public:
  CClassifierStore();
  Bool Load ( const char* dirPath );
};

class CRayCatalog
{
public:
    Bool Load( const char* filePath );
};

class CTemplateCatalog
{
public:
    CTemplateCatalog( std::string cremovalmethod="Median", Float64 mediankernelsize=75.0, Float64 waveletsScales=8, std::string waveletsDFBinPath="");
    Bool Load( const char* filePath );
};

class CProcessFlowContext {
public:
  CProcessFlowContext();
  bool Init( const char* spectrumPath,
	     const char* noisePath,
	     std::string processingID,
	     std::shared_ptr<const CTemplateCatalog> templateCatalog,
	     std::shared_ptr<const CRayCatalog> rayCatalog,
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

