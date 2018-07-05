%module(directors="1") redshift

%include typemaps.i
%include <std_string.i>
%include <std_shared_ptr.i>
%include std_except.i

%shared_ptr(CParameterStore)
%shared_ptr(CLogConsoleHandler)
%shared_ptr(CClassifierStore)
%shared_ptr(CTemplateCatalog)
%shared_ptr(CRayCatalog)
%shared_ptr(CSpectrum)
%shared_ptr(CSpectrumIOReader)
%shared_ptr(CSpectrumIOGenericReader)
%feature("director");

%apply std::string &OUTPUT { std::string& out_str };
%apply Int64 &OUTPUT { Int64& out_int };
%apply Float64 &OUTPUT { Float64& out_float };

%{
        #define SWIG_FILE_WITH_INIT
        #include "RedshiftLibrary/log/log.h"
        #include "RedshiftLibrary/log/consolehandler.h"
        #include "RedshiftLibrary/processflow/parameterstore.h"
        #include "RedshiftLibrary/reliability/zclassifierstore.h"
        #include "RedshiftLibrary/processflow/context.h"
        #include "RedshiftLibrary/processflow/processflow.h"
        #include "RedshiftLibrary/processflow/resultstore.h"
        #include "RedshiftLibrary/ray/catalog.h"
        #include "RedshiftLibrary/spectrum/template/catalog.h"
        #include "RedshiftLibrary/spectrum/io/reader.h"
        #include "RedshiftLibrary/spectrum/io/genericreader.h"
        #include "RedshiftLibrary/spectrum/axis.h"
        #include "RedshiftLibrary/spectrum/fluxaxis.h"
        #include "RedshiftLibrary/spectrum/spectralaxis.h"
        using namespace NSEpic;
%}

%include numpy.i

%init %{
import_array();
%}

typedef int Int32;
typedef short Int16;
typedef signed char Int8 ;
typedef long long Int64 ;
typedef unsigned long long UInt64 ;
typedef unsigned int UInt32 ;
typedef unsigned short UInt16 ;
typedef unsigned char UInt8 ;
typedef float   Float32 ;
typedef double  Float64 ;
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
    void Load( const char* filePath );
};

%catches(std::string, std::runtime_error, ...) CTemplateCatalog::Load;

class CTemplateCatalog
{
public:
    CTemplateCatalog( std::string cremovalmethod="Median", Float64 mediankernelsize=75.0, Float64 waveletsScales=8, std::string waveletsDFBinPath="");
    void Load( const char* filePath );
};

class CProcessFlowContext {
public:
  CProcessFlowContext();
  bool Init(std::shared_ptr<CSpectrum> spectrum,
	    std::string processingID,
	    std::shared_ptr<const CTemplateCatalog> templateCatalog,
	    std::shared_ptr<const CRayCatalog> rayCatalog,
	    std::shared_ptr<CParameterStore> paramStore,
	    std::shared_ptr<CClassifierStore> zqualStore  );
  CDataStore& GetDataStore();
};

%catches(std::string, std::runtime_error, ...) CProcessFlow::Process;

class CProcessFlow {
public:
  CProcessFlow();
  void Process( CProcessFlowContext& ctx );
};

class CDataStore
{
public:
  CDataStore( COperatorResultStore& resultStore, CParameterStore& parameStore );
  Void SaveRedshiftResult( const std::string& dir );
  Void SaveReliabilityResult( const std::string& dir );
  Void SaveAllResults( const std::string& dir, const std::string opt ) const;

};

class COperatorResultStore
{
public:
  COperatorResultStore();
};

%catches(std::string, ...) CSpectrum::LoadSpectrum;

class CSpectrum
{
 public:
  CSpectrum();
  CSpectrumFluxAxis& GetFluxAxis();
  CSpectrumSpectralAxis& GetSpectralAxis();
  void LoadSpectrum(const char* spectrumFilePath, const char* noiseFilePath);

};

%catches(std::string, ...) CSpectrumIOReader::Read;

class CSpectrumIOReader
{
 public:
  CSpectrumIOReader();
  virtual ~CSpectrumIOReader();
  virtual void Read(const char* filePath, CSpectrum& s) = 0;
};

class CSpectrumIOGenericReader : public CSpectrumIOReader
{
 public:
  CSpectrumIOGenericReader();
  virtual ~CSpectrumIOGenericReader();
  virtual void Read( const char* filePath, CSpectrum& s );
};

class CSpectrumAxis
{
 public:
  CSpectrumAxis();
  Float64* GetSamples();
  virtual void SetSize( UInt32 s );
};

class CSpectrumSpectralAxis : public CSpectrumAxis {
 public:
  CSpectrumSpectralAxis();
%apply (double* IN_ARRAY1, int DIM1) {( const Float64* samples, UInt32 n)};
  CSpectrumSpectralAxis( const Float64* samples, UInt32 n, Bool isLogScale  );
%clear (double* IN_ARRAY1, int DIM1);
};


class CSpectrumFluxAxis : public CSpectrumAxis
{
 public:
  CSpectrumFluxAxis();
%apply (double* IN_ARRAY1, int DIM1) {( const Float64* samples, UInt32 n)};
  CSpectrumFluxAxis( const Float64* samples, UInt32 n );
%clear (double* IN_ARRAY1, int DIM1);
  void SetSize( UInt32 s );
};
