%module(directors="1") redshift

%include typemaps.i
%include std_string.i
%include std_shared_ptr.i
%include std_except.i

%shared_ptr(CClassifierStore)
%shared_ptr(CLog)
%shared_ptr(CLogConsoleHandler)
%shared_ptr(CLogHandler)
%shared_ptr(CParameterStore)
%shared_ptr(CRayCatalog)
%shared_ptr(CSingleton<CLog>)
%shared_ptr(CSpectrum)
%shared_ptr(CSpectrumAxis)
%shared_ptr(CSpectrumFluxAxis)
%shared_ptr(CSpectrumIOGenericReader)
%shared_ptr(CSpectrumIOReader)
%shared_ptr(CSpectrumSpectralAxis)
%shared_ptr(CTemplate)
%shared_ptr(CTemplateCatalog)

%feature("director");
%feature("nodirector") CSpectrumFluxAxis;

%{
#define SWIG_FILE_WITH_INIT
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/version.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/common/singleton.h"
#include "RedshiftLibrary/log/consolehandler.h"
#include "RedshiftLibrary/log/filehandler.h"
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

// %include "../RedshiftLibrary/RedshiftLibrary/common/datatypes.h"
typedef	double Float64;
typedef unsigned int UInt32;

namespace NSEpic {
}

using namespace NSEpic;

const char* get_version();

class CLog {
public:
  enum ELevel   {
     nLevel_Critical = 100,
     nLevel_Error = 90,
     nLevel_Warning = 80,
     nLevel_Info = 70,
     nLevel_Detail = 65,
     nLevel_Debug = 60,
     nLevel_None = 0
   };
  CLog( );

  void LogError( const char* format, ... );
  void LogWarning( const char* format, ... );
  void LogInfo( const char* format, ... );
  void LogDetail( const char* format, ... );
  void LogDebug( const char* format, ... );
  void Indent();
  void UnIndent();
};

class CLogConsoleHandler {
public:
  CLogConsoleHandler( CLog& logger );
  void SetLevelMask( UInt32 mask );
};

class CLogFileHandler {
public:
  CLogFileHandler( CLog& logger, const char* filePath );
  void SetLevelMask( UInt32 mask );
};

template <typename T>
class CRange
{
 public:
  CRange( T begin, T end );
  const T& GetEnd() const;
  const T& GetBegin() const;
};
typedef CRange<Float64> TFloat64Range;
typedef TFloat64Range   TLambdaRange;
%template(TFloat64Range) CRange<Float64>;

%apply std::string &OUTPUT { std::string& out_str };
%apply Int64 &OUTPUT { Int64& out_int };
%apply Float64 &OUTPUT { Float64& out_float };

class CParameterStore {
%rename(Get_String) Get( const std::string& name, std::string& out_str, std::string = "");
%rename(Get_Int64) Get( const std::string& name, Int64& out_int, Int64 defaultValue = 0);
%rename(Get_Float64) Get( const std::string& name, Float64& out_float, Float64 defaultValue = 0);
%rename(Set_String) Set( const std::string& name, const std::string& v);
public:
  CParameterStore();
  void Load( const std::string& path );
  void Save( const std::string& path ) const;
  void Get( const std::string& name, std::string& out_str, std::string defaultValue = "" );
  void Get( const std::string& name, Int64& out_int, Int64 defaultValue = 0 );
  void Get( const std::string& name, Float64& out_float, Float64 defaultValue  = 0 );
  void Set( const std::string& name, const std::string& v );

};

class CClassifierStore {
public:
  CClassifierStore();
  bool Load ( const char* dirPath );
};

class CRayCatalog
{
public:
    void Load( const char* filePath );
    bool Save( const char* filePath );
    void ConvertVacuumToAir();
};

%catches(std::string, std::runtime_error, ...) CTemplateCatalog::Load;

class CTemplateCatalog
{
public:
    CTemplateCatalog( std::string cremovalmethod="Median", Float64 mediankernelsize=75.0, Float64 waveletsScales=8, std::string waveletsDFBinPath="");
    void Load( const char* filePath );
    void Add( std::shared_ptr<CTemplate> r );
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
  COperatorResultStore& GetResultStore();
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
  CDataStore(COperatorResultStore& resultStore, CParameterStore& parameStore);
  void SaveRedshiftResult(const std::string& dir);
  void SaveCandidatesResult(const std::string& dir);
  void SaveReliabilityResult(const std::string& dir);
  void SaveStellarResult(const std::string& dir);
  void SaveClassificationResult(const std::string& dir);
  void SaveAllResults(const std::string& dir, const std::string opt) const;
};

class COperatorResultStore
{
public:
  COperatorResultStore();

};

%catches(std::string, ...) CSpectrum::LoadSpectrum;

class CSpectrum
{
 %rename(CSpectrum_default) CSpectrum();
 public:
  CSpectrum();
  CSpectrum(CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis);
  CSpectrumFluxAxis& GetFluxAxis();
  CSpectrumSpectralAxis& GetSpectralAxis();
  void LoadSpectrum(const char* spectrumFilePath, const char* noiseFilePath);
  TLambdaRange GetLambdaRange() const;
  %apply Float64& OUTPUT { Float64& mean };
  %apply Float64& OUTPUT { Float64& std };

  void  SetName( const char* name );
  const std::string GetName() const;

  const bool IsNoiseValid( Float64 LambdaMin,  Float64 LambdaMax ) const;
  bool GetMeanAndStdFluxInRange(TFloat64Range wlRange,  Float64& mean, Float64 &std) const;
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

%rename(CSpectrumAxis_default) CSpectrumAxis();
%rename(CSpectrumAxis_empty) CSpectrumAxis(UInt32 n);
%rename(CSpectrumAxis_withSpectrum) CSpectrumAxis(const Float64* samples, UInt32 n);

%apply (double* IN_ARRAY1, int DIM1) {(const Float64* samples, UInt32 n)};
class CSpectrumAxis
{
 public:
  CSpectrumAxis();
  CSpectrumAxis( UInt32 n );
  CSpectrumAxis(const Float64* samples, UInt32 n );
  Float64* GetSamples();
  UInt32 GetSamplesCount() const;
  virtual void SetSize( UInt32 s );
};
%clear (const Float64* samples, UInt32 n);

%apply (double* IN_ARRAY1, int DIM1) {(const Float64* samples, UInt32 n)};
class CSpectrumSpectralAxis : public CSpectrumAxis {
 public:
  // CSpectrumSpectralAxis(); // needs %rename
  CSpectrumSpectralAxis( const Float64* samples, UInt32 n );
};
%clear (const Float64* samples, UInt32 n);


//%apply (double* IN_ARRAY1, int DIM1) {(const Float64* samples, UInt32 n)};

%rename(CSpectrumFluxAxis_default) CSpectrumFluxAxis();
%rename(CSpectrumFluxAxis_empty) CSpectrumFluxAxis(UInt32 n);
%rename(CSpectrumFluxAxis_withSpectrum) CSpectrumFluxAxis(const Float64* samples, UInt32 n);
%rename(CSpectrumFluxAxis_withError) CSpectrumFluxAxis( double* samples,
							UInt32 n,
 							double* error,
							UInt32 m );

%apply (double* IN_ARRAY1, int DIM1) {(const Float64* samples, UInt32 n)}
%apply (double* IN_ARRAY1, int DIM1) {(const double* samples, UInt32 n),
                                      (const double* error, UInt32 m)}

class CSpectrumFluxAxis : public CSpectrumAxis
{
 public:
  // CSpectrumFluxAxis(); // needs %rename
  CSpectrumFluxAxis();
  CSpectrumFluxAxis( UInt32 n );
  CSpectrumFluxAxis( const Float64* samples, UInt32 n );
  CSpectrumFluxAxis( const double* samples, UInt32 n,
  		     const double* error, UInt32 m );
  void SetSize( UInt32 s );
};
//%clear (const Float64* samples, UInt32 n);
//%clear (const Float64* _samples, const Float64* _samples, UInt32 n);

class CTemplate : public CSpectrum
{
 public:
  CTemplate( const std::string& name, const std::string& category,
	     CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis);
  bool Save( const char* filePath ) const;
};

