%module(directors="1") redshift

%include typemaps.i
%include std_string.i
%include std_shared_ptr.i
%include std_except.i

%shared_ptr(CClassifierStore)
%shared_ptr(CLogConsoleHandler)
%shared_ptr(CLogHandler)
%shared_ptr(CLog)
%shared_ptr(CSingleton<CLog>)
%shared_ptr(CParameterStore)
%shared_ptr(CRayCatalog)
%shared_ptr(CSpectrum)
%shared_ptr(CSpectrumAxis)
%shared_ptr(CSpectrumFluxAxis)
%shared_ptr(CSpectrumSpectralAxis)
%shared_ptr(CSpectrumIOGenericReader)
%shared_ptr(CSpectrumIOReader)
%shared_ptr(CTemplateCatalog)

%apply std::string &OUTPUT { std::string& out_str };
%apply NSEpic::Int64 &OUTPUT { NSEpic::Int64& out_int };
%apply NSEpic::Float64 &OUTPUT { NSEpic::Float64& out_float };

%feature("director");

%{
        #define SWIG_FILE_WITH_INIT
        #include "RedshiftLibrary/common/datatypes.h"
        #include "RedshiftLibrary/common/range.h"
        #include "RedshiftLibrary/log/log.h"
        #include <RedshiftLibrary/common/singleton.h>
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

%include "../RedshiftLibrary/RedshiftLibrary/common/datatypes.h"

using namespace NSEpic;

/* typedef int Int32; */
/* typedef short Int16; */
/* typedef signed char Int8 ; */
/* typedef long long Int64 ; */
/* typedef unsigned long long UInt64 ; */
/* typedef unsigned int UInt32 ; */
/* typedef unsigned short UInt16 ; */
/* typedef unsigned char UInt8 ; */
/* typedef float   Float32 ; */
/* typedef double  Float64 ; */
/* typedef char Char; */
/* typedef unsigned char Byte; */
/* typedef void Void; */
/* typedef unsigned int Bool; */
/* typedef const char* String; */

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

class CParameterStore {
%rename(Get_String) Get( const std::string& name, std::string& out_str, std::string = "");
%rename(Get_Int64) Get( const std::string& name, Int64& out_int, Int64 defaultValue = 0);
%rename(Get_Float64) Get( const std::string& name, Float64& out_float, Float64 defaultValue = 0);
%rename(Set_String) Set( const std::string& name, std::string& out_str, std::string = "");
public:
  CParameterStore();
  Bool Load( const std::string& path );
  Bool Save( const std::string& path ) const;
  Bool Get( const std::string& name, std::string& out_str, std::string defaultValue = "" );
  Bool Get( const std::string& name, Int64& out_int, Int64 defaultValue = 0 );
  Bool Get( const std::string& name, Float64& out_float, Float64 defaultValue  = 0 );
  Bool Set( const std::string& name, const std::string& v );

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
  void SaveRedshiftResult( const std::string& dir );
  void SaveReliabilityResult( const std::string& dir );
  void SaveAllResults( const std::string& dir, const std::string opt ) const;

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
  const Bool IsNoiseValid( Float64 LambdaMin,  Float64 LambdaMax ) const;
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

%apply (double* IN_ARRAY1, int DIM1) {(const Float64* samples, UInt32 n)};
class CSpectrumAxis
{
 public:
  // CSpectrumAxis(); // needs %rename
  CSpectrumAxis(const Float64* samples, UInt32 n );
  Float64* GetSamples();
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
//%rename(CSpectrumFluxAxis_default) CSpectrumFluxAxis(const Float64* samples, UInt32 n);
//%rename(CSpectrumFluxAxis_withError) CSpectrumFluxAxis( const Float64* _samples, UInt32 n,
// 							const Float64* _error, UInt32 m );
//%rename(CSpectrumFluxAxis_withError) CSpectrumFluxAxis( double* , int , double* , int );
%apply (double* IN_ARRAY1, int DIM1) {(const NSEpic::Float64* _samples, NSEpic::UInt32 n),
                                      (const NSEpic::Float64* _error, NSEpic::UInt32 m)};

class CSpectrumFluxAxis : public CSpectrumAxis
{
 public:
  // CSpectrumFluxAxis(); // needs %rename
  //CSpectrumFluxAxis( const Float64* samples, UInt32 n );
  CSpectrumFluxAxis( const Float64* _samples, UInt32 n, const Float64* _error, UInt32 m );
  void SetSize( UInt32 s );
};
//%clear (const Float64* samples, UInt32 n);
//%clear (const Float64* _samples, const Float64* _samples, UInt32 n);

