%module(directors="1") redshift

%include typemaps.i
%include std_string.i
%include std_shared_ptr.i
%include std_except.i
%include std_vector.i
%include exception.i

%shared_ptr(CClassifierStore)
%shared_ptr(CLog)
%shared_ptr(CLogConsoleHandler)
%shared_ptr(CLogHandler)
%shared_ptr(CParameterStore)
%shared_ptr(COperatorResultStore)
%shared_ptr(CRayCatalog)
%shared_ptr(CLSF)
%shared_ptr(CLSFGaussianConstantWidth)
%shared_ptr(CLSFGaussianVariableWidth)
%shared_ptr(CSpectrum)
%shared_ptr(CSpectrumAxis)
%shared_ptr(CSpectrumFluxAxis)
%shared_ptr(CSpectrumIOGenericReader)
%shared_ptr(CSpectrumIOReader)
%shared_ptr(CSpectrumSpectralAxis)
%shared_ptr(CTemplate)
%shared_ptr(CTemplateCatalog)
%shared_ptr(CClassificationResult)
%shared_ptr(CPdfMargZLogResult)
%shared_ptr(TCandidateZ)
%shared_ptr(TExtremaResult)
%shared_ptr(TTplCombinationResult)
%shared_ptr(TLineModelResult)
%shared_ptr(CModelSpectrumResult)
%shared_ptr(CModelFittingResult)
%shared_ptr(CSpectraFluxResult)
%shared_ptr(TLSFArguments)
%shared_ptr(TLSFGaussianVarWidthArgs)
%shared_ptr(TLSFGaussianConstantWidthArgs)
%shared_ptr(TLSFGaussianConstantResolutionArgs)
%shared_ptr(TLSFGaussianNISPVSSPSF201707Args)
%shared_ptr(CLineModelSolution)


%feature("director");
%feature("nodirector") CSpectrumFluxAxis;

%{
#define SWIG_FILE_WITH_INIT
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/pyconv.h"
#include "RedshiftLibrary/version.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/common/singleton.h"
#include "RedshiftLibrary/log/consolehandler.h"
#include "RedshiftLibrary/log/filehandler.h"
#include "RedshiftLibrary/common/exception.h"
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
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/spectrum/LSFConstantWidth.h"
#include "RedshiftLibrary/method/solvedescription.h"
#include "RedshiftLibrary/method/classificationresult.h"
#include "RedshiftLibrary/operator/pdfMargZLogResult.h"
#include "RedshiftLibrary/statistics/pdfcandidatesz.h"
#include "RedshiftLibrary/statistics/pdfcandidateszresult.h"
#include "RedshiftLibrary/operator/extremaresult.h"
#include "RedshiftLibrary/linemodel/linemodelextremaresult.h"
#include "RedshiftLibrary/operator/tplCombinationExtremaResult.h"
#include "RedshiftLibrary/linemodel/modelfittingresult.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/operator/spectraFluxResult.h"
using namespace NSEpic;
static PyObject* pParameterException;
static PyObject* pGlobalException;
static PyObject* pSolveException;
static PyObject* pAmzException;
 %}

%include numpy.i

%init %{
  import_array();
  pParameterException = PyErr_NewException("redshift_.ParameterException", 0, 0);
  Py_INCREF(pParameterException);
  PyModule_AddObject(m, "ParameterException", pParameterException);
  pGlobalException = PyErr_NewException("redshift_.GlobalException", 0, 0);
  Py_INCREF(pGlobalException);
  PyModule_AddObject(m, "GlobalException", pGlobalException);
  pSolveException = PyErr_NewException("redshift_.SolveException", 0, 0);
  Py_INCREF(pSolveException);
  PyModule_AddObject(m, "SolveException", pSolveException);
  pAmzException = PyErr_NewException("redshift_.AmzException", 0, 0);
  Py_INCREF(pAmzException);
  PyModule_AddObject(m, "AmzException", pAmzException);

%}

%{

#define CATCH_PE(Exception) \
    catch(const Exception &e) \
    { \
       SWIG_Python_Raise(SWIG_NewPointerObj(new Exception(e), \
            SWIGTYPE_p_##Exception,SWIG_POINTER_OWN), \
            #Exception, SWIGTYPE_p_##Exception); \
       SWIG_fail; \
    } \
      
  /**/

// should be in "derived first" order
#define FOR_EACH_EXCEPTION(ACTION) \
   ACTION(ParameterException)       \
   ACTION(GlobalException) \
   ACTION(SolveException) \
   ACTION(AmzException) \
/**/
%}

%exception {
    try {
        $action
    }
    FOR_EACH_EXCEPTION(CATCH_PE)     
    catch (const std::exception & e)
    {
        SWIG_exception(SWIG_RuntimeError, (std::string("C++ std::exception: ") + e.what()).c_str());
    }
    catch (...)
    {
        SWIG_exception(SWIG_UnknownError, "C++ anonymous exception");
    }
}

%exceptionclass AmzException; 
// %include "../RedshiftLibrary/RedshiftLibrary/common/datatypes.h"
typedef double Float64;
typedef long long Int64;
typedef int Int32;
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

  static CLog& GetInstance();

  void LogError( const char* format, ... );
  void LogWarning( const char* format, ... );
  void LogInfo( const char* format, ... );
  void LogDetail( const char* format, ... );
  void LogDebug( const char* format, ... );
  void Indent();
  void UnIndent();

private:
   CLog();
   ~CLog();
};

class CLogConsoleHandler {
public:
  CLogConsoleHandler();
  void SetLevelMask( UInt32 mask );
};

class CLogFileHandler {
public:
  CLogFileHandler( const char* filePath );
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
typedef std::vector<std::string> TScopeStack;
typedef std::vector<Float64> TFloat64List;
typedef std::vector<Int32> TInt32List;

%template(TFloat64Range) CRange<Float64>;
%template(TFloat64List) std::vector<Float64>;
%template(TInt32List) std::vector<Int32>;

%apply std::string &OUTPUT { std::string& out_str };
%apply Int32 &OUTPUT { Int32& out_int };
%apply Int64 &OUTPUT { Int64& out_long };
%apply Float64 &OUTPUT { Float64& out_float };

class PC
{
 public:
  %rename(Get_Float64Array) get(const TFloat64List& vec,double ** ARGOUTVIEW_ARRAY1, int * DIM1);
  static void get(const TFloat64List& vec,double ** ARGOUTVIEW_ARRAY1, int * DIM1);
  %rename(Get_Int32Array) get(const TInt32List& vec,int ** ARGOUTVIEW_ARRAY1, int * DIM1);
  static void get(const TInt32List& vec,int ** ARGOUTVIEW_ARRAY1, int * DIM1);
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
    CTemplateCatalog( std::string cremovalmethod="Median", Float64 mediankernelsize=75.0, Float64 waveletsScales=8.0, std::string waveletsDFBinPath="" );
    void Add( std::shared_ptr<CTemplate> r);
};


class COperatorResult
{

public:

    COperatorResult();
    virtual ~COperatorResult();

    const std::string& getType();

};


%include "method/classificationresult.i"
%include "operator/pdfMargZLogResult.i"
%include "statistics/pdfcandidatesz.i"
%include "operator/extremaresult.i"
%include "operator/tplCombinationExtremaResult.i"
%include "linemodel/linemodelextremaresult.i"
%include "linemodel/modelfittingresult.i"
%include "operator/modelspectrumresult.i"
%include "linemodel/linemodelsolution.i"


class CSpectraFluxResult : public COperatorResult
{

public:

    CSpectraFluxResult();
    virtual ~CSpectraFluxResult();

    TFloat64List   fluxes;
    TFloat64List   wavel;

};

class CProcessFlowContext {
public:
  CProcessFlowContext();
  void Init(std::shared_ptr<CSpectrum> spectrum,
            std::shared_ptr<CTemplateCatalog> templateCatalog,
            std::shared_ptr<CRayCatalog> galaxy_rayCatalog,
            std::shared_ptr<CRayCatalog> qso_rayCatalog);
  std::shared_ptr<COperatorResultStore> GetResultStore();
  std::shared_ptr<const CParameterStore> LoadParameterStore(const std::string& paramsJSONString);
  void testResultStore();
};


class CParameterStore : public CScopeStore
{
  public:
    CParameterStore(const TScopeStack& stack);
    void FromString(const std::string& json);
    template<typename T> T Get(const std::string& name) const;
};

%template(Get_string) CParameterStore::Get<std::string>;

class CProcessFlow {
public:
  CProcessFlow();
  void Process( CProcessFlowContext& ctx );
};


class COperatorResultStore
{

 public:
  COperatorResultStore(const TScopeStack& scopeStack);

  std::shared_ptr<const CClassificationResult> GetClassificationResult(const std::string& objectType,
                                                                       const std::string& method,
                                                                       const std::string& name ) const;
  
  std::shared_ptr<const CPdfMargZLogResult> GetPdfMargZLogResult(const std::string& objectType,
								    const std::string& method,
								    const std::string& name ) const;

  std::shared_ptr<const TLineModelResult> GetLineModelResult(const std::string& objectType,
							     const std::string& method,
							     const std::string& name ,
							     const int& rank
							     ) const;
  std::shared_ptr<const TTplCombinationResult> GetTplCombinationResult(const std::string& objectType,
										 const std::string& method,
										 const std::string& name ,
										 const int& rank
										 ) const;
  std::shared_ptr<const TExtremaResult> GetExtremaResult(const std::string& objectType,
										 const std::string& method,
										 const std::string& name ,
										 const int& rank
									       ) const;

  std::shared_ptr<const CModelFittingResult> GetModelFittingResult(const std::string& objectType,
								   const std::string& method,
								   const std::string& name ,
								   const int& rank
								   ) const  ;

  std::shared_ptr<const CModelFittingResult> GetModelFittingResult(const std::string& objectType,
								   const std::string& method,
								   const std::string& name 
								   ) const  ;

  std::shared_ptr<const CLineModelSolution> GetLineModelSolution(const std::string& objectType,
								 const std::string& method,
								 const std::string& name 
								 ) const  ;

  std::shared_ptr<const CModelSpectrumResult> GetModelSpectrumResult(const std::string& objectType,
								     const std::string& method,
								     const std::string& name ,
								     const int& rank
								     ) const  ;

  std::shared_ptr<const CSpectraFluxResult> GetSpectraFluxResult(const std::string& objectType,
								   const std::string& method,
								   const std::string& name ,
								   const int& rank
								   ) const  ;
  
  const std::string&  GetGlobalResultType(const std::string& objectType,
                                          const std::string& method,
                                          const std::string& name ) const;

    const std::string&  GetCandidateResultType(const std::string& objectType,
					     const std::string& method,
					     const std::string& name ,
					     const std::string& dataset) const;

      bool HasCandidateDataset(const std::string& objectType,
			   const std::string& method,
			   const std::string& name ,
			   const std::string& dataset) const;

      bool HasDataset(const std::string& objectType,
			   const std::string& method,
			   const std::string& name ) const;

  int getNbRedshiftCandidates(const std::string& objectType,
			      const std::string& method) const;

  void test();

};

%catches(std::string, ...) CSpectrum::LoadSpectrum;

class CSpectrum
{
 %rename(CSpectrum_default) CSpectrum();
 public:
  CSpectrum();
  CSpectrum(CSpectrumSpectralAxis spectralAxis, CSpectrumFluxAxis fluxAxis);
  CSpectrum(CSpectrumSpectralAxis spectralAxis, CSpectrumFluxAxis fluxAxis, const std::shared_ptr<CLSF>& lsf);
  std::shared_ptr<const CLSF> GetLSF() const;
  void SetLSF(const std::shared_ptr<CLSF>& lsf);
  CSpectrumFluxAxis& GetFluxAxis();
  CSpectrumSpectralAxis& GetSpectralAxis();
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
%rename(CSpectrumFluxAxis_withError) CSpectrumFluxAxis( const double* samples, UInt32 n,
                                                        const double* error, UInt32 m );

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
	     CSpectrumSpectralAxis spectralAxis, CSpectrumFluxAxis fluxAxis);
  bool Save( const char* filePath ) const;
};

class CLSF
{
 enum TLSFType {
      GaussianConstantWidth,
      GaussianConstantResolution,
      GaussianNISPSIM2016,
      GaussianNISPVSSPSF201707,
      GaussianVariableWidth
    };
 public:
  CLSF(TLSFType name);
  virtual ~CLSF();
  virtual Float64 GetWidth(Float64 lambda) const=0;
  virtual bool IsValid() const=0;
};

class CLSFGaussianConstantWidth : public CLSF
{
 public:
  CLSFGaussianConstantWidth(const Float64 sigma=0.0);
  ~CLSFGaussianConstantWidth();
  Float64 GetWidth(Float64 lambda) const;
  bool IsValid() const;
};

class CLSFGaussianVariableWidth : public CLSF
{
 public:
  CLSFGaussianVariableWidth(const std::shared_ptr<const TLSFGaussianVarWidthArgs>& args);
  ~CLSFGaussianVariableWidth();
  Float64 GetWidth(Float64 lambda) const;
  bool IsValid() const;
};

class CLSFFactory : public CSingleton<CLSFFactory>
{
  public:
    std::shared_ptr<CLSF> Create(const std::string &name,
                                const std::shared_ptr<const TLSFArguments>& args);
    static CLSFFactory& GetInstance();                              
  private:
      friend class CSingleton<CLSFFactory>;

      CLSFFactory();
      ~CLSFFactory() = default;
};
  typedef enum ErrorCode
    {
      INVALID_SPECTRA_FLUX=	0,
      INVALID_NOISE	,
      SMALL_WAVELENGTH_RANGE ,
      NEGATIVE_CONTINUUMFIT	,
      BAD_CONTINUUMFIT	,
      NULL_AMPLITUDES	,
      PEAK_NOT_FOUND_PDF	,
      MAX_AT_BORDER_PDF	,
      UNKNOWN_PARAMETER  ,
      BAD_PARAMETER_VALUE,
      UNKNOWN_ATTRIBUTE ,
      BAD_LINECATALOG
    } ErrorCode;

typedef struct{
    virtual ~TLSFArguments(){};
    TLSFArguments()=default;
    TLSFArguments(const TLSFArguments & other) = default; 
    TLSFArguments(TLSFArguments && other) = default; 
    TLSFArguments& operator=(const TLSFArguments& other) = default;  
    TLSFArguments& operator=(TLSFArguments&& other) = default; 
}TLSFArguments;

struct TLSFGaussianVarWidthArgs : virtual TLSFArguments
{
  TFloat64List lambdas; 
  TFloat64List width;
  TLSFGaussianVarWidthArgs(TFloat64List _lambdas, TFloat64List _width);
};

struct TLSFGaussianConstantWidthArgs : virtual TLSFArguments
{
  Float64 width;
  TLSFGaussianConstantWidthArgs(const std::shared_ptr<const CParameterStore>& parameterStore);
};

struct TLSFGaussianConstantResolutionArgs : virtual TLSFArguments
{
  Float64 resolution;
  TLSFGaussianConstantResolutionArgs(const std::shared_ptr<const CParameterStore>& parameterStore);
};

struct TLSFGaussianNISPVSSPSF201707Args : virtual TLSFArguments
{
  Float64 sourcesize;
  TLSFGaussianNISPVSSPSF201707Args(const std::shared_ptr<const CParameterStore>& parameterStore);
};

class AmzException : public std::exception
{

 public:
  AmzException(ErrorCode ec,std::string message);
  ~AmzException();
 
  const char* getStackTrace() const;
  ErrorCode getErrorCode();
  virtual const char* what() ;
};


class GlobalException: public AmzException
{
 public:
  GlobalException(ErrorCode ec,std::string message);
  GlobalException(const GlobalException& e);
  ~GlobalException();
};


class SolveException: public AmzException
{
 public:
  SolveException(ErrorCode ec,std::string message);
  SolveException(const SolveException& e);
  ~SolveException();
};


class ParameterException: public AmzException
{
 public:
  ParameterException(ErrorCode ec,std::string message);
  ParameterException(const ParameterException& e);
  ~ParameterException();  
};

class CSolveDescription
{
 public:
  CSolveDescription(){}
  ~CSolveDescription(){}
    
  static const std::string GetDescription(const std::string& method);
};


