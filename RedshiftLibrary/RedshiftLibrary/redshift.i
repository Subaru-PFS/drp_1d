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
%module(directors="1") redshift

%include typemaps.i
%include std_string.i
%include std_shared_ptr.i
%include std_except.i
%include std_vector.i
%include std_map.i
%include exception.i

%shared_ptr(CClassifierStore)
%shared_ptr(CLog)
%shared_ptr(CLogConsoleHandler)
%shared_ptr(CLogHandler)
%shared_ptr(CParameterStore)
%shared_ptr(COperatorResultStore)
%shared_ptr(CLineCatalog)
%shared_ptr(CLSF)
%shared_ptr(CLSFGaussianConstantWidth)
%shared_ptr(CLSFGaussianVariableWidth)
%shared_ptr(CSpectrum)
%shared_ptr(CSpectrumAxis)
%shared_ptr(CSpectrumFluxAxis)
%shared_ptr(CSpectrumNoiseAxis)
%shared_ptr(CSpectrumSpectralAxis)
%shared_ptr(CTemplate)
%shared_ptr(CTemplateCatalog)
%shared_ptr(CClassificationResult)
%shared_ptr(CReliabilityResult)
%shared_ptr(CPdfMargZLogResult)
%shared_ptr(TCandidateZ)
%shared_ptr(TExtremaResult)
%shared_ptr(TTplCombinationResult)
%shared_ptr(TLineModelResult)
%shared_ptr(CModelSpectrumResult)
%shared_ptr(CSpectraFluxResult)
%shared_ptr(TLSFArguments)
%shared_ptr(TLSFGaussianVarWidthArgs)
%shared_ptr(TLSFGaussianConstantWidthArgs)
%shared_ptr(TLSFGaussianConstantResolutionArgs)
%shared_ptr(TLSFGaussianNISPVSSPSF201707Args)
%shared_ptr(CLineModelSolution)
%shared_ptr(CPhotometricData)
%shared_ptr(CPhotometricBand)
%shared_ptr(std::map<std::string, CPhotometricBand>) // needed for CPhotBandCatalog (the base classes in the hierarchy must be declared as shared_ptr as well)
%shared_ptr(CPhotBandCatalog)
%shared_ptr(CLineCatalogsTplShape)
%shared_ptr(CLineRatioCatalog)
%shared_ptr(CFlagLogResult)
%shared_ptr(CFlagWarning)
%shared_ptr(CSpectrumFluxCorrectionMeiksin)
%shared_ptr(CSpectrumFluxCorrectionCalzetti) 
%feature("director");
%feature("nodirector") CSpectrumFluxAxis;

%{
#define SWIG_FILE_WITH_INIT
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/pyconv.h"
#include "RedshiftLibrary/version.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/common/singleton.h"
#include "RedshiftLibrary/log/consolehandler.h"
#include "RedshiftLibrary/log/filehandler.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/processflow/resultstore.h"
#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/line/lineRatioCatalog.h"
#include "RedshiftLibrary/line/catalogsTplShape.h"
#include "RedshiftLibrary/line/lineprofile.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/spectrum/fluxaxis.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"
#include "RedshiftLibrary/method/solvedescription.h"
#include "RedshiftLibrary/method/classificationresult.h"
#include "RedshiftLibrary/method/reliabilityresult.h"
#include "RedshiftLibrary/operator/pdfMargZLogResult.h"
#include "RedshiftLibrary/operator/flagResult.h"
#include "RedshiftLibrary/statistics/pdfcandidatesz.h"
#include "RedshiftLibrary/statistics/pdfcandidateszresult.h"
#include "RedshiftLibrary/operator/extremaresult.h"
#include "RedshiftLibrary/linemodel/linemodelextremaresult.h"
#include "RedshiftLibrary/operator/tplCombinationExtremaResult.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/operator/spectraFluxResult.h"
#include "RedshiftLibrary/photometry/photometricdata.h"
#include "RedshiftLibrary/photometry/photometricband.h"
#include "RedshiftLibrary/method/linemodelsolve.h"
#include "RedshiftLibrary/method/templatefittingsolve.h"
#include "RedshiftLibrary/method/linemeassolve.h"
#include "RedshiftLibrary/method/tplcombinationsolve.h"
#include "RedshiftLibrary/method/reliabilitysolve.h"
#include "RedshiftLibrary/method/classificationsolve.h"
#include "RedshiftLibrary/method/linematchingsolve.h"
#include "RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"

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

class CFlagWarning {
public:
  typedef enum WarningCode
    {
      WARNING_NONE=0,
      AIR_VACCUM_CONVERSION_IGNORED, //1
      CRANGE_VALUE_OUTSIDERANGE, //2
      CRANGE_VECTBORDERS_OUTSIDERANGE, //3
      CRANGE_NO_INTERSECTION, //4
      FINDER_NO_PEAKS, //5
      STDESTIMATION_NO_MATCHING, //6
      STDESTIMATION_FAILED, //7
      MULTIROLL_STRTAG_NOTFOUND, //8
      LINEMATCHING_REACHED_ENDLOOP, //9
      FORCE_LOGSAMPLING_FFT,//10
      IGNORELINESSUPPORT_DISABLED_FFT,//11
      FORCE_FROMSPECTRUM_NEG_CONTINUUMAMP,//12
      INVALID_MERIT_VALUES,//13
      AIR_VACCUM_REACHED_MAX_ITERATIONS, //14
      ASYMFIT_NAN_PARAMS,//15
      DELTAZ_COMPUTATION_FAILED, //16
      INVALID_FOLDER_PATH,//17
      TPL_NAME_EMPTY,//18
      RELIABILITY_NEEDS_TENSORFLOW//19
    }WarningCode;

  static CFlagWarning& GetInstance();

  void warning(WarningCode c, std::string message);
  void warning(WarningCode c, const char* format, ... );
  Int32 getBitMask();
  void resetFlag();

private:
  CFlagWarning();
  ~CFlagWarning();
};

class CLogConsoleHandler {
public:
  CLogConsoleHandler();
  void SetLevelMask( Int32 mask );
};

class CLogFileHandler {
public:
  CLogFileHandler( const char* filePath );
  void SetLevelMask( Int32 mask );
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
typedef std::vector<std::string> TStringList;


%template(TFloat64Range) CRange<Float64>;
%template(TFloat64List) std::vector<Float64>;
%template(TInt32List) std::vector<Int32>;
%template(TStringList) std::vector<std::string>;
%template(VecTFloat64List) std::vector<  std::vector<Float64> >;

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
  %rename(Get_AxisSampleList) getasl(const TAxisSampleList& vec,double ** ARGOUTVIEW_ARRAY1, int * DIM1);
  static void getasl(const TAxisSampleList& vec,double ** ARGOUTVIEW_ARRAY1, int * DIM1);

};

class CLineCatalog
{
public:
  CLineCatalog();
  CLineCatalog(Float64 nSigmaSupport);
  void AddLineFromParams(const std::string& name,
			const Float64& position,
			const std::string& type,
			const std::string& force,
			const std::string& profile,
			const TAsymParams& asymParams,
			const std::string& groupName,
			const Float64& nominalAmplitude,
			const std::string& velocityGroup,
			const Float64& velocityOffset,
			const bool& enableVelocityFit,
			const Int32& id,
			const std::string& str_id);

  void setLineAmplitude(const std::string& str_id,const Float64& nominalAmplitude);
  void setAsymProfileAndParams(const std::string& profile, TAsymParams params);


};

typedef struct {
  TAsymParams(Float64 sigma,Float64 alpha,Float64 delta);
  Float64 sigma, alpha, delta;
} TAsymParams;

class CLineRatioCatalog : public CLineCatalog
{
 public:
  CLineRatioCatalog(const std::string& name, const CLineCatalog& lineCatalog);
  ~CLineRatioCatalog();
  void addVelocity(const std::string& name, Float64 value);
  void setPrior(Float64 prior);

  void setIsmIndex(Float64 ismIndex);

};

class CLineCatalogsTplShape
{

public:
  CLineCatalogsTplShape();
  void addLineRatioCatalog(const CLineRatioCatalog &lr_catalog);

};

%catches(std::string, std::runtime_error, ...) CTemplateCatalog::Load;

class CTemplateCatalog
{
public:
    CTemplateCatalog();
    void Add( std::shared_ptr<CTemplate> r);
};


class COperatorResult
{

public:

    COperatorResult();
    virtual ~COperatorResult();

    const std::string& getType();

};

%template(TMapFloat64) std::map<std::string, Float64>;

%include "method/classificationresult.i"
%include "method/reliabilityresult.i"
%include "operator/pdfMargZLogResult.i"
%include "operator/flagResult.i"
%include "statistics/pdfcandidatesz.i"
%include "operator/extremaresult.i"
%include "operator/tplCombinationExtremaResult.i"
%include "linemodel/linemodelextremaresult.i"
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

  static CProcessFlowContext& GetInstance();

  void Init();
  void setLineCatalog(const std::string& objectType, const std::string& method, const std::shared_ptr<CLineCatalog> &catalog); 
  void setLineRatioCatalogCatalog(const std::string& objectType, const std::shared_ptr<CLineCatalogsTplShape> &catalog); 
  void setTemplateCatalog(const std::shared_ptr<CTemplateCatalog> &templateCatalog){ m_TemplateCatalog = templateCatalog;}
  void setPhotBandCatalog(const std::shared_ptr<CPhotBandCatalog> &photBandCatalog){ m_photBandCatalog = photBandCatalog;}
  void setSpectrum(const std::shared_ptr<CSpectrum> &spectrum){ m_Spectrum = spectrum;}
  void setfluxCorrectionMeiksin(const std::shared_ptr<CSpectrumFluxCorrectionMeiksin> &igmcorrectionMeiksin){m_igmcorrectionMeiksin = igmcorrectionMeiksin;}
  void setfluxCorrectionCalzetti(const std::shared_ptr<CSpectrumFluxCorrectionCalzetti> &ismcorrectionCalzetti){m_ismcorrectionCalzetti = ismcorrectionCalzetti;}
  void reset();


  const std::shared_ptr<COperatorResultStore> &GetResultStore();
  std::shared_ptr<const CParameterStore> LoadParameterStore(const std::string& paramsJSONString);
 
  TScopeStack  m_ScopeStack;

 private:
    CProcessFlowContext();
    ~CProcessFlowContext();
    
};



class CParameterStore : public CScopeStore
{
  public:
    CParameterStore(const TScopeStack& stack);
    void FromString(const std::string& json);
    template<typename T> T Get(const std::string& name) const;
};

%template(Get_string) CParameterStore::Get<std::string>;

class COperatorResultStore
{

 public:
  COperatorResultStore(const TScopeStack& scopeStack);

  std::shared_ptr<const CClassificationResult> GetClassificationResult(const std::string& objectType,
                                                                       const std::string& method,
                                                                       const std::string& name ) const;

  std::shared_ptr<const CReliabilityResult> GetReliabilityResult(const std::string& objectType,
									 const std::string& method,
                                                                       const std::string& name ) const;

  std::shared_ptr<const CPdfMargZLogResult> GetPdfMargZLogResult(const std::string& objectType,
								    const std::string& method,
								    const std::string& name ) const;

  std::shared_ptr<const CFlagLogResult> GetFlagResult(const std::string& objectType,
								    const std::string& method,
								    const std::string& name ) const;

  std::shared_ptr<const TLineModelResult> GetLineModelResult(const std::string& objectType,
							     const std::string& method,
							     const std::string& name ,
							     const int& rank,
                   bool firstpassResults
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


  std::shared_ptr<const CLineModelSolution> GetLineModelSolution(const std::string& objectType,
								 const std::string& method,
								 const std::string& name,
								     const int& rank 
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

  std::shared_ptr<const CModelSpectrumResult> GetModelSpectrumResult(const std::string& objectType,
								     const std::string& method,
								     const std::string& name 
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

  void StoreFlagResult( const std::string& name, Int32  result );
};


class CSpectrum
{
 %rename(CSpectrum_default) CSpectrum();
 public:
  CSpectrum();
  CSpectrum(CSpectrumSpectralAxis spectralAxis, CSpectrumFluxAxis fluxAxis);
  CSpectrum(CSpectrumSpectralAxis spectralAxis, CSpectrumFluxAxis fluxAxis, const std::shared_ptr<const CLSF>& lsf);
  std::shared_ptr<const CLSF> GetLSF() const;
  void SetLSF(const std::shared_ptr<const CLSF>& lsf);
  void SetPhotData(const std::shared_ptr<const CPhotometricData>& photData);
  CSpectrumFluxAxis& GetFluxAxis();
  CSpectrumSpectralAxis& GetSpectralAxis();
  const CSpectrumNoiseAxis&  GetErrorAxis() const;
  TLambdaRange GetLambdaRange() const;
  %apply Float64& OUTPUT { Float64& mean };
  %apply Float64& OUTPUT { Float64& std };

  void  SetName( const char* name );
  const std::string GetName() const;

  const bool IsNoiseValid( Float64 LambdaMin,  Float64 LambdaMax ) const;
  bool GetMeanAndStdFluxInRange(TFloat64Range wlRange,  Float64& mean, Float64 &std) const;
};


%rename(CSpectrumAxis_default) CSpectrumAxis();
%rename(CSpectrumAxis_empty) CSpectrumAxis(Int32 n);
%rename(CSpectrumAxis_withSpectrum) CSpectrumAxis(const Float64* samples, Int32 n);

%apply (double* IN_ARRAY1, int DIM1) {(const Float64* samples, Int32 n)};
class CSpectrumAxis
{
 public:
  CSpectrumAxis();
  CSpectrumAxis( Int32 n );
  CSpectrumAxis(const Float64* samples, Int32 n );
  Float64* GetSamples();
  const TAxisSampleList& GetSamplesVector() const;
  Int32 GetSamplesCount() const;
  virtual void SetSize( Int32 s );
};
%clear (const Float64* samples, Int32 n);

%apply (double* IN_ARRAY1, int DIM1) {(const Float64* samples, Int32 n)};
class CSpectrumSpectralAxis : public CSpectrumAxis {
 public:
  // CSpectrumSpectralAxis(); // needs %rename
  CSpectrumSpectralAxis( const Float64* samples, Int32 n, std::string AirVacuum="" );
};
%clear (const Float64* samples, Int32 n);

%rename(CSpectrumFluxAxis_default) CSpectrumFluxAxis();
%rename(CSpectrumFluxAxis_empty) CSpectrumFluxAxis(Int32 n);
%rename(CSpectrumFluxAxis_withSpectrum) CSpectrumFluxAxis(const Float64* samples, Int32 n);
%rename(CSpectrumFluxAxis_withError) CSpectrumFluxAxis( const double* samples, Int32 n,
                                                        const double* error, Int32 m );

%apply (double* IN_ARRAY1, int DIM1) {(const Float64* samples, Int32 n)}
%apply (double* IN_ARRAY1, int DIM1) {(const double* samples, Int32 n),
                                      (const double* error, Int32 m)}

class CSpectrumFluxAxis : public CSpectrumAxis
{
 public:
  // CSpectrumFluxAxis(); // needs %rename
  CSpectrumFluxAxis();
  CSpectrumFluxAxis( Int32 n );
  CSpectrumFluxAxis( const Float64* samples, Int32 n );
  CSpectrumFluxAxis( const double* samples, Int32 n,
  		     const double* error, Int32 m );
  void SetSize( Int32 s );

};

class  CSpectrumNoiseAxis : public CSpectrumAxis
{
 public:
  CSpectrumNoiseAxis();
};

//%clear (const Float64* samples, Int32 n);
//%clear (const Float64* _samples, const Float64* _samples, Int32 n);

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
  virtual ~CLSF();
  virtual Float64 GetWidth(Float64 lambda) const=0;
  virtual bool IsValid() const=0;
protected:
  CLSF();
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

class CPhotometricData
{
public:
    CPhotometricData(const TStringList & name,  const TFloat64List & flux, const TFloat64List & fluxerr);

    Float64 GetFlux(const std::string &name) const;
    Float64 GetFluxErr(const std::string &name) const;
    Float64 GetFluxOverErr2(const std::string & name) const;
    TStringList GetNameList() const;
};

%apply (double* IN_ARRAY1, int DIM1) {(const Float64 * trans, Int32 n1),(const Float64 * lambda, Int32 n2 )}
class CPhotometricBand
{
public:
    CPhotometricBand() = default;
    CPhotometricBand(const Float64 * trans, Int32 n1, const Float64 * lambda, Int32 n2  ); //for swig binding to numpy array
    const TFloat64List & GetTransmission() const;
    const TFloat64List & GetWavelength() const;
    Float64 GetMinLambda() const;
    Float64 IntegrateFlux(const TFloat64List &flux) const;
};
%clear (const Float64 * trans, Int32 n1);
%clear (const Float64 * lambda, Int32 n2);

%template(TMapPhotBand) std::map<std::string, CPhotometricBand>;

class CPhotBandCatalog : public std::map<std::string, CPhotometricBand> 
{
public:
    void Add(const std::string & name, const CPhotometricBand & filter);
    TStringList GetNameList() const;
    TStringList GetNameListSortedByLambda() const;
};
 //struct ErrorCode {
  typedef enum ErrorCode{
    //enum{
      INTERNAL_ERROR=0,
      EXTERNAL_LIB_ERROR,
      INVALID_SPECTRA_FLUX,
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
      BAD_LINECATALOG,
      BAD_LOGSAMPLEDSPECTRUM,
      BAD_COUNTMATCH,
      BAD_TEMPLATECATALOG,
      INVALID_SPECTRUM,
      OVERLAPRATE_NOTACCEPTABLE,
      DZ_NOT_COMPUTABLE
    }ErrorCode;
//}; 
class AmzException : public std::exception
{

 public: 
  AmzException(ErrorCode ec,std::string message);
  virtual ~AmzException();
 
  const char* getStackTrace() const;
  ErrorCode getErrorCode();
  virtual const char* what() ;
};


class GlobalException: public AmzException
{
 public:
  GlobalException(ErrorCode ec,std::string message);
  GlobalException(const GlobalException& e);
  virtual ~GlobalException();
};


class SolveException: public AmzException
{
 public:
  SolveException(ErrorCode ec,std::string message);
  SolveException(const SolveException& e);
  virtual ~SolveException();
};


class ParameterException: public AmzException
{
 public:
  ParameterException(ErrorCode ec,std::string message);
  ParameterException(const ParameterException& e);
  virtual ~ParameterException();  
};

class CSolveDescription
{
 public:
  CSolveDescription(){}
  ~CSolveDescription(){}
    
  static const std::string GetDescription(const std::string& method);
};

class CSolve{
 public:
  CSolve()=delete;
    void Compute();
};

class CObjectSolve{
 public:
  CSolve()=delete;
    void Compute();
};

  class CClassificationSolve:public CSolve
  {

  public:

    CClassificationSolve(TScopeStack &scope,std::string objectType);

  };
  class CReliabilitySolve:public CSolve
  {

  public:

    CReliabilitySolve(TScopeStack &scope,std::string objectType);
  };
  class CLineModelSolve:public CObjectSolve
  {

  public:

    CLineModelSolve(TScopeStack &scope,std::string objectType);
  };
  class CLineMeasSolve:public CObjectSolve
  {

  public:

    CLineMeasSolve(TScopeStack &scope,std::string objectType);
  };

  class CTemplateFittingSolve : public CObjectSolve
{
  public:

    CTemplateFittingSolve(TScopeStack &scope,std::string objectType);
  };

class CTplcombinationSolve : public CObjectSolve
{

 public:
  CTplcombinationSolve(TScopeStack &scope,std::string objectType);
};

class CLineMatchingSolve: public CObjectSolve
{
public:

    CLineMatchingSolve(TScopeStack &scope,std::string objectType);
};

typedef struct {
  MeiksinCorrection(TFloat64List _lbda, std::vector<TFloat64List> _fluxcorr);
  TFloat64List lbda; // wavelength
  std::vector<TFloat64List > fluxcorr; // 7 flux correction lists
}MeiksinCorrection;

%template(VecMeiksinCorrection) std::vector< MeiksinCorrection >;
class CSpectrumFluxCorrectionMeiksin
{
  public:
    CSpectrumFluxCorrectionMeiksin(std::vector<MeiksinCorrection> meiksinCorrectionCurves);
};

typedef struct {
  CalzettiCorrection(TFloat64List _lbda, TFloat64List _fluxcorr);
  TFloat64List lbda;
  TFloat64List fluxcorr; 
}CalzettiCorrection;
class CSpectrumFluxCorrectionCalzetti
{
  public:
    CSpectrumFluxCorrectionCalzetti(CalzettiCorrection _calzettiCorr, Float64 ebmv_start, Float64 ebmv_step, Float64 ebmv_n);
};