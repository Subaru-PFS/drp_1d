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
%shared_ptr(CScopeStack)
%shared_ptr(CParameterStore)
%shared_ptr(COperatorResultStore)
%shared_ptr(COperatorResult)
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
%shared_ptr(CLogZPdfResult)
%shared_ptr(TCandidateZ)
%shared_ptr(TExtremaResult)
%shared_ptr(TTplCombinationResult)
%shared_ptr(TLineModelResult)
%shared_ptr(CModelSpectrumResult)
%shared_ptr(CModelPhotValueResult)
%shared_ptr(TLSFArguments)
%shared_ptr(TLSFGaussianVarWidthArgs)
%shared_ptr(TLSFGaussianConstantWidthArgs)
%shared_ptr(TLSFGaussianConstantResolutionArgs)
%shared_ptr(TLSFGaussianNISPVSSPSF201707Args)
%shared_ptr(CLineModelSolution)
%shared_ptr(CLineModelSolveResult)
%shared_ptr(CPhotometricData)
%shared_ptr(CPhotometricBand)
%shared_ptr(std::map<std::string, CPhotometricBand>) // needed for CPhotBandCatalog (the base classes in the hierarchy must be declared as shared_ptr as well)
%shared_ptr(CPhotBandCatalog)
%shared_ptr(CLineCatalogsTplRatio)
%shared_ptr(CLineRatioCatalog)
%shared_ptr(CFlagLogResult)
%shared_ptr(CFlagWarning)
%shared_ptr(CSpectrumFluxCorrectionMeiksin)
%shared_ptr(CSpectrumFluxCorrectionCalzetti) 
%shared_ptr(TZgridListParams)
%feature("director");
%feature("nodirector") CSpectrumFluxAxis;

%{
#define SWIG_FILE_WITH_INIT
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/zgridparam.h"
#include "RedshiftLibrary/version.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/common/singleton.h"
#include "RedshiftLibrary/log/consolehandler.h"
#include "RedshiftLibrary/log/filehandler.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/processflow/resultstore.h"
#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/line/lineRatioCatalog.h"
#include "RedshiftLibrary/line/catalogsTplRatio.h"
#include "RedshiftLibrary/line/lineprofile.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/spectrum/fluxaxis.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"
#include "RedshiftLibrary/method/classificationresult.h"
#include "RedshiftLibrary/method/reliabilityresult.h"
#include "RedshiftLibrary/method/linemodelsolveresult.h"
#include "RedshiftLibrary/operator/logZPdfResult.h"
#include "RedshiftLibrary/operator/flagResult.h"
#include "RedshiftLibrary/statistics/pdfcandidatesz.h"
#include "RedshiftLibrary/statistics/pdfcandidateszresult.h"
#include "RedshiftLibrary/operator/extremaresult.h"
#include "RedshiftLibrary/linemodel/linemodelextremaresult.h"
#include "RedshiftLibrary/operator/tplCombinationExtremaResult.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/operator/modelphotvalueresult.h"
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
static PyObject* pAmzException;
 %}

%include numpy.i

%init %{
  import_array();
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

typedef Float64 Sample;
typedef unsigned char UInt8;
typedef UInt8 Mask;

const char* get_version();

static const std::string undefStr;

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

  void LogInfo(const std::string &s);
  void LogDetail(const std::string &s);
  void LogDebug(const std::string &s);

private:
   CLog();
   ~CLog();
};

%pythonprepend CFlagWarning::warning(WarningCode c, const std::string &message) %{
    c = c.bit
%} 

%pythonappend CFlagWarning::getBitMask() %{
    val = WarningCode(val)
%}

class CFlagWarning {
public:

  static CFlagWarning& GetInstance();

  void warning(WarningCode c, const std::string &message);
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

class CMask {
 public:
  explicit CMask(Int32 weightsCount, Int32 defaultValue = 0);
  //Mask &operator[](const Int32 i);
  const TMaskList& getMaskList() const;
};
%extend CMask {
   Mask __getitem__(unsigned int i) {
      return (*($self))[i];
    }
    Mask __setitem__(unsigned int i, int m) {
      (*($self))[i] = m;
      return (*($self))[i];
    }
}



typedef CRange<Float64> TFloat64Range;
typedef TFloat64Range   TLambdaRange;
typedef std::vector<std::string> TScopeStack;
typedef std::vector<Float64> TFloat64List;
typedef std::vector<Float32> TFloat32List;
typedef std::vector<Int32> TInt32List;
typedef std::vector<std::string> TStringList;
typedef std::vector<Mask> TMaskList;
typedef std::vector<bool> TBoolList;
typedef std::vector<Sample> TAxisSampleList;
typedef std::map<std::string, TFloat64List> TMapTFloat64List;

%template(TFloat64Range) CRange<Float64>;
%template(TFloat64List) std::vector<Float64>;
%template(TInt32List) std::vector<Int32>;
%template(TStringList) std::vector<std::string>;
%template(VecTFloat64List) std::vector<  std::vector<Float64> >;
%template(TBoolList) std::vector<bool>;
%template(TMaskList) std::vector<Mask>;
%template(TMapTFloat64List) std::map<std::string, TFloat64List>;

%apply std::string &OUTPUT { std::string& out_str };
%apply Int32 &OUTPUT { Int32& out_int };
%apply Int64 &OUTPUT { Int64& out_long };
%apply Float64 &OUTPUT { Float64& out_float };
%apply bool &OUTPUT { bool& out_bool };

%inline %{
  template <typename Tvec, typename Tarr>
  void stdVectorToNumpy(const std::vector<Tvec> &vec, Tarr **ARGOUTVIEWM_ARRAY1, int *DIM1) {
    size_t const byte_size = sizeof(Tarr) * vec.size();
    *ARGOUTVIEWM_ARRAY1 = (Tarr *) malloc(byte_size);
    for (std::size_t i = 0; i < vec.size(); i++)
      (*ARGOUTVIEWM_ARRAY1)[i] = vec[i];
    *DIM1 = vec.size();
  }

  template  <typename Tarr>
  void stdVectorToNumpy(const std::vector<Tarr> &vec, Tarr **ARGOUTVIEWM_ARRAY1, int *DIM1) {
    size_t const byte_size = sizeof(Tarr) * vec.size();
    *ARGOUTVIEWM_ARRAY1 = (Tarr *) malloc(byte_size);
      std::memcpy(*ARGOUTVIEWM_ARRAY1, vec.data(), byte_size);
    *DIM1 = vec.size();
  }
  
  template <typename Tarr>
  void stdMapVectorToNumpy(const std::map<std::string, std::vector<Tarr>> &map, const std::string &key, 
                        Tarr **ARGOUTVIEWM_ARRAY1, int *DIM1) {
    stdVectorToNumpy<Tarr>(map.at(key), ARGOUTVIEWM_ARRAY1, DIM1);                          
  }
%}

%template(stdVectorToNumpy) stdVectorToNumpy<Float64>;
%template(stdVectorToNumpy) stdVectorToNumpy<Int32>;
%template(stdVectorToNumpy) stdVectorToNumpy<Mask>;
// convert c++ bool type to short (numpy_typemaps for c++ bool not implemented)
%template(stdVectorToNumpy) stdVectorToNumpy<bool, short>;

%template(stdMapVectorToNumpy) stdMapVectorToNumpy<Float64>;

%define VectorExtendToNumpy(name)
%pythoncode %{
  name.to_numpy = lambda self: stdVectorToNumpy(self)
%}
%enddef

VectorExtendToNumpy(TFloat64List)
VectorExtendToNumpy(TInt32List)
VectorExtendToNumpy(TMaskList)
VectorExtendToNumpy(TBoolList)

%pythoncode %{
  TMapTFloat64List.to_numpy = lambda self, key: stdMapVectorToNumpy(self, key)
%}

class CLineCatalog : public CLineCatalogBase<>
{
public:
  CLineCatalog();
  CLineCatalog(Float64 nSigmaSupport);
  void AddLineFromParams(const std::string& name,
			Float64 position,
			const std::string& type,
			const std::string& force,
			const std::string& profile,
			const TAsymParams& asymParams,
			const std::string& groupName,
			Float64 nominalAmplitude,
			const std::string& velocityGroup,
			Float64 velocityOffset,
			bool enableVelocityFit,
			Int32 line_id,
      const std::string& strId,
      const std::shared_ptr<CSpectrumFluxCorrectionMeiksin>& igmcorrection=nullptr);

  void setLineAmplitude(Int32 line_id, Float64 nominalAmplitude);
  void setAsymProfileAndParams(const std::string& profile, TAsymParams params);
  void convertLineProfiles2SYMIGM(
      const std::shared_ptr<CSpectrumFluxCorrectionMeiksin> &igmcorrection);
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

class CLineCatalogsTplRatio
{

public:
  CLineCatalogsTplRatio();
  void addLineRatioCatalog(const CLineRatioCatalog &lr_catalog);

};


class CTemplateCatalog
{
public:
    CTemplateCatalog();
    void Add( std::shared_ptr<CTemplate> r);
};


class COperatorResult
{

public:
  COperatorResult(const std::string &type) : m_type(type){};
  virtual ~COperatorResult();

  const std::string &getType();

};

%template(TMapFloat64) std::map<std::string, Float64>;
%template(TZGridListParams) std::vector<CZGridParam>;

%include "method/classificationresult.i"
%include "method/reliabilityresult.i"
%include "method/linemodelsolveresult.i"
%include "operator/flagResult.i"
%include "statistics/pdfcandidatesz.i"
%include "operator/logZPdfResult.i"
%include "common/zgridparam.i"
%include "operator/extremaresult.i"
%include "operator/tplCombinationExtremaResult.i"
%include "linemodel/linemodelextremaresult.i"
%include "operator/modelspectrumresult.i"
%include "operator/modelphotvalueresult.i"
%include "linemodel/linemodelsolution.i"
%include "linemodel/continuummodelsolution.i"
%include "common/errorcodes.i"
%include "common/warningcodes.i"

class CProcessFlowContext {
public:

  static CProcessFlowContext& GetInstance();

  void Init();
  void setLineCatalog(const std::string& objectType, const std::string& method, const std::shared_ptr<CLineCatalog> &catalog); 
  void setLineRatioCatalogCatalog(const std::string& objectType, const std::shared_ptr<CLineCatalogsTplRatio> &catalog); 
  void setTemplateCatalog(const std::shared_ptr<CTemplateCatalog> &templateCatalog);
  void setPhotBandCatalog(const std::shared_ptr<CPhotBandCatalog> &photBandCatalog);
  void addSpectrum(const std::shared_ptr<CSpectrum> &spectrum);
  void setFluxCorrectionMeiksin(const std::shared_ptr<CSpectrumFluxCorrectionMeiksin> &igmcorrectionMeiksin);
  void setFluxCorrectionCalzetti(const std::shared_ptr<CSpectrumFluxCorrectionCalzetti> &ismcorrectionCalzetti);
  void reset();

  const std::shared_ptr<COperatorResultStore> &GetResultStore();
  std::shared_ptr<const CParameterStore> LoadParameterStore(const std::string& paramsJSONString);
 
  const std::string &GetCurrentCategory() const;
  const std::string &GetCurrentStage() const;
  const std::string &GetCurrentMethod() const;

  std::shared_ptr<CScopeStack>  m_ScopeStack;

 private:
    CProcessFlowContext();
    ~CProcessFlowContext();
    
};

enum class ScopeType { UNDEFINED, SPECTRUMMODEL, STAGE, METHOD };

%pythonappend CScopeStack::get_current_type() const %{
    val = ScopeType(val)
%}

class CScopeStack {
public:
  CScopeStack();
  CScopeStack(const TScopeStack &scope);
  CScopeStack(TScopeStack &&scope);
  CScopeStack(std::initializer_list<std::string> init);
  std::size_t size() const;
  bool empty() const;
  TStringList::const_iterator begin() const;
  TStringList::const_iterator end() const;
  std::string const &front() const;
  std::string const &back() const;
  std::string const &at(std::size_t pos) const;

  void push_back(const std::string &value, ScopeType type = ScopeType::UNDEFINED);
  void push_back(std::string &&value, ScopeType type = ScopeType::UNDEFINED);
  void pop_back();

  bool has_type(ScopeType type) const;
  const std::string &get_type_value(ScopeType type) const;
  size_t get_type_level(ScopeType type) const;
  ScopeType get_current_type() const;
};

%extend CScopeStack {
  std::string __getitem__(size_t pos) {
    return (*($self))[pos];
  }
  // std::string __setitem__(size_t pos, const std::string &value) {
  //     (*($self))[pos] = value;
  //     return (*($self))[pos];
  //   }
}




class CParameterStore : public CScopeStore
{
  public:
    CParameterStore(const std::shared_ptr<const CScopeStack> &stack);
    void FromString(const std::string& json);
    template<typename T> T Get(const std::string& name) const;
};

%template(Get_string) CParameterStore::Get<std::string>;

class COperatorResultStore
{

 public:
  COperatorResultStore(const std::shared_ptr<const CScopeStack> &scopeStack);

  std::shared_ptr<const CClassificationResult> GetClassificationResult(const std::string& objectType,
                                                                       const std::string& stage,
                                                                       const std::string& method,
                                                                       const std::string& name ) const;

  std::shared_ptr<const CReliabilityResult> GetReliabilityResult(const std::string& objectType, const std::string& stage,
									 const std::string& method,
                                                                       const std::string& name ) const;

  std::shared_ptr<const CLogZPdfResult> GetLogZPdfResult(const std::string& objectType, const std::string& stage,
								    const std::string& method,
								    const std::string& name ) const;

  std::shared_ptr<const CFlagLogResult> GetFlagLogResult(const std::string& objectType, const std::string& stage,
							 const std::string& method,
							 const std::string& name ) const;

  std::shared_ptr<const TLineModelResult> GetLineModelResult(const std::string& objectType, const std::string& stage,
							     const std::string& method,
							     const std::string& name ,
							     const std::string &dataset,
							     const int& rank,
                   bool firstpassResults
							     ) const;
  std::shared_ptr<const TTplCombinationResult> GetTplCombinationResult(const std::string& objectType, const std::string& stage,
										 const std::string& method,
										 const std::string& name ,
								     const std::string &dataset,
										 const int& rank
										 ) const;
  std::shared_ptr<const TExtremaResult> GetExtremaResult(const std::string& objectType, const std::string& stage,
										 const std::string& method,
										 const std::string& name ,
                     const std::string &dataset,
										 const int& rank,
                     bool firstpassResults
									       ) const;


  std::shared_ptr<const CLineModelSolution> GetLineModelSolution(const std::string& objectType, const std::string& stage,
								 const std::string& method,
								 const std::string& name,
								 const std::string &dataset,
								 const int& rank 
								 ) const  ;

    std::shared_ptr<const CLineModelSolution> GetLineModelSolution(const std::string& objectType, const std::string& stage,
								 const std::string& method,
								 const std::string& name 
								 ) const  ;

  std::shared_ptr<const CModelSpectrumResult> GetModelSpectrumResult(const std::string& objectType, const std::string& stage,
								     const std::string& method,
								     const std::string& name ,
								     const std::string &dataset,
								     const int& rank
								     ) const  ;

  std::shared_ptr<const CModelSpectrumResult> GetModelSpectrumResult(const std::string& objectType, const std::string& stage,
								     const std::string& method,
								     const std::string& name 
								     ) const  ;
  std::shared_ptr<const CModelPhotValueResult> GetModelPhotValueResult(const std::string &objectType, const std::string& stage,
                                            const std::string &method,
                                            const std::string &name,
                                            const std::string &dataset,
                                            const int &rank) const;

  std::shared_ptr<const CLineModelSolveResult>
  GetLineModelSolveResult(const std::string &spectrumModel, const std::string &stage,
			  const std::string &method, const std::string &name) const;

  const std::string&  GetGlobalResultType(const std::string& objectType, const std::string& stage,
                                          const std::string& method,
                                          const std::string& name ) const;

    const std::string&  GetCandidateResultType(const std::string& objectType, const std::string& stage,
					     const std::string& method,
					     const std::string& name ,
					     const std::string& dataset) const;

      bool HasCandidateDataset(const std::string& objectType, const std::string& stage,
			   const std::string& method,
			   const std::string& name ,
			   const std::string& dataset) const;

      bool HasDataset(const std::string& objectType, const std::string& stage,
			   const std::string& method,
			   const std::string& name ) const;

  int getNbRedshiftCandidates(const std::string& objectType, const std::string& stage,
			      const std::string& method) const;

  void StoreScopedFlagResult( const std::string& name);

  void StoreScopedGlobalResult(const std::string &name,
                               std::shared_ptr<const COperatorResult> result);

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

  void setObsID(const std::string& obsID);

  void ValidateNoise( Float64 LambdaMin,  Float64 LambdaMax ) const;
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
  TAxisSampleList& GetSamplesVector();
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
  virtual Float64 GetWidth(Float64 lambda, bool cliplambda = false) const=0;
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
    TLSFArguments(const std::shared_ptr<const CParameterStore>& parameterStore);
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

%pythonappend AmzException::getErrorCode() const %{
    val = ErrorCode(val)
%}

%pythonprepend AmzException::LogError() const %{
    if self.logged:
      return
    if not args:
      args = (self.__str__(), )

    self.logged = True
%}

class AmzException : public std::exception
{

 public: 

  AmzException(ErrorCode ec, const std::string &message, const char * filename_,
		  const char * method_, int line_);

  virtual ~AmzException();
 
  ErrorCode getErrorCode() const;
  virtual const char* what() const noexcept override;
  const std::string &getMessage() const;

  const std::string &getFileName() const;
  void setFilename(std::string const &filename_);
  const std::string &getMethod() const;
  void setMethod(std::string const &method_);
  int getLine() const;
  void setLine(int line_);

  void LogError(const std::string &msg = std::string()) const;

};

%pythoncode %{
  def _AmzException__str__(self):
    msg = f"{self.getErrorCode().name}: {self.getMessage()}"
    filename = self.getFileName()
    if filename:
      msg += f" [{filename}:{self.getLine()}:{self.getMethod()}]"
    return msg

  AmzException.__str__ = _AmzException__str__
  AmzException.logged = False
%}

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

    CClassificationSolve();

  };
  class CReliabilitySolve:public CSolve
  {

  public:

    CReliabilitySolve();
  };
  class CLineModelSolve:public CObjectSolve
  {

  public:

    CLineModelSolve();
  };
  class CLineMeasSolve:public CObjectSolve
  {

  public:

    CLineMeasSolve();
  };

  class CTemplateFittingSolve : public CObjectSolve
{
  public:

    CTemplateFittingSolve();
  };

class CTplCombinationSolve : public CObjectSolve
{

 public:
  CTplCombinationSolve();
};

class CLineMatchingSolve: public CObjectSolve
{
public:

    CLineMatchingSolve();
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
    CSpectrumFluxCorrectionMeiksin(std::vector<MeiksinCorrection> meiksinCorrectionCurves, TFloat64List zbins);
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

//code that runs after the cpp mapping takes place, it transfroms the cpp enum into python enum

%pythoncode %{
from enum import Enum, Flag
def redo(prefix, flag=False):
    tmpD = {k:v for k,v in globals().items() if k.startswith(prefix + '_')}
    for k,v in tmpD.items():
        del globals()[k]
    if flag:    
      tmpD = {k[len(prefix)+1:]:1<<v for k,v in tmpD.items()}
      flag = Flag(prefix,tmpD)
      for f in flag:
        f.bit = f.value.bit_length() - 1
      # python 3.9:  
      flag.bit_list = property(lambda self: [f.bit for f in flag if self.__contains__(f)])
      # python 3.11:
      # flag.bit_list = property(lambda self: [f.bit for f in self])
      # python 3.9:
      flag.name_list = property(lambda self: [f.name for f in flag if self.__contains__(f)]) 
      # python 3.11:
      # flag.name_list = flag.name.split('|') # only since python 3.11 (name is defined for aliases)
      globals()[prefix] = flag
    else:
      tmpD = {k[len(prefix)+1:]:v for k,v in tmpD.items()}
      globals()[prefix] = Enum(prefix,tmpD)
redo('ErrorCode')
redo('WarningCode', flag=True)
redo('ScopeType')
del redo  # cleaning up the namespace
del Enum
del Flag
%}

