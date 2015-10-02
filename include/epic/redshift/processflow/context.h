#ifndef _REDSHIFT_PROCESSFLOW_CONTEXT_
#define _REDSHIFT_PROCESSFLOW_CONTEXT_

#include <epic/core/common/ref.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/processflow/result.h>
#include <epic/redshift/processflow/resultstore.h>
#include <epic/redshift/processflow/parameterstore.h>

#include <map>
#include <string>


namespace NSEpic
{

class CSpectrum;
class CContinuum;
class CPeakStore;
class CRayCatalog;
class CTemplateCatalog;
class CRayCatalog;

/**
 * Store all data concerning computation and processign of a given spectrum
 *
 */
class CProcessFlowContext : public COperatorResultStore
{

    DEFINE_MANAGED_OBJECT( CProcessFlowContext )

public:

    enum EMethod
    {
         nMethod_BlindSolve = 0,
         nMethod_Correlation,
         nMethod_Chisquare,
         nMethod_LineMatching,
         nMethod_LineMatching2,
         nMethod_LineModel,
         nMethod_FullSolve,
         nMethod_DecisionalTree7,
         nMethod_DecisionalTreeA,
         nMethod_Count,
         nMethod_None = -1
    };

    struct SParam
    {
        SParam();
        TTemplateCategoryList   templateCategoryList;
        TFloat64Range           lambdaRange;
        TFloat64Range           redshiftRange;
        Float64                 redshiftStep;
        Float64                 overlapThreshold;
        Int32                   smoothWidth;
        EMethod                 method;
        Int32                   correlationExtremumCount;
    };


    CProcessFlowContext();
    ~CProcessFlowContext();

    bool Init( const char* spectrumPath, const char* noisePath, const char* tempalteCatalogPath, const char* rayCatalogPath, const SParam& params  );
    bool Init( const char* spectrumPath, const char* noisePath, const CTemplateCatalog& templateCatalog, const CRayCatalog& rayCatalog, const SParam& params  );

    const CSpectrum&                GetSpectrum() const;
    const CSpectrum&                GetSpectrumWithoutContinuum() const;
    const CTemplateCatalog&         GetTemplateCatalog() const;
    const CRayCatalog&              GetRayCatalog() const;
    const SParam&                   GetParams() const;

    static std::string              GetMethodName( EMethod method );

    Bool                            GetScopedParam( const char* name, TFloat64List& v, const TFloat64List& defaultValue = TFloat64List() );
    Bool                            GetScopedParam( const char* name, TInt64List& v, const TInt64List& defaultValue = TInt64List() );
    Bool                            GetScopedParam( const char* name, TBoolList& v, const TBoolList& defaultValue = TBoolList() );
    Bool                            GetScopedParam( const char* name, Float64& v, Float64 defaultValue  = 0 );
    Bool                            GetScopedParam( const char* name, Int64& v, Int64 defaultValue = 0 );
    Bool                            GetScopedParam( const char* name, Bool& v, Bool defaultValue = true );

private:

    CMethodParameterStore           m_ParameterStore;

    CRef<CSpectrum>                 m_Spectrum;
    CRef<CSpectrum>                 m_SpectrumWithoutContinuum;

    CRef<CTemplateCatalog>          m_TemplateCatalog;
    CRef<CRayCatalog>               m_RayCatalog;

    SParam                          m_Params;


};


}

#endif
