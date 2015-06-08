#ifndef _REDSHIFT_PROCESSFLOW_CONTEXT_
#define _REDSHIFT_PROCESSFLOW_CONTEXT_

#include <epic/core/common/ref.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/operator/result.h>

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
class CProcessFlowContext : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( CProcessFlowContext )

public:

    typedef std::vector< CTemplate::ECategory >                 TTemplateCategoryList;
    typedef std::map< std::string, CConstRef<COperatorResult> > TResultsMap;
    typedef std::map< std::string, TResultsMap>                 TPerTemplateResultsMap;

    enum EMethod
    {
        nMethod_BlindSolve = 1,
        nMethod_LineMatching = 2,
        nMethod_FullSolve = 3,
        nMethod_DecisionalTree7 = 4,
    };

    enum EDTREEPATH
    {
        nDtreePath_None = -1,
        nDtreePath_BlindSolve = 1,
        nDtreePath_OnlyFit = 2,
        nDtreePath_FullSolve = 3,
        nDtreePath_OnlyCorrelation = 4,
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
    Int32                           m_dtreepath;
    Float64                         m_dtreepathnum;

    Void  StorePerTemplateResult( const CTemplate& t, const char* name, const COperatorResult& result );
    Void  StoreGlobalResult( const char* name, const COperatorResult& result );

    const COperatorResult*  GetPerTemplateResult( const CTemplate& t, const char* name ) const;
    TOperatorResultMap      GetPerTemplateResult( const char* name ) const;
    const COperatorResult*  GetGlobalResult( const char* name ) const;

private:

    void StoreResult( TResultsMap& map, const char* name, const COperatorResult& result );

    CRef<CSpectrum>                 m_Spectrum;
    CRef<CSpectrum>                 m_SpectrumWithoutContinuum;
    CRef<CTemplateCatalog>          m_TemplateCatalog;
    CRef<CRayCatalog>               m_RayCatalog;


    SParam                          m_Params;
    std::string                     m_SpectrumName;


    TPerTemplateResultsMap          m_PerTemplateResults;
    TResultsMap                     m_GlobalResults;

};


}

#endif
