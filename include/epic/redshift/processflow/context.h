#ifndef _REDSHIFT_PROCESSFLOW_CONTEXT_
#define _REDSHIFT_PROCESSFLOW_CONTEXT_

#include <epic/core/common/ref.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/operator/result.h>
#include <epic/redshift/operator/resultstore.h>

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
        nMethod_LineMatching,
        nMethod_FullSolve,
        nMethod_DecisionalTree7,
        nMethod_Count,
        nMethod_None = -1
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

    static std::string              GetMethodName( EMethod method );
    void                            SaveRedshift( const char* dir );

private:


    CRef<CSpectrum>                 m_Spectrum;
    CRef<CSpectrum>                 m_SpectrumWithoutContinuum;
    CRef<CTemplateCatalog>          m_TemplateCatalog;
    CRef<CRayCatalog>               m_RayCatalog;

    SParam                          m_Params;


};


}

#endif
