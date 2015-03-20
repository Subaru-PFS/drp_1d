#ifndef _REDSHIFT_PROCESSFLOW_CONTEXT_
#define _REDSHIFT_PROCESSFLOW_CONTEXT_

#include <epic/core/common/ref.h>
#include <epic/core/common/managedobject.h>
#include <epic/redshift/common/redshifts.h>
#include <epic/redshift/spectrum/template/template.h>

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

class CProcessFlowContext : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( CProcessFlowContext )

public:

    enum ESearchCriterion
    {
        nSearchCriterion_Maximized,
        nSearchCriterion_Minimized
    };

    struct SCorrelationResult
    {
        CRedshifts      SelectedRedshifts;
        TFloat64List    SelectedMerits;

        CRedshifts      Redshifts;
        TFloat64List    CorrelationValues;
    };

    typedef std::map< std::string, SCorrelationResult >     TCorrelationResults;
    typedef std::vector< CTemplate::ECategory >             TTemplateCategoryList;

    struct SParam
    {
        SParam();
        TTemplateCategoryList   templateCategoryList;
        TFloat64Range           lambdaRange;
        TFloat64Range           redshiftRange;
        Float64                 redshiftStep;
        Float64                 overlapThreshold;
        Int32                   smoothWidth;
    };


    CProcessFlowContext();
    ~CProcessFlowContext();

    bool Init( const char* spectrumPath, const char* noisePath, const char* tempalteCatalogPath, const char* rayCatalogPath, const SParam& params  );
    bool Init( const char* spectrumPath, const char* noisePath, CTemplateCatalog& templateCatalog, CRayCatalog& rayCatalog, const SParam& params  );

    CSpectrum&                      GetSpectrum();
    CSpectrum&                      GetSpectrumWithoutContinuum();
    CTemplateCatalog&               GetTemplateCatalog();
    CRayCatalog&                    GetRayCatalog();
    const TFloat64Range&            GetLambdaRange() const;
    const TFloat64Range&            GetRedshiftRange() const;
    const TTemplateCategoryList&    GetTemplateCategoryList() const;
    Float64                         GetOverlapThreshold() const;
    Float64                         GetRedshiftStep() const;

    Bool                            AddResults( const CTemplate& tpl, const CRedshifts& selectedRedshifts, const TFloat64List& selectedMerits, const CRedshifts& allRedshifts, const TFloat64List& allMerits  );
    Bool                            GetBestCorrelationResult( Float64& redshift, Float64& merit, std::string& tplName, ESearchCriterion criterion ) const;

    Bool                            DumpCorrelationResultsToCSV( const char* outputDirName ) const;

private:

    CRef<CSpectrum>                 m_Spectrum;
    CRef<CSpectrum>                 m_SpectrumWithoutContinuum;
    CRef<CTemplateCatalog>          m_TemplateCatalog;
    CRef<CRayCatalog>               m_RayCatalog;
    TFloat64Range                   m_LambdaRanges;
    TFloat64Range                   m_RedshiftRange;
    Float64                         m_OverlapThreshold;
    Float64                         m_RedshiftStep;
    TTemplateCategoryList           m_TemplateCategoryList;
    std::string                     m_SpectrumName;

    TCorrelationResults             m_CorrelationResult;


};


}

#endif
