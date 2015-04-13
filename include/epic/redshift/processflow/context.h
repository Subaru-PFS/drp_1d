#ifndef _REDSHIFT_PROCESSFLOW_CONTEXT_
#define _REDSHIFT_PROCESSFLOW_CONTEXT_

#include <epic/core/common/ref.h>
#include <epic/core/common/managedobject.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/operator/operator.h>

#include <epic/redshift/ray/matching.h>

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

    struct SResults
    {
        TFloat64List                SelectedRedshifts;
        TFloat64List                SelectedMerits;
        TFloat64List                SelectedCorrelations;

        COperator::TStatusList      SelectedStatus;

        TFloat64List                AllRedshifts;
        TFloat64List                AllCorrelation;
    };

    typedef std::map< std::string, SResults >       TResultsMap;
    typedef std::vector< CTemplate::ECategory >     TTemplateCategoryList;

    struct SRayMatchingResult
    {
        Float64     BestRedshift;
        Int32       BestRedshiftMatchingNumber;
        TRedshiftSolutionSetList  MatchingSolutions;
    };
    typedef std::map< std::string, SRayMatchingResult >     TRayMatchingResults;



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
    bool Init( const char* spectrumPath, const char* noisePath, const CTemplateCatalog& templateCatalog, const CRayCatalog& rayCatalog, const SParam& params  );

    const CSpectrum&                GetSpectrum();
    const CSpectrum&                GetSpectrumWithoutContinuum();
    const CTemplateCatalog&         GetTemplateCatalog();
    const CRayCatalog&              GetRayCatalog();
    const TFloat64Range&            GetLambdaRange() const;
    const TFloat64Range&            GetRedshiftRange() const;
    const TTemplateCategoryList&    GetTemplateCategoryList() const;
    Float64                         GetOverlapThreshold() const;
    Float64                         GetRedshiftStep() const;

    Bool                            AddResults( const CTemplate& tpl,
                                                const TFloat64List& selectedRedshifts, const TFloat64List& selectedCorrelation,
                                                const TFloat64List& selectedMerits, const COperator::TStatusList& selectedStatus,
                                                const TFloat64List& allRedshifts, const TFloat64List& allCorrelation );

    const TResultsMap&              GetResults() const;
    Bool                            AddMeritResults( const CTemplate& tpl,
                                          const TFloat64List& selectedRedshifts,
                                          const TFloat64List& selectedMerits, const COperator::TStatusList& selectedMeritsStatus,
                                          const TFloat64List& redshifts);
    Bool                            SetRayDetectionResult(CRayCatalog& detectedRayCatalog);
    CRayCatalog&                    GetDetectedRayCatalog();
    Bool                            SetRayMatchingResult(const TRedshiftSolutionSetList &allresults, Float64 bestRedshift, Int32 bestRedshiftMatchingNumber);
    Bool                            GetBestRayMatchingResult(Float64& bestRedshift, Float64& bestRedshiftMatchingNumber) const;
    Bool                            GetBestCorrelationResult( Float64& redshift, Float64& merit, std::string& tplName ) const;

    Bool                            DumpCorrelationResultsToCSV( const char* outputDirName ) const;
    Bool                            GetIntermediateResults(std::string& corrStr, std::string& fitStr);

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

    TResultsMap                     m_Results;

    CRef<CRayCatalog>               m_DetectedRayCatalog;
    SRayMatchingResult              m_RayMatchingResult;

};


}

#endif
