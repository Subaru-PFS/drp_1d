#ifndef _REDSHIFT_METHOD_CHISQUARELOGSOLVE_
#define _REDSHIFT_METHOD_CHISQUARELOGSOLVE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/templatefittingresult.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/operator/templatefittinglog.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

/**
 * \ingroup Redshift
 */
class CMethodTemplateFittingLogSolve
{

 public:


    CMethodTemplateFittingLogSolve( std::string calibrationPath="" );
    ~CMethodTemplateFittingLogSolve();

    const std::string GetDescription();

    std::shared_ptr<CTemplateFittingSolveResult> Compute(CDataStore& resultStore,
                                                   const CSpectrum& spc,
                                                   const CTemplateCatalog& tplCatalog,
                                                   const TStringList& tplCategoryList,
                                                   const TFloat64Range& lambdaRange,
                                                   const TFloat64List& redshifts,
                                                   Float64 overlapThreshold,
                                                   std::vector<CMask> maskList,
                                                   const std::string outputPdfRelDir,
                                                   const Float64 radius,
                                                   std::string spcComponent="raw",
                                                   std::string opt_interp="lin",
                                                   std::string opt_extinction="no",
                                                   std::string opt_dustFit="no");

    Bool ExtractCandidateResults(CDataStore& store, std::vector<Float64> zcandidates_unordered_list);


private:

    Bool Solve(CDataStore& resultStore,
               const CSpectrum& spc,
               const CTemplate& tpl,
               const TFloat64Range& lambdaRange, 
               const TFloat64List& redshifts,
               Float64 overlapThreshold,
               std::vector<CMask> maskList,
               CTemplateFittingSolveResult::EType spctype=CTemplateFittingSolveResult::nType_raw,
               std::string opt_interp="lin",
               std::string opt_extinction="no",
               std::string opt_dustFitting="no");

    Int32 CombinePDF(CDataStore& store,
                     std::string scopeStr,
                     std::string opt_combine,
                     std::shared_ptr<NSEpic::CPdfMargZLogResult> postmargZResult);


    std::string m_opt_pdfcombination;
    std::string m_opt_saveintermediateresults;
    Bool m_opt_enableSaveIntermediateTemplateFittingResults=false;
    std::string m_opt_spclogrebin;

    COperatorTemplateFittingLog* m_templateFittingOperator;
    Float64 m_radius;

};


}

#endif // _REDSHIFT_METHOD_CHISQUARELOGSOLVE_
