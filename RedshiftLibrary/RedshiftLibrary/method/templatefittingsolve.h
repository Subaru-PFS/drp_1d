#ifndef _REDSHIFT_METHOD_CHISQUARE2SOLVE_
#define _REDSHIFT_METHOD_CHISQUARE2SOLVE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/templatefittingresult.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/operator/templatefitting.h>
#include <RedshiftLibrary/operator/pdfz.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>
#include <RedshiftLibrary/operator/modelcontinuumfittingresult.h>
#include <RedshiftLibrary/operator/modelspectrumresult.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

/**
 * \ingroup Redshift
 */
class CMethodTemplateFittingSolve
{

 public:

    CMethodTemplateFittingSolve() = default;
    ~CMethodTemplateFittingSolve() = default;

    const std::string GetDescription() const;

    std::shared_ptr<CTemplateFittingSolveResult> Compute(CDataStore& resultStore,
                                                   const CSpectrum& spc,
                                                   const CTemplateCatalog& tplCatalog,
                                                   const TStringList& tplCategoryList,
                                                   const TFloat64Range& lambdaRange,
                                                   const TFloat64List& redshifts,
                                                   Float64 overlapThreshold,
                                                   std::vector<CMask> maskList,
                                                   const std::string outputPdfRelDir,
                                                   const Float64 redshiftSeparation,
                                                   std::string spcComponent="raw" ,
                                                   std::string opt_interp="lin",
                                                   std::string opt_extinction="no",
                                                   std::string opt_dustFit="no");

    

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

    ChisquareArray BuildChisquareArray(const CDataStore& store, const std::string & scopeStr) const;
    void SaveSpectrumResults(CDataStore &dataStore) const ;
    
    COperatorTemplateFitting m_templateFittingOperator;
    std::vector<std::shared_ptr<CModelSpectrumResult>>  m_savedModelSpectrumResults;
    std::vector<std::shared_ptr<CModelContinuumFittingResult>> m_savedModelContinuumFittingResults;

    std::string m_opt_pdfcombination;
    Float64 m_redshiftSeparation;
    Int64 m_opt_maxCandidate;
    std::string m_opt_saveintermediateresults;
    Bool m_opt_enableSaveIntermediateTemplateFittingResults=false;
};


}

#endif // _REDSHIFT_METHOD_CHISQUARE2SOLVE_
