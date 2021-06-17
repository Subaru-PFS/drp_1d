#ifndef _REDSHIFT_METHOD_TPLCOMBINATIONSOLVE_
#define _REDSHIFT_METHOD_TPLCOMBINATIONSOLVE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/templatefittingsolveresult.h>
#include <RedshiftLibrary/method/tplcombinationsolveresult.h>

#include <RedshiftLibrary/method/solve.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/operator/tplcombination.h>
#include <RedshiftLibrary/operator/pdfz.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CResultStore;

/**
 * \ingroup Redshift
 */
  class CMethodTplcombinationSolve : public CSolve
{

 public:

    enum EType
    {
             nType_raw = 1,
             nType_continuumOnly = 2,
             nType_noContinuum = 3,
             nType_all = 4,
    };

  CMethodTplcombinationSolve(TScopeStack &scope,std::string objectType);

    const std::string GetDescription() const;

    std::shared_ptr<CSolveResult> compute(std::shared_ptr<const CInputContext> inputContext,
                                        std::shared_ptr<COperatorResultStore> resultStore,
                                        TScopeStack &scope);

private:

  Bool Solve(std::shared_ptr<COperatorResultStore> resultStore,
               const CSpectrum& spc,
               const CTemplateCatalog& tplCatalog,
               const TStringList& tplCategoryList,
               const TFloat64Range& lambdaRange,
               const TFloat64List& redshifts,
               Float64 overlapThreshold,
               std::vector<CMask> maskList,
               EType spctype=nType_raw,
               std::string opt_interp="lin",
               std::string opt_extinction="no",
               std::string opt_dustFitting="no");
    
    ChisquareArray BuildChisquareArray(std::shared_ptr<COperatorResultStore> store, const std::string & scopeStr) const;
    void StoreExtremaResults( std::shared_ptr<COperatorResultStore> resultStore, 
                              std::shared_ptr<const ExtremaResult> & extremaResult) const;
    std::shared_ptr<const ExtremaResult> 
    SaveExtremaResult(std::shared_ptr<const COperatorResultStore> store,
                                               const std::string & scopeStr,
                                               const TCandidateZbyRank & ranked_zCandidates,
                                               const CSpectrum& spc,
                                               const CTemplateCatalog& tplCatalog,
                                               const TStringList& tplCategoryList,
                                               const TFloat64Range& lambdaRange,
                                               Float64 overlapThreshold,
                                               std::string opt_interp);
    COperatorTplcombination m_tplcombinationOperator;

    std::string m_opt_pdfcombination;
    Float64 m_redshiftSeparation;
    Int64 m_opt_maxCandidate;
    std::string m_opt_saveintermediateresults;
    Bool m_opt_enableSaveIntermediateChisquareResults=false;
};


}

#endif
