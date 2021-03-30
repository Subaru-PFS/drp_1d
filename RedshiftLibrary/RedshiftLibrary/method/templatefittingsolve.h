#ifndef _REDSHIFT_METHOD_TEMPLATEFITTINGSOLVE_
#define _REDSHIFT_METHOD_TEMPLATEFITTINGSOLVE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/solve.h>
#include <RedshiftLibrary/processflow/resultstore.h>
#include <RedshiftLibrary/processflow/inputcontext.h>
#include <RedshiftLibrary/method/templatefittingsolveresult.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/operator/templatefittingBase.h>
#include <RedshiftLibrary/operator/pdfz.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

/**
 * \ingroup Redshift
 */
  class CMethodTemplateFittingSolve : public CSolve
{

 public:

    enum EType
    {
             nType_raw = 1,
             nType_continuumOnly = 2,
             nType_noContinuum = 3,
             nType_all = 4,
    };


  CMethodTemplateFittingSolve(TScopeStack &scope,std::string objectType);



  std::shared_ptr<CSolveResult> compute(const CInputContext &inputContext,
                                        COperatorResultStore &resultStore,
                                        TScopeStack &scope);

  /*
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
  */
    

private:

    Bool Solve(COperatorResultStore& resultStore,
               const CSpectrum& spc,
               const CTemplate& tpl,
               const TFloat64Range& lambdaRange,
               const TFloat64List& redshifts,
               Float64 overlapThreshold,
               std::vector<CMask> maskList,
               EType spctype=nType_raw,
               std::string opt_interp="lin",
               std::string opt_extinction="no",
               std::string opt_dustFitting="no");

    ChisquareArray BuildChisquareArray(const COperatorResultStore& store, const std::string & scopeStr) const;

    std::shared_ptr<const CExtremaResult>  SaveExtremaResult(   const COperatorResultStore& store,
                                                                const std::string & scopeStr,
                                                                const TCandidateZbyRank & ranked_zCandidates,
                                                                const CSpectrum& spc,
                                                                const CTemplateCatalog& tplCatalog,
                                                                const TStringList& tplCategoryList,
                                                                const TFloat64Range& lambdaRange,
                                                                Float64 overlapThreshold,
                                                                std::string opt_interp,
                                                                std::string opt_extinction
                                                                );

    void StoreExtremaResults(COperatorResultStore &dataStore, std::shared_ptr<const CExtremaResult> & ExtremaResult) const ;
    
  std::shared_ptr<COperatorTemplateFittingBase> m_templateFittingOperator;

    std::string m_opt_pdfcombination;
    Float64 m_redshiftSeparation;
    Int64 m_opt_maxCandidate;
    std::string m_opt_saveintermediateresults;
    Bool m_opt_enableSaveIntermediateTemplateFittingResults=false;
  std::string m_opt_spclogrebin;

};


}

#endif
