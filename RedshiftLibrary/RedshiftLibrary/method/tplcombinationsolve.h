#ifndef _REDSHIFT_METHOD_TPLCOMBINATIONSOLVE_
#define _REDSHIFT_METHOD_TPLCOMBINATIONSOLVE_


#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/tplcombinationsolveresult.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/operator/tplcombination.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

/**
 * \ingroup Redshift
 */
class CMethodTplcombinationSolve
{

 public:


    CMethodTplcombinationSolve( std::string calibrationPath="" );
    ~CMethodTplcombinationSolve();
    const std::string GetDescription();

    std::shared_ptr<CTplCombinationSolveResult> Compute(CDataStore& resultStore,
                                                    const CSpectrum& spc,
                                                    const CSpectrum& spcWithoutCont,
                                                    const CTemplateCatalog& tplCatalog,
                                                    const TStringList& tplCategoryList,
                                                    const TFloat64Range& lambdaRange,
                                                    const TFloat64List& redshifts,
                                                    Float64 overlapThreshold,
                                                    std::vector<CMask> maskList,
                                                    const std::string outputPdfRelDir,
                                                    std::string spcComponent="raw" ,
                                                    std::string opt_interp="lin",
                                                    std::string opt_extinction="no",
                                                    std::string opt_dustFit="no");



private:

    Bool Solve(CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                                   const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold , std::vector<CMask> maskList, Int32 spctype=CChisquare2SolveResult::nType_raw, std::string opt_interp="lin", std::string opt_extinction="no", std::string opt_dustFitting="no");
    Int32 CombinePDF(CDataStore& store,
                     std::string scopeStr,
                     std::string opt_combine,
                     std::shared_ptr<CPdfMargZLogResult> postmargZResult);


    COperatorTplCombination* m_tplcombinationOperator;


    std::string m_opt_pdfcombination;
    std::string m_opt_saveintermediateresults;
    Bool m_opt_enableSaveIntermediateChisquareResults=false;

};


}

#endif


