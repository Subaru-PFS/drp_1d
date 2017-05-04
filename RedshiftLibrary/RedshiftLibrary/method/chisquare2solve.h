#ifndef _REDSHIFT_METHOD_CHISQUARE2SOLVE_
#define _REDSHIFT_METHOD_CHISQUARE2SOLVE_


#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/chisquare2solveresult.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/operator/chisquare2.h>


namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

/**
 * \ingroup Redshift
 */
class CMethodChisquare2Solve
{

 public:


    CMethodChisquare2Solve( std::string calibrationPath="" );
    ~CMethodChisquare2Solve();
    const std::string GetDescription();

    std::shared_ptr<const CChisquare2SolveResult> Compute(CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                        const CTemplateCatalog& tplCatalog, const TStringList& tplCategoryList,
                                        const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold, std::vector<CMask> maskList, std::string spcComponent="raw" , std::string opt_interp="lin", std::string opt_extinction="no", std::string opt_dustFit="no");



private:

    Bool Solve(CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                                   const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold , std::vector<CMask> maskList, Int32 spctype=CChisquare2SolveResult::nType_raw, std::string opt_interp="lin", std::string opt_extinction="no", std::string opt_dustFitting="no");
    COperatorChiSquare2* m_chiSquareOperator;
};


}

#endif // _REDSHIFT_METHOD_CHISQUARE2SOLVE_


