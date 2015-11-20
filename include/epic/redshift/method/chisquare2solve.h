#ifndef _REDSHIFT_METHOD_CHISQUARE2SOLVE_
#define _REDSHIFT_METHOD_CHISQUARE2SOLVE_


#include <epic/core/common/datatypes.h>
#include <epic/redshift/method/chisquare2solveresult.h>
#include <epic/redshift/spectrum/template/template.h>

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


    CMethodChisquare2Solve();
    ~CMethodChisquare2Solve();
    const std::string GetDescription();

    std::shared_ptr<const CChisquare2SolveResult> Compute(   CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                        const CTemplateCatalog& tplCatalog, const TStringList& tplCategoryList,
                                        const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold, std::string spcComponent="raw" );



private:

    Bool Solve(CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                                   const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold , Int32 spctype=CChisquare2SolveResult::nType_raw);
};


}

#endif // _REDSHIFT_METHOD_CHISQUARE2SOLVE_


