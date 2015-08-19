#ifndef _REDSHIFT_METHOD_CHISQUARE2SOLVE_
#define _REDSHIFT_METHOD_CHISQUARE2SOLVE_


#include <epic/core/common/managedobject.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/method/chisquare2solveresult.h>
#include <epic/redshift/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class COperatorResultStore;

class CMethodChisquare2Solve : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( CMethodChisquare2Solve )

    public:


    CMethodChisquare2Solve();
    ~CMethodChisquare2Solve();
    const CChisquare2SolveResult *Compute(   COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                        const CTemplateCatalog& tplCatalog, const TTemplateCategoryList& tplCategoryList,
                                        const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold  );



private:

    Bool Solve(COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                                   const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold , Int32 spctype=CChisquare2SolveResult::nType_raw);
};


}

#endif // _REDSHIFT_METHOD_CHISQUARE2SOLVE_


