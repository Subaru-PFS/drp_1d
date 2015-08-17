#ifndef _REDSHIFT_METHOD_CHISQUARESOLVE_
#define _REDSHIFT_METHOD_CHISQUARESOLVE_


#include <epic/core/common/managedobject.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/method/chisquaresolveresult.h>
#include <epic/redshift/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class COperatorResultStore;

class CMethodChisquareSolve : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( CMethodChisquareSolve )

public:

    CMethodChisquareSolve();
    ~CMethodChisquareSolve();
    const CChisquareSolveResult *Compute(   COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                        const CTemplateCatalog& tplCatalog, const TTemplateCategoryList& tplCategoryList,
                                        const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold  );



private:

    Bool Solve( COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                                   const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold );
};


}

#endif // _REDSHIFT_METHOD_CHISQUARESOLVE_

