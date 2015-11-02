#ifndef _REDSHIFT_OPERATOR_DTREEBSOLVE_
#define _REDSHIFT_OPERATOR_DTREEBSOLVE_

#include <epic/core/common/managedobject.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/method/dtreebsolveresult.h>
#include <epic/redshift/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class COperatorResultStore;

class COperatorDTreeBSolve : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( COperatorDTreeBSolve )

public:

    COperatorDTreeBSolve();
    ~COperatorDTreeBSolve();

    const CDTreeBSolveResult* Compute(COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                        const CTemplateCatalog& tplCatalog, const TTemplateCategoryList& tplCategoryList, const CRayCatalog &restRayCatalog,
                                        const TFloat64Range& lambdaRange, const TFloat64List& redshifts );


private:


    Bool Solve(COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                              const CTemplateCatalog& tplCatalog, const TTemplateCategoryList& tplCategoryList, const CRayCatalog &restRayCatalog,
                              const TFloat64Range& lambdaRange, const TFloat64List& redshifts );

    TTemplateCategoryList getFilteredTplCategory(TTemplateCategoryList tplCategoryListIn, CTemplate::ECategory CategoryFilter);
};


}

#endif
