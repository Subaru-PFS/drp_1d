#ifndef _REDSHIFT_OPERATOR_CORRELATIONSOLVE_
#define _REDSHIFT_OPERATOR_CORRELATIONSOLVE_

#include <epic/core/common/managedobject.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/method/correlationsolveresult.h>
#include <epic/redshift/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class COperatorResultStore;

class COperatorCorrelationSolve : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( COperatorCorrelationSolve )

public:

    COperatorCorrelationSolve();
    ~COperatorCorrelationSolve();

    const CCorrelationSolveResult* Compute(   COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                        const CTemplateCatalog& tplCatalog, const TTemplateCategoryList& tplCategoryList,
                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep,
                                        Float64 overlapThreshold  );


private:

    Bool Solve( COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                                   const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, Float64 overlapThreshold );
};


}

#endif
