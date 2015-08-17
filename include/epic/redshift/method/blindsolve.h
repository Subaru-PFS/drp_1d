#ifndef _REDSHIFT_OPERATOR_BLINDSOLVE_
#define _REDSHIFT_OPERATOR_BLINDSOLVE_

#include <epic/core/common/managedobject.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/method/blindsolveresult.h>
#include <epic/redshift/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class COperatorResultStore;

class COperatorBlindSolve : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( COperatorBlindSolve )

public:

    COperatorBlindSolve();
    ~COperatorBlindSolve();

    const CBlindSolveResult* Compute(   COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                        const CTemplateCatalog& tplCatalog, const TTemplateCategoryList& tplCategoryList,
                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep,
                                        Int32 correlationExtremumCount, Float64 overlapThreshold  );


private:

    Bool BlindSolve( COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                                   const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, Int32 correlationExtremumCount,
                                   Float64 overlapThreshold );
};


}

#endif
