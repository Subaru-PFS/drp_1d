#ifndef _REDSHIFT_OPERATOR_LINEMODELMATCHINGSOLVE_
#define _REDSHIFT_OPERATOR_LINEMODELMATCHINGSOLVE_

#include <epic/core/common/managedobject.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/method/linemodelsolveresult.h>
#include <epic/redshift/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

class CLineModelSolve : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( CLineModelSolve )

public:

    CLineModelSolve();
    ~CLineModelSolve();

    const std::string GetDescription();


    const CLineModelSolveResult* Compute(CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CRayCatalog& restraycatalog,
                                           const TFloat64Range& lambdaRange, const TFloat64List& redshifts );

    Bool Solve(CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CRayCatalog& restraycatalog,
                                 const TFloat64Range& lambdaRange, const TFloat64List& redshifts);

private:


};


}

#endif
