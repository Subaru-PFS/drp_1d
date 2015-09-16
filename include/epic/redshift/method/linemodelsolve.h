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
class COperatorResultStore;

class CLineModelSolve : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( CLineModelSolve )

public:

    CLineModelSolve();
    ~CLineModelSolve();

    const CLineModelSolveResult* Compute(  COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CRayCatalog& restraycatalog,
                                           const TFloat64Range& lambdaRange, const TFloat64List& redshifts );

    Bool Solve( COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CRayCatalog& restraycatalog,
                                 const TFloat64Range& lambdaRange, const TFloat64List& redshifts );
private:

    Float64 m_winsize;

};


}

#endif
