#ifndef _REDSHIFT_OPERATOR_LINEMATCHINGSOLVE_
#define _REDSHIFT_OPERATOR_LINEMATCHINGSOLVE_

#include <epic/core/common/managedobject.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/method/linematchingsolveresult.h>
#include <epic/redshift/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class COperatorResultStore;

class COperatorLineMatchingSolve : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( COperatorLineMatchingSolve )

public:

    COperatorLineMatchingSolve();
    ~COperatorLineMatchingSolve();

    const CLineMatchingSolveResult* Compute(COperatorResultStore& resultStore, const CSpectrum& spc,
                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, const CRayCatalog &restRayCatalog);


private:

    // Peak Detection
    Float64 m_winsize;
    Float64 m_cut;
    Float64 m_strongcut;

    // Line Matching
    Int32 m_minMatchNum;
    Float64 m_tol;

};


}

#endif
