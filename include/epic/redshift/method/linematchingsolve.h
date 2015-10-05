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
class CDataStore;

class COperatorLineMatchingSolve : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( COperatorLineMatchingSolve )

public:

    COperatorLineMatchingSolve();
    ~COperatorLineMatchingSolve();

    const CLineMatchingSolveResult* Compute(CDataStore& resultStore, const CSpectrum& spc,
                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, const CRayCatalog &restRayCatalog);


private:

    // Peak Detection
    Float64 m_winsize;
    Float64 m_minsize;
    Float64 m_maxsize;
    Float64 m_cut;
    Float64 m_strongcut;
    Float64 m_enlargeRate;

    // Line Matching
    Int32 m_minMatchNum;
    Float64 m_tol;

};


}

#endif
