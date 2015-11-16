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

/**
 * \ingroup Redshift
 */
class COperatorLineMatchingSolve : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( COperatorLineMatchingSolve )

public:

    COperatorLineMatchingSolve();
    ~COperatorLineMatchingSolve();


    const std::string GetDescription();

    const CLineMatchingSolveResult* Compute(CDataStore& resultStore, const CSpectrum& spc,
                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, const CRayCatalog &restRayCatalog);


private:


};


}

#endif
