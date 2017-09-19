#ifndef _REDSHIFT_OPERATOR_LINEMATCHINGSOLVE_
#define _REDSHIFT_OPERATOR_LINEMATCHINGSOLVE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/linematchingsolveresult.h>
#include <RedshiftLibrary/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

/**
 * \ingroup Redshift
 */
class COperatorLineMatchingSolve
{

public:

    COperatorLineMatchingSolve();
    ~COperatorLineMatchingSolve();


    const std::string GetDescription();

    std::shared_ptr<CLineMatchingSolveResult> Compute(CDataStore& resultStore, const CSpectrum& spc,
                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, const CRayCatalog &restRayCatalog);


private:


};


}

#endif
