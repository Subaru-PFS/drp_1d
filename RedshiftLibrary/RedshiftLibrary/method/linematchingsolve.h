#ifndef _REDSHIFT_METHOD_LINEMATCHINGSOLVE_
#define _REDSHIFT_METHOD_LINEMATCHINGSOLVE_

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
class CMethodLineMatchingSolve
{

public:

    CMethodLineMatchingSolve();
    ~CMethodLineMatchingSolve();

    const std::string GetDescription();

    std::shared_ptr<CLineMatchingSolveResult> Compute(CDataStore& resultStore,
                                                      const CSpectrum& spc,
                                                      const TFloat64Range& lambdaRange,
                                                      const TFloat64Range& redshiftsRange,
                                                      Float64 redshiftStep,
                                                      const CRayCatalog &restRayCatalog);


private:


};


}

#endif
