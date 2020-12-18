#ifndef _REDSHIFT_METHOD_CORRELATIONSOLVE_
#define _REDSHIFT_METHOD_CORRELATIONSOLVE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/correlationsolveresult.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

/**
 * \ingroup Redshift
 */
class CMethodCorrelationSolve
{

public:

    CMethodCorrelationSolve();
    ~CMethodCorrelationSolve();

    const std::string GetDescription();

    std::shared_ptr<CCorrelationSolveResult> Compute( CDataStore& resultStore,
                                                      const CSpectrum& spc,
                                                      const CTemplateCatalog& tplCatalog,
                                                      const TStringList& tplCategoryList,
                                                      const TFloat64Range& lambdaRange,
                                                      const TFloat64Range& redshiftsRange,
                                                      Float64 redshiftStep,
                                                      Float64 overlapThreshold=-1.0 );


private:

    Bool Solve( CDataStore& resultStore,
                const CSpectrum& spc,
                const CTemplate& tpl,
                const TFloat64Range& lambdaRange,
                const TFloat64Range& redshiftsRange,
                Float64 redshiftStep,
                Float64 overlapThreshold );
};


}

#endif
