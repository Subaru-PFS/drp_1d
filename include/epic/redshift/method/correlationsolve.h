#ifndef _REDSHIFT_OPERATOR_CORRELATIONSOLVE_
#define _REDSHIFT_OPERATOR_CORRELATIONSOLVE_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/method/correlationsolveresult.h>
#include <epic/redshift/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

/**
 * \ingroup Redshift
 */
class COperatorCorrelationSolve
{

public:

    COperatorCorrelationSolve();
    ~COperatorCorrelationSolve();
    const std::string GetDescription();

    std::shared_ptr<const CCorrelationSolveResult>  Compute( CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                        const CTemplateCatalog& tplCatalog, const TStringList& tplCategoryList,
                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep,
                                        Float64 overlapThreshold=-1.0  );


private:

    Bool Solve( CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                                   const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, Float64 overlapThreshold );
};


}

#endif
