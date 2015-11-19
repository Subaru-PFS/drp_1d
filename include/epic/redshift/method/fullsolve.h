
#ifndef _REDSHIFT_OPERATOR_FULLSOLVE_
#define _REDSHIFT_OPERATOR_FULLSOLVE_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/method/fullsolveresult.h>
#include <epic/redshift/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

class COperatorFullSolve
{

public:

    COperatorFullSolve();
    ~COperatorFullSolve();

    std::shared_ptr<const CFullSolveResult> Compute( CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                        const CTemplateCatalog& tplCatalog, const TStringList& tplCategoryList,
                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, Float64 overlapThreshold  );


private:

    Bool SolveBrute( CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                                   const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, Float64 overlapThreshold );
};


}

#endif
