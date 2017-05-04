
#ifndef _REDSHIFT_OPERATOR_FULLSOLVE_
#define _REDSHIFT_OPERATOR_FULLSOLVE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/fullsolveresult.h>
#include <RedshiftLibrary/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

/**
 * \ingroup Redshift
 */
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
