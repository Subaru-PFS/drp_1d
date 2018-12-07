#ifndef _REDSHIFT_OPERATOR_BLINDSOLVE_
#define _REDSHIFT_OPERATOR_BLINDSOLVE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/blindsolveresult.h>
#include <RedshiftLibrary/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

/**
 * \ingroup Redshift
 */
class COperatorBlindSolve
{

public:

    COperatorBlindSolve();
    ~COperatorBlindSolve();

    const std::string GetDescription();

    std::shared_ptr<CBlindSolveResult> Compute(   CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                        const CTemplateCatalog& tplCatalog, const TStringList& tplCategoryList,
                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep,
                                        Int32 correlationExtremumCount=-1, Float64 overlapThreshold=-1.0  );


private:

    Bool BlindSolve( CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                                   const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, Int32 correlationExtremumCount,
                                   Float64 overlapThreshold );
};


}

#endif
