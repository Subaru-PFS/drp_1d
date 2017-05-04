#ifndef _REDSHIFT_OPERATOR_DTREEBSOLVE_
#define _REDSHIFT_OPERATOR_DTREEBSOLVE_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/method/dtreebsolveresult.h>
#include <epic/redshift/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class COperatorResultStore;

/**
 * \ingroup Redshift
 */
class COperatorDTreeBSolve
{

public:

    COperatorDTreeBSolve( std::string calibrationPath="" );
    ~COperatorDTreeBSolve();

    const std::string GetDescription();

    std::shared_ptr<const CDTreeBSolveResult> Compute(CDataStore &resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                        const CTemplateCatalog& tplCatalog, const TStringList& tplCategoryList, const CRayCatalog &restRayCatalog,
                                        const TFloat64Range& lambdaRange, const TFloat64List& redshifts );


private:

    Bool Solve(CDataStore &resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                              const CTemplateCatalog& tplCatalog, const TStringList& tplCategoryList, const CRayCatalog &restRayCatalog,
                              const TFloat64Range& lambdaRange, const TFloat64List& redshifts );

    Bool GetCombinedRedshift(CDataStore& store);
    TFloat64List GetBestRedshiftChi2List(CDataStore& store, std::string scopeStr, Float64 &minmerit, TFloat64List &zList);

    std::string m_calibrationPath;

};


}

#endif