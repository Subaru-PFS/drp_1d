#ifndef _REDSHIFT_METHOD_DTREEBSOLVE_
#define _REDSHIFT_METHOD_DTREEBSOLVE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/dtreebsolveresult.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class COperatorResultStore;

/**
 * \ingroup Redshift
 */
class CMethodDTreeBSolve
{

public:

    CMethodDTreeBSolve( std::string calibrationPath="" );
    ~CMethodDTreeBSolve();

    const std::string GetDescription();

    std::shared_ptr<CDTreeBSolveResult> Compute(CDataStore& resultStore,
                                                const CSpectrum& spc,
                                                const CTemplateCatalog& tplCatalog,
                                                const TStringList& tplCategoryList,
                                                const CRayCatalog& restRayCatalog,
                                                const TFloat64Range& lambdaRange,
                                                const TFloat64List& redshifts,
                                                const Float64 radius);


private:

    Bool Solve(CDataStore& resultStore,
               const CSpectrum& spc,
               const CTemplateCatalog& tplCatalog,
               const TStringList& tplCategoryList,
               const CRayCatalog& restRayCatalog,
               const TFloat64Range& lambdaRange,
               const TFloat64List& redshifts);

    Bool GetCombinedRedshift(CDataStore& store);

    TFloat64List GetBestRedshiftChi2List(CDataStore& store, std::string scopeStr, Float64& minmerit, TFloat64List& zList);

    std::string m_calibrationPath;
    Float64 m_radius;
};


}

#endif
