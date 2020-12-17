#ifndef _REDSHIFT_METHOD_DTREECSOLVE_
#define _REDSHIFT_METHOD_DTREECSOLVE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/dtreecsolveresult.h>
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
class CMethodDTreeCSolve
{

public:

    CMethodDTreeCSolve( std::string calibrationPath="" );
    ~CMethodDTreeCSolve();

    const std::string GetDescription();

    std::shared_ptr<CDTreeCSolveResult> Compute(CDataStore& resultStore,
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
               const TFloat64List& redshifts,
               std::string& scopeStr);

    Bool GetCombinedRedshift(CDataStore& store, std::string scopeStr);

    TFloat64List GetBestRedshiftChi2List(CDataStore& store, std::string scopeStr, Float64& minmerit, TFloat64List& zList);

    TFloat64List GetMargChi2List(CDataStore& store, std::string scopeStr, Float64& minmerit, TFloat64List& zList);

    TFloat64List GetChi2ListForGivenTemplateName(CDataStore& store, std::string scopeStr, TFloat64List givenRedshifts, std::vector<std::string> givenTplNames);

    std::string m_calibrationPath;
    Float64 m_radius;
};


}

#endif
