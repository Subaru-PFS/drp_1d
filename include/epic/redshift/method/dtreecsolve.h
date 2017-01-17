#ifndef _REDSHIFT_OPERATOR_DTREECSOLVE_
#define _REDSHIFT_OPERATOR_DTREECSOLVE_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/method/dtreecsolveresult.h>
#include <epic/redshift/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class COperatorResultStore;

/**
 * \ingroup Redshift
 */
class COperatorDTreeCSolve
{

public:

    COperatorDTreeCSolve(std::string calibrationPath);
    ~COperatorDTreeCSolve();

    const std::string GetDescription();

    std::shared_ptr<const CDTreeCSolveResult> Compute(CDataStore &resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                        const CTemplateCatalog& tplCatalog, const TStringList& tplCategoryList, const CRayCatalog &restRayCatalog,
                                        const TFloat64Range& lambdaRange, const TFloat64List& redshifts );


private:


    Bool Solve(CDataStore &resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                              const CTemplateCatalog& tplCatalog, const TStringList& tplCategoryList, const CRayCatalog &restRayCatalog,
                              const TFloat64Range& lambdaRange, const TFloat64List& redshifts );

    Bool GetCombinedRedshift(CDataStore& store, std::string scopeStr);
    TFloat64List GetBestRedshiftChi2List(CDataStore& store, std::string scopeStr, Float64 &minmerit, TFloat64List &zList);
    TFloat64List GetChi2ListForGivenTemplateName(CDataStore& store, std::string scopeStr, TFloat64List givenRedshifts, std::vector<std::string> givenTplNames);

    std::string m_calibrationPath;
};


}

#endif
