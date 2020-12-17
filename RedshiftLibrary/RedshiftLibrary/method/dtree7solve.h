#ifndef _REDSHIFT_METHOD_DTREE7SOLVE_
#define _REDSHIFT_METHOD_DTREE7SOLVE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/dtree7solveresult.h>
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
class CMethodDTree7Solve
{


public:

    CMethodDTree7Solve();
    ~CMethodDTree7Solve();

    const std::string GetDescription();

    std::shared_ptr<CDTree7SolveResult> Compute(CDataStore& dataStore,
                                                const CSpectrum& spc,
                                                const CTemplateCatalog& tplCatalog,
                                                const TStringList& tplCategoryList,
                                                const CRayCatalog& restRayCatalog,
                                                const TFloat64Range& lambdaRange,
                                                const TFloat64Range& redshiftRange,
                                                Float64 redshiftStep,
                                                const Float64 radius,
                                                Int32 correlationExtremumCount=-1,
                                                Float64 overlapThreshold=-1.0);


private:

    // Peak Detection
    Float64 m_winsize;
    Float64 m_cut;
    Float64 m_strongcut;

    // Line Matching
    Int64 m_minMatchNum;
    Float64 m_tol;

    // Tree path
    Float64 m_dtreepathnum;
    Float64 m_radius;


    Bool SolveDecisionalTree7(CDataStore& dataStore,
                              const CSpectrum& spc,
                              const CTemplateCatalog& tplCatalog,
                              const TStringList& tplCategoryList,
                              const CRayCatalog& restRayCatalog,
                              const TFloat64Range& lambdaRange,
                              const TFloat64Range& redshiftRange,
                              Float64 redshiftStep,
                              Int32 correlationExtremumCount,
                              Float64 overlapThreshold);

    TStringList getFilteredTplCategory(const TStringList& tplCategoryListIn, const std::string& CategoryFilter);

};


}

#endif
