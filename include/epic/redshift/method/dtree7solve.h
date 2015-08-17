#ifndef _REDSHIFT_OPERATOR_DTREE7SOLVE_
#define _REDSHIFT_OPERATOR_DTREE7SOLVE_

#include <epic/core/common/managedobject.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/method/dtree7solveresult.h>
#include <epic/redshift/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class COperatorResultStore;

class COperatorDTree7Solve : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( COperatorDTree7solve )

public:

    COperatorDTree7Solve();
    ~COperatorDTree7Solve();

    const CDTree7SolveResult* Compute(COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                        const CTemplateCatalog& tplCatalog, const TTemplateCategoryList& tplCategoryList, const CRayCatalog &restRayCatalog,
                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftRange, Float64 redshiftStep,
                                        Int32 correlationExtremumCount, Float64 overlapThreshold  );


private:
    // Peak Detection
    Float64 m_winsize;
    Float64 m_cut;
    Float64 m_strongcut;

    // Line Matching
    Int32 m_minMatchNum;
    Float64 m_tol;

    // tree path
    Float64 m_dtreepathnum;


    Bool SolveDecisionalTree7(COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                              const CTemplateCatalog& tplCatalog, const TTemplateCategoryList& tplCategoryList, const CRayCatalog &restRayCatalog,
                              const TFloat64Range& lambdaRange, const TFloat64Range& redshiftRange, Float64 redshiftStep,
                              Int32 correlationExtremumCount, Float64 overlapThreshold );

    TTemplateCategoryList getFilteredTplCategory(TTemplateCategoryList tplCategoryListIn, CTemplate::ECategory CategoryFilter);
};


}

#endif
