#ifndef _REDSHIFT_OPERATOR_DTREEASOLVE_
#define _REDSHIFT_OPERATOR_DTREEASOLVE_

#include <epic/core/common/managedobject.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/method/dtreeasolveresult.h>
#include <epic/redshift/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

class COperatorDTreeASolve : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( COperatorDTreeASolve )

public:

    COperatorDTreeASolve();
    ~COperatorDTreeASolve();

    const CDTreeASolveResult* Compute(CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                        const CTemplateCatalog& tplCatalog, const TStringList& tplCategoryList, const CRayCatalog &restRayCatalog,
                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftRange, Float64 redshiftStep,
                                        Int32 correlationExtremumCount, Float64 overlapThreshold  );


private:
    // Peak Detection
    Float64 m_winsize;
    Float64 m_minsize;
    Float64 m_maxsize;
    Float64 m_cut;
    Float64 m_strongcut;
    Float64 m_enlargeRate;

    // Line Matching
    Int32 m_minMatchNum;
    Float64 m_tol;


    Bool Solve(CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                              const CTemplateCatalog& tplCatalog, const TStringList& tplCategoryList, const CRayCatalog &restRayCatalog,
                              const TFloat64Range& lambdaRange, const TFloat64Range& redshiftRange, Float64 redshiftStep,
                              Int32 correlationExtremumCount, Float64 overlapThreshold );

    TStringList getFilteredTplCategory( const TStringList& tplCategoryListIn, const std::string& CategoryFilter);
};


}

#endif