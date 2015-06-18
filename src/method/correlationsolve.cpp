#include <epic/redshift/method/correlationsolve.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/operator/resultstore.h>

using namespace NSEpic;
using namespace std;

IMPLEMENT_MANAGED_OBJECT( COperatorCorrelationSolve )

COperatorCorrelationSolve::COperatorCorrelationSolve()
{

}

COperatorCorrelationSolve::~COperatorCorrelationSolve()
{

}

const CCorrelationSolveResult* COperatorCorrelationSolve::Compute(  COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                                        const CTemplateCatalog& tplCatalog, const TTemplateCategoryList& tplCategoryList,
                                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep,
                                                        Float64 overlapThreshold )
{
    Bool storeResult = false;

    COperatorResultStore::CAutoScope resultScope( resultStore, "correlationsolve" );

    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        CTemplate::ECategory category = tplCategoryList[i];

        for( UInt32 j=0; j<tplCatalog.GetTemplateCount( category ); j++ )
        {
            const CTemplate& tpl = tplCatalog.GetTemplate( category, j );
            const CTemplate& tplWithoutCont = tplCatalog.GetTemplateWithoutContinuum( category, j );

            Solve( resultStore, spc, spcWithoutCont, tpl, tplWithoutCont, lambdaRange, redshiftsRange, redshiftStep, overlapThreshold );

            storeResult = true;
        }
    }


    if( storeResult )
    {
        CCorrelationSolveResult*  solveResult = new CCorrelationSolveResult();
        return solveResult;
    }

    return NULL;
}

Bool COperatorCorrelationSolve::Solve( COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                               const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, Float64 overlapThreshold )
{
    CSpectrum s = spc;
    s.GetSpectralAxis().ConvertToLogScale();

    CTemplate t = tpl;
    t.GetSpectralAxis().ConvertToLogScale();


    // Create redshift initial list by spanning redshift acdross the given range, with the given delta
    TFloat64List redshifts = redshiftsRange.SpreadOver( redshiftStep );
    DebugAssert( redshifts.size() > 0 );

    // Compute correlation factor at each of those redshifts
    COperatorCorrelation correlation;
    CRef<CCorrelationResult> correlationResult = (CCorrelationResult*) correlation.Compute( spcWithoutCont, tplWithoutCont, lambdaRange, redshifts, overlapThreshold );

    if( !correlationResult )
    {
        return false;
    }

    resultStore.StorePerTemplateResult( tpl, "correlation", *correlationResult );


    return true;
}
