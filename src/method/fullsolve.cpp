
#include <epic/redshift/method/fullsolve.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/processflow/resultstore.h>

using namespace NSEpic;
using namespace std;

IMPLEMENT_MANAGED_OBJECT( COperatorFullSolve )

COperatorFullSolve::COperatorFullSolve()
{

}

COperatorFullSolve::~COperatorFullSolve()
{

}

const CFullSolveResult* COperatorFullSolve::Compute(  COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                                        const CTemplateCatalog& tplCatalog, const TTemplateCategoryList& tplCategoryList,
                                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, Float64 overlapThreshold )
{
    Bool storeResult = false;

    COperatorResultStore::CAutoScope resultScope( resultStore, "fullsolve" );

    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        CTemplate::ECategory category = tplCategoryList[i];

        for( UInt32 j=0; j<tplCatalog.GetTemplateCount( category ); j++ )
        {
            const CTemplate& tpl = tplCatalog.GetTemplate( category, j );
            const CTemplate& tplWithoutCont = tplCatalog.GetTemplateWithoutContinuum( category, j );

            SolveBrute( resultStore, spc, spcWithoutCont, tpl, tplWithoutCont, lambdaRange, redshiftsRange, redshiftStep, overlapThreshold );

            storeResult = true;
        }
    }


    if( storeResult )
    {
        CFullSolveResult*  fullSolveResult = new CFullSolveResult();
        return fullSolveResult;
    }

    return NULL;
}

Bool COperatorFullSolve::SolveBrute( COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                               const TFloat64Range& lambdaRange, const TFloat64Range& redshiftRange, Float64 redshiftStep,
                               Float64 overlapThreshold )
{


    CSpectrum s = spc;
    s.GetSpectralAxis().ConvertToLogScale();

    //CTemplate t = tpl;
    //t.GetSpectralAxis().ConvertToLogScale();

    // Create redshift initial list by spanning redshift acdross the given range, with the given delta
    TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
    DebugAssert( redshifts.size() > 0 );


    // Compute correlation factor at each of those redshifts
    COperatorCorrelation correlation;
    CRef<CCorrelationResult> correlationResult = (CCorrelationResult*) correlation.Compute( spcWithoutCont, tplWithoutCont, lambdaRange, redshifts, overlapThreshold );

    if( !correlationResult )
    {
        //Log.LogInfo( "Failed to compute correlation results");
        return false;
    }else{
        // Store results
        resultStore.StorePerTemplateResult( tpl, "correlation", *correlationResult );
    }

    COperatorChiSquare meritChiSquare;
    CRef<CChisquareResult> chisquareResult = (CChisquareResult*)meritChiSquare.Compute( spc, tpl, lambdaRange, redshifts, overlapThreshold );
    if( !chisquareResult )
    {
        //Log.LogInfo( "Failed to compute chi square value");
        return false;
    }else{
        // Store results
        resultStore.StorePerTemplateResult( tpl, "chisquare", *chisquareResult );
    }


    return true;
}

//todo, adapt this in the fullsolve method object
/* // Deactivated: Todo: has to be moved in operator Fullsolve (same as blindsolve)
Bool CProcessFlow::FullSolve( CProcessFlowContext& ctx, const CTemplate& tpl, const CTemplate& tplWithoutCont )
{
    Log.LogInfo( "Process FullSolve deactivated, running FullSolveBrute instead... )");

    const CSpectrum& spc = ctx.GetSpectrum();
    const CSpectrum& spcWithoutCont = ctx.GetSpectrumWithoutContinuum();

    CSpectrum s = spc;
    s.GetSpectralAxis().ConvertToLogScale();

    CTemplate t = tpl;
    t.GetSpectralAxis().ConvertToLogScale();

    // Create redshift initial list by spanning redshift acdross the given range, with the given delta
    TFloat64List redshifts = ctx.GetParams().redshiftRange.SpreadOver( ctx.GetParams().redshiftStep );
    DebugAssert( redshifts.size() > 0 );

    // Compute correlation factor at each of those redshifts
    COperatorChicorr chicorr;
    CChisquareResult* chisquareResult = new CChisquareResult();
    CCorrelationResult* correlationResult = new CCorrelationResult();
    chicorr.Compute( spc, spcWithoutCont, tpl, tplWithoutCont, ctx.GetParams().lambdaRange, redshifts, ctx.GetParams().overlapThreshold, correlationResult, chisquareResult);


    // Store results
    ctx.StorePerTemplateResult( tpl, "blindsolve.correlation", *correlationResult );
    ctx.StorePerTemplateResult( tpl, "blindsolve.merit", *chisquareResult );

    CRef<CBlindSolveResult>  chisquareResults = new CBlindSolveResult();
    ctx.StorePerTemplateResult( tpl, "blindsolve", *chisquareResults );
    return true;

    return FullSolveBrute( ctx, tpl, tplWithoutCont );
}
*/
