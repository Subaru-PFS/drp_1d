#include <epic/redshift/operator/blindsolve.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/operator/resultstore.h>

using namespace NSEpic;
using namespace std;

IMPLEMENT_MANAGED_OBJECT( COperatorBlindSolve )

COperatorBlindSolve::COperatorBlindSolve()
{

}

COperatorBlindSolve::~COperatorBlindSolve()
{

}

const CBlindSolveResult* COperatorBlindSolve::Compute(  COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                                        const CTemplateCatalog& tplCatalog, const TTemplateCategoryList& tplCategoryList,
                                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep,
                                                        Int32 correlationExtremumCount, Float64 overlapThreshold )
{
    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        CTemplate::ECategory category = tplCategoryList[i];
        if( category == CTemplate::nCategory_Star )
        {
        }
        else
        {
            for( UInt32 j=0; j<tplCatalog.GetTemplateCount( category ); j++ )
            {
                const CTemplate& tpl = tplCatalog.GetTemplate( category, j );
                const CTemplate& tplWithoutCont = tplCatalog.GetTemplateWithoutContinuum( category, j );

                BlindSolve( resultStore, spc, spcWithoutCont, tpl, tplWithoutCont, lambdaRange, redshiftsRange, redshiftStep, correlationExtremumCount, overlapThreshold );
            }

            CRef<CBlindSolveResult>  blindSolveResult = new CBlindSolveResult();
            resultStore.StoreGlobalResult(  "blindsolve", *blindSolveResult );

            return blindSolveResult;
        }
    }

    return NULL;
}

Bool COperatorBlindSolve::BlindSolve( COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                               const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, Int32 correlationExtremumCount,
                               Float64 overlapThreshold )
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

    resultStore.StorePerTemplateResult( tpl, "blindsolve.correlation", *correlationResult );

    // Find redshifts extremum
    TPointList extremumList;
    CExtremum extremum( redshiftsRange, correlationExtremumCount);
    extremum.Find( correlationResult->Redshifts, correlationResult->Correlation, extremumList );

    if( extremumList.size() == 0 )
    {
        return false;
    }

    // Compute merit function
    TFloat64List extremumRedshifts( extremumList.size() );
    TFloat64List extremumCorrelation( extremumList.size() );
    for( Int32 i=0; i<extremumList.size(); i++ )
    {
        extremumRedshifts[i] = extremumList[i].X;
        extremumCorrelation[i] = extremumList[i].Y;
    }

    COperatorChiSquare meritChiSquare;
    CRef<CCorrelationResult> chisquareResult = (CCorrelationResult*)meritChiSquare.Compute( spc, tpl, lambdaRange, extremumRedshifts, overlapThreshold );
    if( !chisquareResult )
    {
        return false;
    }

    // Store results
    resultStore.StorePerTemplateResult( tpl, "blindsolve.merit", *chisquareResult );

    return true;
}
