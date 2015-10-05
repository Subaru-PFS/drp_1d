#include <epic/redshift/method/correlationsolve.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/processflow/datastore.h>

using namespace NSEpic;
using namespace std;

IMPLEMENT_MANAGED_OBJECT( COperatorCorrelationSolve )

COperatorCorrelationSolve::COperatorCorrelationSolve()
{

}

COperatorCorrelationSolve::~COperatorCorrelationSolve()
{

}

const CCorrelationSolveResult* COperatorCorrelationSolve::Compute(  CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                                        const CTemplateCatalog& tplCatalog, const TTemplateCategoryList& tplCategoryList,
                                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep,
                                                        Float64 overlapThreshold )
{
    Bool storeResult = false;

    CDataStore::CAutoScope resultScope( resultStore, "correlationsolve" );

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

Bool COperatorCorrelationSolve::Solve( CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
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
    CRef<CCorrelationResult> result = (CCorrelationResult*) correlation.Compute( spcWithoutCont, tplWithoutCont, lambdaRange, redshifts, overlapThreshold );

    if( !result )
    {
        return false;
    }

    // extrema
    Int32 extremumCount = 10;
    TPointList extremumList;
    TFloat64Range range( result->Redshifts[0], result->Redshifts[result->Redshifts.size()-1] );
    CExtremum extremum( range, extremumCount);
    extremum.Find( result->Redshifts, result->Correlation, extremumList );
    // Refine Extremum with a second maximum search around the z candidates:
    // This corresponds to the finer xcorrelation in EZ Pandora (in standard_DP fctn in SolveKernel.py)
    Float64 radius = 0.001;
    for( Int32 i=0; i<extremumList.size(); i++ )
    {
        Float64 x = extremumList[i].X;
        Float64 left_border = max(range.GetBegin(), x-radius);
        Float64 right_border=min(range.GetEnd(), x+radius);

        TPointList extremumListFine;
        TFloat64Range rangeFine = TFloat64Range( left_border, right_border );
        CExtremum extremumFine( rangeFine , 1);
        extremumFine.Find( result->Redshifts, result->Correlation, extremumListFine );
        extremumList[i] = extremumListFine[0];
    }
    // store extrema results
    result->Extrema.resize( extremumCount );
    for( Int32 i=0; i<extremumList.size(); i++ )
    {

        result->Extrema[i] = extremumList[i].X;
    }

    resultStore.StorePerTemplateResult( tpl, "correlation", *result );


    return true;
}
