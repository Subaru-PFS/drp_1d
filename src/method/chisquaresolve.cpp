
#include <epic/redshift/method/chisquaresolve.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/operator/resultstore.h>

using namespace NSEpic;
using namespace std;

IMPLEMENT_MANAGED_OBJECT( CMethodChisquareSolve )

CMethodChisquareSolve::CMethodChisquareSolve()
{

}

CMethodChisquareSolve::~CMethodChisquareSolve()
{

}


const CChisquareSolveResult* CMethodChisquareSolve::Compute(  COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                                        const CTemplateCatalog& tplCatalog, const TTemplateCategoryList& tplCategoryList,
                                                        const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold )
{
    Bool storeResult = false;

    COperatorResultStore::CAutoScope resultScope( resultStore, "chisquaresolve" );

    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        CTemplate::ECategory category = tplCategoryList[i];

        for( UInt32 j=0; j<tplCatalog.GetTemplateCount( category ); j++ )
        {
            const CTemplate& tpl = tplCatalog.GetTemplate( category, j );
            const CTemplate& tplWithoutCont = tplCatalog.GetTemplateWithoutContinuum( category, j );

            Solve( resultStore, spc, spcWithoutCont, tpl, tplWithoutCont, lambdaRange, redshifts, overlapThreshold );

            storeResult = true;
        }
    }


    if( storeResult )
    {
        CChisquareSolveResult*  ChisquareSolveResult = new CChisquareSolveResult();
        return ChisquareSolveResult;
    }

    return NULL;
}

Bool CMethodChisquareSolve::Solve( COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                               const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold )
{

    // Compute merit function
    COperatorChiSquare chiSquare;
    CRef<CChisquareResult>  chisquareResult = (CChisquareResult*)chiSquare.Compute( spc, tpl, lambdaRange, redshifts, overlapThreshold );
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
