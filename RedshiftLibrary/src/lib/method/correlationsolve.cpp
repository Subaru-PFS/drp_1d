#include <RedshiftLibrary/method/correlationsolve.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/operator/correlation.h>
#include <RedshiftLibrary/operator/chisquare.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/processflow/datastore.h>

using namespace NSEpic;
using namespace std;


CMethodCorrelationSolve::CMethodCorrelationSolve()
{

}

CMethodCorrelationSolve::~CMethodCorrelationSolve()
{

}

std::shared_ptr<CCorrelationSolveResult>  CMethodCorrelationSolve::Compute(  CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                                        const CTemplateCatalog& tplCatalog, const TStringList& tplCategoryList,
                                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep,
                                                        Float64 overlapThreshold )
{
    Bool storeResult = false;

    CDataStore::CAutoScope resultScope( resultStore, "correlationsolve" );

    if(overlapThreshold==-1.0){
        resultStore.GetScopedParam( "overlapThreshold", overlapThreshold, 1.0 );
    }

    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        std::string category = tplCategoryList[i];

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
        return std::shared_ptr<CCorrelationSolveResult>( new CCorrelationSolveResult() );
    }

    return NULL;
}

Bool CMethodCorrelationSolve::Solve( CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                               const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, Float64 overlapThreshold )
{
    CSpectrum s = spc;
    s.GetSpectralAxis().ConvertToLogScale();

    CTemplate t = tpl;
    t.GetSpectralAxis().ConvertToLogScale();


    // Create redshift initial list by spanning redshift acdross the given range, with the given delta
    TFloat64List redshifts = redshiftsRange.SpreadOver( redshiftStep ); //TODO: this should be done in processflow, not in the method itself.
    DebugAssert( redshifts.size() > 0 );

    // prepare the unused masks
    std::vector<CMask> maskList;

    // Compute correlation factor at each of those redshifts
    COperatorCorrelation correlation;
    auto result = dynamic_pointer_cast<CCorrelationResult>( correlation.Compute( spcWithoutCont, tplWithoutCont, lambdaRange, redshifts, overlapThreshold, maskList ) );

    if( !result )
    {
        return false;
    }

    // extrema
    Int32 extremumCount = 10;
    Float64 radius = 0.005;
    TPointList extremumList;
    TFloat64Range range( result->Redshifts[0], result->Redshifts[result->Redshifts.size()-1] );
    CExtremum extremum( range, extremumCount, 2*radius);
    extremum.Find( result->Redshifts, result->Correlation, extremumList);
    
    // store extrema results
    result->Extrema.resize( extremumCount );
    for( Int32 i=0; i<extremumList.size(); i++ )
    {

        result->Extrema[i] = extremumList[i].X;
    }

    resultStore.StoreScopedPerTemplateResult( tpl, "correlation", result );


    return true;
}
