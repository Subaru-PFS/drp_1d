#include <RedshiftLibrary/method/blindsolve.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/operator/correlation.h>
#include <RedshiftLibrary/operator/chisquare.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/processflow/datastore.h>

using namespace NSEpic;
using namespace std;


CMethodBlindSolve::CMethodBlindSolve()
{

}

CMethodBlindSolve::~CMethodBlindSolve()
{

}

const std::string CMethodBlindSolve::GetDescription()
{
    std::string desc;

    desc = "Method blindsolve:\n";

    desc.append("\tparam: blindsolve.overlapThreshold = <float value>\n");
    desc.append("\tparam: blindsolve.correlationExtremumCount = <float value>\n");

    return desc;

}

std::shared_ptr<CBlindSolveResult> CMethodBlindSolve::Compute( CDataStore& resultStore,
                                                               const CSpectrum& spc,
                                                               const CTemplateCatalog& tplCatalog,
                                                               const TStringList& tplCategoryList,
                                                               const TFloat64Range& lambdaRange,
                                                               const TFloat64Range& redshiftsRange,
                                                               Float64 redshiftStep,
                                                               Int32 correlationExtremumCount,
                                                               Float64 overlapThreshold )
{
    Bool storeResult = false;

    CDataStore::CAutoScope resultScope( resultStore, "blindsolve" );

    if(correlationExtremumCount==-1){
        Float64 count=0.0;
        resultStore.GetScopedParam( "correlationExtremumCount", count, 5.0 );
        correlationExtremumCount = (Int32)count;
    }
    if(overlapThreshold==-1.0){
        resultStore.GetScopedParam( "overlapThreshold", overlapThreshold, 1.0 );
    }

    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        std::string category = tplCategoryList[i];

        for( UInt32 j=0; j<tplCatalog.GetTemplateCount( category ); j++ )
        {
            const CTemplate& tpl = tplCatalog.GetTemplate( category, j );

            BlindSolve( resultStore, spc, tpl, lambdaRange, redshiftsRange, redshiftStep, correlationExtremumCount, overlapThreshold );

            storeResult = true;
        }
    }


    if( storeResult )
    {
        return  std::shared_ptr<CBlindSolveResult>( new CBlindSolveResult() );
    }

    return NULL;
}

Bool CMethodBlindSolve::BlindSolve( CDataStore& resultStore,
                                    const CSpectrum& spc,
                                    const CTemplate& tpl,
                                    const TFloat64Range& lambdaRange,
                                    const TFloat64Range& redshiftsRange,
                                    Float64 redshiftStep,
                                    Int32 correlationExtremumCount,
                                    Float64 overlapThreshold )
{
    CSpectrumSpectralAxis spcSpectralAxis = spc.GetSpectralAxis();
    spcSpectralAxis.ConvertToLogScale();
    CSpectrumFluxAxis spcFluxAxis = spc.GetWithoutContinuumFluxAxis();

    CSpectrumSpectralAxis tplSpectralAxis = tpl.GetSpectralAxis();
    tplSpectralAxis.ConvertToLogScale();
    CSpectrumFluxAxis tplFluxAxis = tpl.GetWithoutContinuumFluxAxis();

    CSpectrum spcWithoutCont(spcSpectralAxis, spcFluxAxis);
    CTemplate tplWithoutCont(tpl.GetName(), tpl.GetCategory(), tplSpectralAxis, tplFluxAxis);

    // Create redshift initial list by spanning redshift across the given range, with the given delta
    TFloat64List redshifts = redshiftsRange.SpreadOver( redshiftStep ); //TODO: this should be done in processflow, not in the method itself.
    DebugAssert( redshifts.size() > 0 );

    // prepare the unused masks
    std::vector<CMask> maskList;

    // Compute correlation factor at each of those redshifts
    COperatorCorrelation correlation;
    std::shared_ptr<const CCorrelationResult> correlationResult = dynamic_pointer_cast<const CCorrelationResult> ( correlation.Compute( spcWithoutCont, tplWithoutCont, lambdaRange, redshifts, overlapThreshold, maskList ) );

    if( !correlationResult )
    {
        return false;
    }

    resultStore.StoreScopedPerTemplateResult( tpl, "correlation", correlationResult );
    

    Float64 radius = 0.005;
    TPointList extremumList;
    CExtremum extremum( redshiftsRange, correlationExtremumCount, 2*radius);
    extremum.Find( correlationResult->Redshifts, correlationResult->Correlation, extremumList);
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
    auto chisquareResult = dynamic_pointer_cast<const CChisquareResult>( meritChiSquare.Compute( spc, tpl, lambdaRange, extremumRedshifts, overlapThreshold, maskList ) );
    if( !chisquareResult )
    {
        return false;
    }

    // Store results
    resultStore.StoreScopedPerTemplateResult( tpl, "merit", chisquareResult );

    return true;
}
