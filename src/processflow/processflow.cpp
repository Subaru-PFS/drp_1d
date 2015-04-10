#include <epic/redshift/processflow/processflow.h>

#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/continuum/standard.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/processflow/context.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/continuum/median.h>
#include <epic/core/log/log.h>
#include <epic/core/debug/assert.h>

#include <stdio.h>


using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CProcessFlow )


CProcessFlow::CProcessFlow()
{

}

CProcessFlow::~CProcessFlow()
{

}

Bool CProcessFlow::Process( CProcessFlowContext& ctx )
{
    const CProcessFlowContext::TTemplateCategoryList& tplCategoryList = ctx.GetTemplateCategoryList();
    if( std::find( tplCategoryList.begin(), tplCategoryList.end(), CTemplate::nCategory_Emission )==tplCategoryList.end() )
    {
        return ProcessWithoutEL( ctx );
    }
    else
    {
        return ProcessWithEL( ctx );
    }

    return false;
}


Bool CProcessFlow::ProcessWithoutEL( CProcessFlowContext& ctx )
{
    const CTemplateCatalog& templateCatalog = ctx.GetTemplateCatalog();
    const CProcessFlowContext::TTemplateCategoryList& templateCategotyList = ctx.GetTemplateCategoryList();

    for( UInt32 i=0; i<templateCategotyList.size(); i++ )
    {
        if( ctx.GetTemplateCategoryList()[i] == CTemplate::nCategory_Star )
        {
        }
        else
        {
            for( UInt32 j=0; j<templateCatalog.GetTemplateCount( (CTemplate::ECategory) templateCategotyList[i] ); j++ )
            {
                const CTemplate& tpl = templateCatalog.GetTemplate( (CTemplate::ECategory) templateCategotyList[i], j );
                const CTemplate& tplWithoutCont = templateCatalog.GetTemplateWithoutContinuum( (CTemplate::ECategory) templateCategotyList[i], j );


                BlindSolve( ctx, tpl, tplWithoutCont );
            }
        }
    }

    return true;
}


Bool CProcessFlow::ProcessWithEL( CProcessFlowContext& ctx )
{
    const CTemplateCatalog& templateCatalog = ctx.GetTemplateCatalog();
    const CProcessFlowContext::TTemplateCategoryList& templateCategotyList = ctx.GetTemplateCategoryList();

    for( UInt32 i=0; i<templateCategotyList.size(); i++ )
    {

        if( ctx.GetTemplateCategoryList()[i] == CTemplate::nCategory_Star )
        {
            Log.LogError( "Processing not yet implemented for Star template" );
        }
        else
        {
            for( UInt32 j=0; j<templateCatalog.GetTemplateCount( (CTemplate::ECategory) templateCategotyList[i] ); j++ )
            {
                const CTemplate& tpl = templateCatalog.GetTemplate( (CTemplate::ECategory) templateCategotyList[i], j );
                const CTemplate& tplWithoutCont = templateCatalog.GetTemplateWithoutContinuum( (CTemplate::ECategory) templateCategotyList[i], j );

                BlindSolve( ctx, tpl, tplWithoutCont );
            }
        }
    }

    return true;
}


Bool CProcessFlow::BlindSolve( CProcessFlowContext& ctx, const CTemplate& tpl, const CTemplate& tplWithoutCont )
{
    const CSpectrum& spc = ctx.GetSpectrum();
    const CSpectrum& spcWithoutCont = ctx.GetSpectrumWithoutContinuum();
    const TFloat64Range& lambdaRange = ctx.GetLambdaRange();

    Bool retVal = true;

    CSpectrum s = spc;
    s.GetSpectralAxis().ConvertToLogScale();

    // Create redshift initial list by spanning redshift acdross the given range, with the given delta
    TFloat64List redshifts = ctx.GetRedshiftRange().SpreadOver( ctx.GetRedshiftStep() );
    DebugAssert( redshifts.size() > 0 );

    // Compute correlation factor at each of those redshifts
    COperatorCorrelation correlation;
    correlation.Compute( spcWithoutCont, tplWithoutCont, lambdaRange, redshifts, ctx.GetOverlapThreshold() );

    // Find redshifts extremum
    Int32 nExtremums = 5;
    TPointList extremumList;
    CExtremum extremum( ctx.GetRedshiftRange() , nExtremums);
    extremum.Find( redshifts.data(), correlation.GetResults().data(), redshifts.size(), extremumList );

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
    retVal = meritChiSquare.Compute( spc, tpl, lambdaRange, extremumRedshifts, ctx.GetOverlapThreshold() );
    if( !retVal )
    {
        Log.LogInfo( "Failed to compute chi square value");
        return false;
    }

    // Store results
    ctx.AddResults( tpl,
                      extremumRedshifts, extremumCorrelation,
                      meritChiSquare.GetResults(), meritChiSquare.GetStatus(),
                      redshifts, correlation.GetResults() );


    return true;
}

