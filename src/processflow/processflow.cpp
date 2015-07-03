#include <epic/redshift/processflow/processflow.h>

#include <epic/redshift/continuum/standard.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/template/catalog.h>

#include <epic/redshift/processflow/context.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/continuum/median.h>
#include <epic/core/log/log.h>
#include <epic/core/debug/assert.h>
#include <epic/redshift/operator/peakdetection.h>
#include <epic/redshift/gaussianfit/gaussianfit.h>
#include <epic/redshift/ray/ray.h>
#include <epic/redshift/ray/catalog.h>
#include <epic/redshift/operator/raymatching.h>
#include <epic/redshift/common/median.h>

#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/chicorr.h>
#include <epic/redshift/method/blindsolveresult.h>
#include <epic/redshift/operator/raydetectionresult.h>
#include <epic/redshift/operator/raydetection.h>
#include <epic/redshift/method/blindsolve.h>
#include <epic/redshift/method/blindsolveresult.h>
#include <epic/redshift/operator/chisquare.h>

#include <epic/redshift/operator/peakdetectionresult.h>
#include <epic/redshift/operator/raydetectionresult.h>
#include <epic/redshift/operator/raymatchingresult.h>
#include <epic/redshift/operator/raymatchingresult.h>

#include <epic/redshift/method/chisquaresolve.h>
#include <epic/redshift/method/fullsolve.h>
#include <epic/redshift/method/correlationsolve.h>
#include <epic/redshift/method/linematchingsolve.h>
#include <epic/redshift/method/dtree7solve.h>
#include <epic/redshift/method/dtree7solveresult.h>
#include <epic/redshift/method/dtreeasolve.h>
#include <epic/redshift/method/dtreeasolveresult.h>
#include <epic/redshift/method/linematching2solve.h>


#include <stdio.h>
#include <float.h>

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
    if(ctx.GetParams().method  == CProcessFlowContext::nMethod_Correlation)
        return Correlation( ctx );

    if(ctx.GetParams().method  == CProcessFlowContext::nMethod_LineMatching)
        return LineMatching( ctx );

    if(ctx.GetParams().method  == CProcessFlowContext::nMethod_LineMatching2)
        return LineMatching2( ctx );

    if(ctx.GetParams().method  == CProcessFlowContext::nMethod_BlindSolve)
        return Blindsolve( ctx );

    if(ctx.GetParams().method  == CProcessFlowContext::nMethod_FullSolve)
        return Fullsolve( ctx );

    if(ctx.GetParams().method  == CProcessFlowContext::nMethod_DecisionalTree7)
        return DecisionalTree7( ctx );

    if(ctx.GetParams().method  == CProcessFlowContext::nMethod_DecisionalTreeA)
        return DecisionalTreeA( ctx );

    return false;
}


Bool CProcessFlow::Blindsolve( CProcessFlowContext& ctx, CTemplate::ECategory CategoryFilter)
{
    Log.LogInfo( "Process blindsolve (LambdaRange: %f-%f:%f)",
                 ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());


    // Remove Star category, and filter the list with regard to input variable CategoryFilter
    TTemplateCategoryList   filteredTemplateCategoryList;
    for( UInt32 i=0; i<ctx.GetParams().templateCategoryList.size(); i++ )
    {
        CTemplate::ECategory category = ctx.GetParams().templateCategoryList[i];
        if( category == CTemplate::nCategory_Star )
        {
        }
        else if(CategoryFilter == NSEpic::CTemplate::nCategory_None || CategoryFilter == category)
        {
            filteredTemplateCategoryList.push_back( category );
        }
    }


    COperatorBlindSolve blindSolve;
    CConstRef<CBlindSolveResult> blindsolveResult = blindSolve.Compute( ctx, ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), filteredTemplateCategoryList,
                                                                        ctx.GetParams().lambdaRange, ctx.GetParams().redshiftRange, ctx.GetParams().redshiftStep,
                                                                        ctx.GetParams().correlationExtremumCount, ctx.GetParams().overlapThreshold );

    if( blindsolveResult ) {
        ctx.StoreGlobalResult( "redshiftresult", *blindsolveResult );
    }

    return true;
}

Bool CProcessFlow::Correlation( CProcessFlowContext& ctx, CTemplate::ECategory CategoryFilter)
{
    Log.LogInfo( "Process Correlation (LambdaRange: %f-%f:%f)",
                 ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());


    // Remove Star category, and filter the list with regard to input variable CategoryFilter
    TTemplateCategoryList   filteredTemplateCategoryList;
    for( UInt32 i=0; i<ctx.GetParams().templateCategoryList.size(); i++ )
    {
        CTemplate::ECategory category = ctx.GetParams().templateCategoryList[i];
        if( category == CTemplate::nCategory_Star )
        {
        }
        else if(CategoryFilter == NSEpic::CTemplate::nCategory_None || CategoryFilter == category)
        {
            filteredTemplateCategoryList.push_back( category );
        }
    }

    COperatorCorrelationSolve solve;
    CConstRef<CCorrelationSolveResult> solveResult = solve.Compute( ctx, ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), filteredTemplateCategoryList,
                                                                        ctx.GetParams().lambdaRange, ctx.GetParams().redshiftRange, ctx.GetParams().redshiftStep, ctx.GetParams().overlapThreshold );

    if( solveResult ) {
        ctx.StoreGlobalResult( "redshiftresult", *solveResult );
    }

    return true;
}

Bool CProcessFlow::Fullsolve( CProcessFlowContext& ctx, CTemplate::ECategory CategoryFilter)
{
    Log.LogInfo( "Process fullsolve (LambdaRange: %f-%f:%f)",
                 ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());


    // Remove Star category, and filter the list with regard to input variable CategoryFilter
    TTemplateCategoryList   filteredTemplateCategoryList;
    for( UInt32 i=0; i<ctx.GetParams().templateCategoryList.size(); i++ )
    {
        CTemplate::ECategory category = ctx.GetParams().templateCategoryList[i];
        if( category == CTemplate::nCategory_Star )
        {
        }
        else if(CategoryFilter == NSEpic::CTemplate::nCategory_None || CategoryFilter == category)
        {
            filteredTemplateCategoryList.push_back( category );
        }
    }


    COperatorFullSolve Solve;
    CConstRef<CFullSolveResult> solveResult = Solve.Compute( ctx, ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), filteredTemplateCategoryList,
                                                                        ctx.GetParams().lambdaRange, ctx.GetParams().redshiftRange, ctx.GetParams().redshiftStep,
                                                                        ctx.GetParams().overlapThreshold );

    if( solveResult ) {
        ctx.StoreGlobalResult( "redshiftresult", *solveResult );
    }

    return true;
}


Bool CProcessFlow::LineMatching( CProcessFlowContext& ctx )
{
    Log.LogInfo( "Process Line Matching (LambdaRange: %f-%f:%f)",
            ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());

    COperatorLineMatchingSolve Solve;
    CConstRef<CLineMatchingSolveResult> solveResult = Solve.Compute(ctx, ctx.GetSpectrum(), ctx.GetParams().lambdaRange, ctx.GetParams().redshiftRange,
                                                                    ctx.GetParams().redshiftStep, ctx.GetRayCatalog() );

    if( solveResult ) {
        ctx.StoreGlobalResult( "redshiftresult", *solveResult );
    }


    return true;
}

Bool CProcessFlow::LineMatching2( CProcessFlowContext& ctx )
{
    Log.LogInfo( "Process Line Matching 2 (LambdaRange: %f-%f:%f)",
            ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());

    COperatorLineMatching2Solve Solve;
    CConstRef<CLineMatching2SolveResult> solveResult = Solve.Compute(ctx, ctx.GetSpectrum(), ctx.GetParams().lambdaRange, ctx.GetParams().redshiftRange,
                                                                    ctx.GetParams().redshiftStep, ctx.GetRayCatalog() );

    if( solveResult ) {
        ctx.StoreGlobalResult( "redshiftresult", *solveResult );
    }


    return true;
}


Bool CProcessFlow::DecisionalTree7( CProcessFlowContext& ctx )
{
    Log.LogInfo( "Process Decisional Tree 7" );


    COperatorDTree7Solve Solve;
    CConstRef<CDTree7SolveResult> solveResult = Solve.Compute( ctx, ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), ctx.GetParams().templateCategoryList, ctx.GetRayCatalog(),
                                                                        ctx.GetParams().lambdaRange, ctx.GetParams().redshiftRange, ctx.GetParams().redshiftStep,
                                                                        ctx.GetParams().correlationExtremumCount , ctx.GetParams().overlapThreshold );

    if( solveResult ) {
        ctx.StoreGlobalResult( "redshiftresult", *solveResult );
    }

    return true;
}

Bool CProcessFlow::DecisionalTreeA( CProcessFlowContext& ctx )
{
    Log.LogInfo( "Process Decisional Tree A" );


    COperatorDTreeASolve Solve;
    CConstRef<CDTreeASolveResult> solveResult = Solve.Compute( ctx, ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), ctx.GetParams().templateCategoryList, ctx.GetRayCatalog(),
                                                                        ctx.GetParams().lambdaRange, ctx.GetParams().redshiftRange, ctx.GetParams().redshiftStep,
                                                                        ctx.GetParams().correlationExtremumCount , ctx.GetParams().overlapThreshold );

    if( solveResult ) {
        ctx.StoreGlobalResult( "redshiftresult", *solveResult );
    }

    return true;
}
