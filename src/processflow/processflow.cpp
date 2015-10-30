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
#include <epic/redshift/method/chisquare2solve.h>
#include <epic/redshift/method/fullsolve.h>
#include <epic/redshift/method/correlationsolve.h>
#include <epic/redshift/method/linematchingsolve.h>
#include <epic/redshift/method/dtree7solve.h>
#include <epic/redshift/method/dtree7solveresult.h>
#include <epic/redshift/method/dtreeasolve.h>
#include <epic/redshift/method/dtreeasolveresult.h>
#include <epic/redshift/method/linematching2solve.h>
#include <epic/redshift/method/linemodelsolve.h>
#include <epic/redshift/method/linemodelsolveresult.h>

#include <boost/algorithm/string.hpp>
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
    std::string methodName;
    ctx.GetParameterStore().Get( "method", methodName );

    boost::algorithm::to_lower(methodName);

    if(methodName  == "correlation" )
        return Correlation( ctx );

    if(methodName  == "chisquare" )
        return Chisquare( ctx );

    if(methodName  == "linematching" )
        return LineMatching( ctx );

    if(methodName  == "linematching2" )
        return LineMatching2( ctx );

    if(methodName  == "linemodel" )
        return LineModelSolve( ctx );

    if(methodName  == "blindsolve" )
        return Blindsolve( ctx );

    if(methodName  == "fullsolve" )
        return Fullsolve( ctx );

    if(methodName  == "decisionaltree7" )
        return DecisionalTree7( ctx );

    if(methodName  == "decisionaltreeA" )
        return DecisionalTreeA( ctx );

    return false;
}


Bool CProcessFlow::Blindsolve( CProcessFlowContext& ctx, const std::string&  CategoryFilter)
{
    Log.LogInfo( "Process blindsolve (LambdaRange: %f-%f:%f)",
                 ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());


    // Remove Star category, and filter the list with regard to input variable CategoryFilter
    TStringList tempalteCategoyList;
    ctx.GetParameterStore().Get( "templateCategoryList", tempalteCategoyList );
    TStringList   filteredTemplateCategoryList;
    for( UInt32 i=0; i<tempalteCategoyList.size(); i++ )
    {
        std::string category = tempalteCategoyList[i];
        if( category == "star" )
        {
        }
        else if(CategoryFilter == "all" || CategoryFilter == category)
        {
            filteredTemplateCategoryList.push_back( category );
        }
    }

    TFloat64Range lambdaRange;
    TFloat64Range redshiftRange;
    Float64       redshiftStep;
    Int64         correlationExtremumCount;
    Float64       overlapThreshold;

    ctx.GetParameterStore().Get( "lambdaRange", lambdaRange );
    ctx.GetParameterStore().Get( "redshiftRange", redshiftRange );
    ctx.GetParameterStore().Get( "redshiftStep", redshiftStep, 0.001 );
    ctx.GetParameterStore().Get( "correlationExtremumCount", correlationExtremumCount, 5.0 );
    ctx.GetParameterStore().Get( "overlapThreshold", overlapThreshold, 1.0 );

    COperatorBlindSolve blindSolve;
    CConstRef<CBlindSolveResult> blindsolveResult = blindSolve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), filteredTemplateCategoryList,
                                                                        lambdaRange, redshiftRange, redshiftStep,
                                                                        (Int32)correlationExtremumCount, overlapThreshold );

    if( blindsolveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", *blindsolveResult );
    }

    return true;
}

Bool CProcessFlow::Correlation( CProcessFlowContext& ctx,  const std::string&  CategoryFilter)
{
    Log.LogInfo( "Process Correlation (LambdaRange: %f-%f:%f)",
                 ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());


    // Remove Star category, and filter the list with regard to input variable CategoryFilter
    TStringList tempalteCategoyList;
    ctx.GetParameterStore().Get( "tempalteCategoyList", tempalteCategoyList );
    TStringList   filteredTemplateCategoryList;
    for( UInt32 i=0; i<tempalteCategoyList.size(); i++ )
    {
        std::string category = tempalteCategoyList[i];
        if( category == "star" )
        {
        }
        else if(CategoryFilter == "all" || CategoryFilter == category)
        {
            filteredTemplateCategoryList.push_back( category );
        }
    }

    TFloat64Range lambdaRange;
    TFloat64Range redshiftRange;
    Float64       redshiftStep;
    Int64         correlationExtremumCount;
    Float64       overlapThreshold;

    ctx.GetParameterStore().Get( "lambdaRange", lambdaRange );
    ctx.GetParameterStore().Get( "redshiftRange", redshiftRange );
    ctx.GetParameterStore().Get( "redshiftStep", redshiftStep );
    ctx.GetParameterStore().Get( "correlationExtremumCount", correlationExtremumCount );
    ctx.GetParameterStore().Get( "overlapThreshold", overlapThreshold );

    COperatorCorrelationSolve solve;
    CConstRef<CCorrelationSolveResult> solveResult = solve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), filteredTemplateCategoryList,
                                                                        lambdaRange, redshiftRange, redshiftStep, overlapThreshold );

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", *solveResult );
    }

    return true;
}


Bool CProcessFlow::Chisquare( CProcessFlowContext& ctx, const std::string& CategoryFilter)
{
    Log.LogInfo( "Process Chisquare (LambdaRange: %f-%f:%f)",
                 ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());


    // Remove Star category, and filter the list with regard to input variable CategoryFilter
    TStringList tempalteCategoyList;
    ctx.GetParameterStore().Get( "tempalteCategoyList", tempalteCategoyList );
    TStringList   filteredTemplateCategoryList;
    for( UInt32 i=0; i<tempalteCategoyList.size(); i++ )
    {
        std::string category = tempalteCategoyList[i];
        if( category == "star" )
        {
        }
        else if(CategoryFilter == "all" || CategoryFilter == category)
        {
            filteredTemplateCategoryList.push_back( category );
        }
    }

    TFloat64Range lambdaRange;
    TFloat64Range redshiftRange;
    Float64       redshiftStep;
    Int64         correlationExtremumCount;
    Float64       overlapThreshold;

    ctx.GetParameterStore().Get( "lambdaRange", lambdaRange );
    ctx.GetParameterStore().Get( "redshiftRange", redshiftRange );
    ctx.GetParameterStore().Get( "redshiftStep", redshiftStep );
    ctx.GetParameterStore().Get( "correlationExtremumCount", correlationExtremumCount );
    ctx.GetParameterStore().Get( "overlapThreshold", overlapThreshold );

    // Create redshift initial list by spanning redshift acdross the given range, with the given delta
    TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
    DebugAssert( redshifts.size() > 0 );

    CMethodChisquare2Solve solve;
    CConstRef<CChisquare2SolveResult> solveResult = solve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), filteredTemplateCategoryList,
                                                                        lambdaRange, redshifts, overlapThreshold );

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", *solveResult );
    }

    return true;
}

Bool CProcessFlow::Fullsolve( CProcessFlowContext& ctx, const std::string& CategoryFilter)
{
    Log.LogInfo( "Process fullsolve (LambdaRange: %f-%f:%f)",
                 ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());


    // Remove Star category, and filter the list with regard to input variable CategoryFilter
    TStringList tempalteCategoyList;
    ctx.GetParameterStore().Get( "tempalteCategoyList", tempalteCategoyList );
    TStringList   filteredTemplateCategoryList;
    for( UInt32 i=0; i<tempalteCategoyList.size(); i++ )
    {
        std::string category = tempalteCategoyList[i];
        if( category == "star" )
        {
        }
        else if(CategoryFilter == "all" || CategoryFilter == category)
        {
            filteredTemplateCategoryList.push_back( category );
        }
    }

    TFloat64Range lambdaRange;
    TFloat64Range redshiftRange;
    Float64       redshiftStep;
    Int64         correlationExtremumCount;
    Float64       overlapThreshold;

    ctx.GetParameterStore().Get( "lambdaRange", lambdaRange );
    ctx.GetParameterStore().Get( "redshiftRange", redshiftRange );
    ctx.GetParameterStore().Get( "redshiftStep", redshiftStep );
    ctx.GetParameterStore().Get( "correlationExtremumCount", correlationExtremumCount );
    ctx.GetParameterStore().Get( "overlapThreshold", overlapThreshold );


    COperatorFullSolve Solve;
    CConstRef<CFullSolveResult> solveResult = Solve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), filteredTemplateCategoryList,
                                                                        lambdaRange, redshiftRange, redshiftStep, overlapThreshold );

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", *solveResult );
    }

    return true;
}


Bool CProcessFlow::LineMatching( CProcessFlowContext& ctx )
{
    Log.LogInfo( "Processing Line Matching (LambdaRange: %f-%f:%f)",
            ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());

    TFloat64Range lambdaRange;
    TFloat64Range redshiftRange;
    Float64       redshiftStep;

    ctx.GetParameterStore().Get( "lambdaRange", lambdaRange );
    ctx.GetParameterStore().Get( "redshiftRange", redshiftRange );
    ctx.GetParameterStore().Get( "redshiftStep", redshiftStep );

    COperatorLineMatchingSolve Solve;
    CConstRef<CLineMatchingSolveResult> solveResult = Solve.Compute(ctx.GetDataStore(), ctx.GetSpectrum(),
                                                                    lambdaRange, redshiftRange,
                                                                    redshiftStep, ctx.GetRayCatalog() );

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", *solveResult );
    }


    return true;
}

Bool CProcessFlow::LineMatching2( CProcessFlowContext& ctx )
{
    TFloat64Range lambdaRange;
    TFloat64Range redshiftRange;
    Float64       redshiftStep;

    ctx.GetParameterStore().Get( "lambdaRange", lambdaRange );
    ctx.GetParameterStore().Get( "redshiftRange", redshiftRange );
    ctx.GetParameterStore().Get( "redshiftStep", redshiftStep );

    TFloat64Range spcLambdaRange;
    ctx.GetSpectrum().GetSpectralAxis().ClampLambdaRange( lambdaRange, spcLambdaRange );

    Log.LogInfo( "Processing Line Matching 2 (LambdaRange: %f-%f:%f)",
            spcLambdaRange.GetBegin(), spcLambdaRange.GetEnd(), ctx.GetSpectrum().GetResolution());

    COperatorLineMatching2Solve Solve;
    CConstRef<CLineMatching2SolveResult> solveResult = Solve.Compute(ctx.GetDataStore(), ctx.GetSpectrum(), spcLambdaRange, redshiftRange,
                                                                    redshiftStep, ctx.GetRayCatalog() );

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", *solveResult );
    }

    return true;
}

Bool CProcessFlow::LineModelSolve( CProcessFlowContext& ctx )
{
    TFloat64Range lambdaRange;
    TFloat64Range redshiftRange;
    Float64       redshiftStep;

    ctx.GetParameterStore().Get( "lambdaRange", lambdaRange );
    ctx.GetParameterStore().Get( "redshiftRange", redshiftRange );
    ctx.GetParameterStore().Get( "redshiftStep", redshiftStep );

    TFloat64Range spcLambdaRange;
    ctx.GetSpectrum().GetSpectralAxis().ClampLambdaRange( lambdaRange, spcLambdaRange );

    Log.LogInfo( "Processing Line Model for spc:%s (LambdaRange: %f-%f:%f)", ctx.GetSpectrum().GetName().c_str(),
            spcLambdaRange.GetBegin(), spcLambdaRange.GetEnd(), ctx.GetSpectrum().GetResolution());

    // Create redshift initial list by spanning redshift acdross the given range, with the given delta
    TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
    DebugAssert( redshifts.size() > 0 );

    CLineModelSolve Solve;
    CConstRef<CLineModelSolveResult> solveResult = Solve.Compute(ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(), ctx.GetRayCatalog(),
                                                                 spcLambdaRange, redshifts);

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", *solveResult );
    }


    return true;
}

Bool CProcessFlow::DecisionalTree7( CProcessFlowContext& ctx )
{
    Log.LogInfo( "Process Decisional Tree 7" );

    TFloat64Range lambdaRange;
    TFloat64Range redshiftRange;
    Float64       redshiftStep;
    Int64         correlationExtremumCount;
    Float64       overlapThreshold;
    TStringList     templateCategoryList;

    ctx.GetParameterStore().Get( "lambdaRange", lambdaRange );
    ctx.GetParameterStore().Get( "redshiftRange", redshiftRange );
    ctx.GetParameterStore().Get( "redshiftStep", redshiftStep );
    ctx.GetParameterStore().Get( "correlationExtremumCount", correlationExtremumCount );
    ctx.GetParameterStore().Get( "overlapThreshold", overlapThreshold );
    ctx.GetParameterStore().Get( "templateCategoryList", templateCategoryList );

    COperatorDTree7Solve Solve;
    CConstRef<CDTree7SolveResult> solveResult = Solve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), templateCategoryList, ctx.GetRayCatalog(),
                                                                        lambdaRange, redshiftRange, redshiftStep,
                                                                        correlationExtremumCount , overlapThreshold );

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", *solveResult );
    }

    return true;
}

Bool CProcessFlow::DecisionalTreeA( CProcessFlowContext& ctx )
{
    Log.LogInfo( "Process Decisional Tree A" );

    TFloat64Range lambdaRange;
    TFloat64Range redshiftRange;
    Float64       redshiftStep;
    Int64         correlationExtremumCount;
    Float64       overlapThreshold;
    TStringList     templateCategoryList;

    ctx.GetParameterStore().Get( "lambdaRange", lambdaRange );
    ctx.GetParameterStore().Get( "redshiftRange", redshiftRange );
    ctx.GetParameterStore().Get( "redshiftStep", redshiftStep );
    ctx.GetParameterStore().Get( "correlationExtremumCount", correlationExtremumCount );
    ctx.GetParameterStore().Get( "overlapThreshold", overlapThreshold );
    ctx.GetParameterStore().Get( "templateCategoryList", templateCategoryList );


    COperatorDTreeASolve Solve;
    CConstRef<CDTreeASolveResult> solveResult = Solve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), templateCategoryList, ctx.GetRayCatalog(),
                                                                        lambdaRange, redshiftRange, redshiftStep,
                                                                        correlationExtremumCount , overlapThreshold );

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", *solveResult );
    }

    return true;
}
