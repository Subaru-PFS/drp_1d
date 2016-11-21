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
#include <epic/redshift/method/dtreebsolve.h>
#include <epic/redshift/method/dtreebsolveresult.h>
#include <epic/redshift/method/dtreecsolve.h>
#include <epic/redshift/method/dtreecsolveresult.h>
#include <epic/redshift/method/linematching2solve.h>
#include <epic/redshift/method/linemodelsolve.h>
#include <epic/redshift/method/linemodelsolveresult.h>
#include <epic/redshift/method/linemodeltplshapesolve.h>
#include <epic/redshift/method/linemodeltplshapesolveresult.h>

#include <boost/algorithm/string.hpp>
#include <stdio.h>
#include <float.h>

using namespace NSEpic;


CProcessFlow::CProcessFlow()
{

}

CProcessFlow::~CProcessFlow()
{

}

Bool CProcessFlow::Process( CProcessFlowContext& ctx )
{
    if(1)
    {
        //Check if the Spectrum is valid on the lambdarange
        TFloat64Range lambdaRange;
        ctx.GetParameterStore().Get( "lambdaRange", lambdaRange );
        const CSpectrumSpectralAxis& spcSpectralAxis = ctx.GetSpectrum().GetSpectralAxis();
        TFloat64Range spcLambdaRange;
        spcSpectralAxis.ClampLambdaRange( lambdaRange, spcLambdaRange );
        const Float64 lmin = spcLambdaRange.GetBegin();
        const Float64 lmax = spcLambdaRange.GetEnd();
        if( !ctx.GetSpectrum().IsFluxValid( lmin, lmax ) ){
            Log.LogError( "Failed to validate spectrum flux: %s, on wavelength range (%.1f ; %.1f)", ctx.GetSpectrum().GetName().c_str(), lmin, lmax );
            return false;
        }
        if( !ctx.GetSpectrum().IsNoiseValid( lmin, lmax ) ){
            Log.LogError( "Failed to validate noise from spectrum: %s, on wavelength range (%.1f ; %.1f)", ctx.GetSpectrum().GetName().c_str(), lmin, lmax );
            return false;
        }
    }



    std::string methodName;
    ctx.GetParameterStore().Get( "method", methodName );

    boost::algorithm::to_lower(methodName);

    if(methodName  == "correlationsolve" )
        return Correlation( ctx );

    if(methodName  == "chisquare2solve" )
        return Chisquare( ctx );

    if(methodName  == "linematching" )
        return LineMatching( ctx );

    if(methodName  == "linematching2" )
        return LineMatching2( ctx );

    if(methodName  == "linemodel" )
        return LineModelSolve( ctx );

    if(methodName  == "linemodeltplshape" )
        return LineModelTplshapeSolve( ctx );

    if(methodName  == "blindsolve" )
        return Blindsolve( ctx );

    if(methodName  == "fullsolve" )
        return Fullsolve( ctx );

    if(methodName  == "amazed0_1" )
        return DecisionalTree7( ctx );

    if(methodName  == "decisionaltreea" )
        return DecisionalTreeA( ctx );

    if(methodName  == "amazed0_2" )
        return DecisionalTreeB( ctx );

    if(methodName  == "amazed0_3" )
        return DecisionalTreeC( ctx );

    Log.LogError("Problem found while parsing the method parameter !");
    return false;
}


Bool CProcessFlow::Blindsolve( CProcessFlowContext& ctx, const std::string&  CategoryFilter)
{
    Log.LogInfo( "Process blindsolve (LambdaRange: %f-%f:%f)",
                 ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());


    // Remove Star category, and filter the list with regard to input variable CategoryFilter
    TStringList templateCategoryList;
    ctx.GetParameterStore().Get( "templateCategoryList", templateCategoryList );
    TStringList   filteredTemplateCategoryList;
    for( UInt32 i=0; i<templateCategoryList.size(); i++ )
    {
        std::string category = templateCategoryList[i];
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

    ctx.GetParameterStore().Get( "lambdaRange", lambdaRange );
    ctx.GetParameterStore().Get( "redshiftRange", redshiftRange );
    ctx.GetParameterStore().Get( "redshiftStep", redshiftStep, 0.0001 );


    COperatorBlindSolve blindSolve;
    std::shared_ptr<const CBlindSolveResult> blindsolveResult = blindSolve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), filteredTemplateCategoryList,
                                                                        lambdaRange, redshiftRange, redshiftStep);

    if( blindsolveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", blindsolveResult );
    }

    return true;
}

Bool CProcessFlow::Correlation( CProcessFlowContext& ctx,  const std::string&  CategoryFilter)
{
    Log.LogInfo( "Process Correlation (LambdaRange: %f-%f:%f)",
                 ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());


    // Remove Star category, and filter the list with regard to input variable CategoryFilter
    TStringList templateCategoryList;
    ctx.GetParameterStore().Get( "templateCategoryList", templateCategoryList );
    TStringList   filteredTemplateCategoryList;
    for( UInt32 i=0; i<templateCategoryList.size(); i++ )
    {
        std::string category = templateCategoryList[i];
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

    ctx.GetParameterStore().Get( "lambdaRange", lambdaRange );
    ctx.GetParameterStore().Get( "redshiftRange", redshiftRange );
    ctx.GetParameterStore().Get( "redshiftStep", redshiftStep );

    COperatorCorrelationSolve solve;
    std::shared_ptr<const CCorrelationSolveResult> solveResult = solve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), filteredTemplateCategoryList,
                                                                        lambdaRange, redshiftRange, redshiftStep );

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", solveResult );
    }

    return true;
}


Bool CProcessFlow::Chisquare( CProcessFlowContext& ctx, const std::string& CategoryFilter)
{
    TFloat64Range lambdaRange;
    TFloat64Range redshiftRange;
    Float64       redshiftStep;

    ctx.GetParameterStore().Get( "lambdaRange", lambdaRange );
    ctx.GetParameterStore().Get( "redshiftRange", redshiftRange );
    ctx.GetParameterStore().Get( "redshiftStep", redshiftStep );

    const CSpectrumSpectralAxis& spcSpectralAxis = ctx.GetSpectrum().GetSpectralAxis();
    TFloat64Range spcLambdaRange;
    spcSpectralAxis.ClampLambdaRange( lambdaRange, spcLambdaRange );

    Log.LogInfo( "Process Chisquare for spc:%s (LambdaRange: %f-%f:%f)", ctx.GetSpectrum().GetName().c_str(),
                 spcLambdaRange.GetBegin(), spcLambdaRange.GetEnd(), ctx.GetSpectrum().GetResolution());


    // Remove Star category, and filter the list with regard to input variable CategoryFilter
    TStringList templateCategoryList;
    ctx.GetParameterStore().Get( "templateCategoryList", templateCategoryList );
    TStringList   filteredTemplateCategoryList;
    for( UInt32 i=0; i<templateCategoryList.size(); i++ )
    {
        std::string category = templateCategoryList[i];
        if( category == "star" )
        {
        }
        else if(CategoryFilter == "all" || CategoryFilter == category)
        {
            filteredTemplateCategoryList.push_back( category );
        }
    }



    // Create redshift initial list by spanning redshift acdross the given range, with the given delta
    TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
    DebugAssert( redshifts.size() > 0 );

    Float64 overlapThreshold;
    ctx.GetParameterStore().Get( "chisquare2solve.overlapThreshold", overlapThreshold, 1.0);
    std::string opt_spcComponent;
    ctx.GetDataStore().GetScopedParam( "chisquare2solve.spectrum.component", opt_spcComponent, "raw" );
    std::string opt_interp;
    ctx.GetDataStore().GetScopedParam( "chisquare2solve.interpolation", opt_interp, "precomputedfinegrid" );
    std::string opt_extinction;
    ctx.GetDataStore().GetScopedParam( "chisquare2solve.extinction", opt_extinction, "no" );

    // prepare the unused masks
    std::vector<CMask> maskList;

    CMethodChisquare2Solve solve;
    std::shared_ptr< const CChisquare2SolveResult> solveResult = solve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), filteredTemplateCategoryList,
                                                                        spcLambdaRange, redshifts, overlapThreshold, maskList, opt_spcComponent, opt_interp, opt_extinction);

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", solveResult );
    }

    return true;
}

Bool CProcessFlow::Fullsolve( CProcessFlowContext& ctx, const std::string& CategoryFilter)
{
    Log.LogInfo( "Process fullsolve (LambdaRange: %f-%f:%f)",
                 ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());


    // Remove Star category, and filter the list with regard to input variable CategoryFilter
    TStringList templateCategoryList;
    ctx.GetParameterStore().Get( "templateCategoryList", templateCategoryList );
    TStringList   filteredTemplateCategoryList;
    for( UInt32 i=0; i<templateCategoryList.size(); i++ )
    {
        std::string category = templateCategoryList[i];
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
    Float64       overlapThreshold;

    ctx.GetParameterStore().Get( "lambdaRange", lambdaRange );
    ctx.GetParameterStore().Get( "redshiftRange", redshiftRange );
    ctx.GetParameterStore().Get( "redshiftStep", redshiftStep );
    ctx.GetParameterStore().Get( "overlapThreshold", overlapThreshold );


    COperatorFullSolve Solve;
    std::shared_ptr<const CFullSolveResult> solveResult = Solve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), filteredTemplateCategoryList,
                                                                        lambdaRange, redshiftRange, redshiftStep, overlapThreshold );

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", solveResult );
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
    std::shared_ptr<const CLineMatchingSolveResult> solveResult = Solve.Compute(ctx.GetDataStore(), ctx.GetSpectrum(),
                                                                    lambdaRange, redshiftRange,
                                                                    redshiftStep, ctx.GetRayCatalog() );

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", solveResult );
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
    std::shared_ptr<const CLineMatching2SolveResult> solveResult = Solve.Compute(ctx.GetDataStore(), ctx.GetSpectrum(), spcLambdaRange, redshiftRange,
                                                                    redshiftStep, ctx.GetRayCatalog() );

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", solveResult );
    }

    return true;
}

Bool CProcessFlow::LineModelSolve( CProcessFlowContext& ctx )
{
    TFloat64Range lambdaRange;
    TFloat64Range redshiftRange;
    Float64       redshiftStep;

    std::string CategoryFilter="all";
    // Remove Star category, and filter the list with regard to input variable CategoryFilter
    TStringList templateCategoryList;
    ctx.GetParameterStore().Get( "templateCategoryList", templateCategoryList );
    TStringList   filteredTemplateCategoryList;
    for( UInt32 i=0; i<templateCategoryList.size(); i++ )
    {
        std::string category = templateCategoryList[i];
        if( category == "star" )
        {
        }
        else if(CategoryFilter == "all" || CategoryFilter == category)
        {
            filteredTemplateCategoryList.push_back( category );
        }
    }

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
    std::shared_ptr<const CLineModelSolveResult> solveResult = Solve.Compute( ctx.GetDataStore(),
                                          ctx.GetSpectrum(),
                                          ctx.GetSpectrumWithoutContinuum(),
                                          ctx.GetTemplateCatalog(),
                                          filteredTemplateCategoryList,
                                          ctx.GetRayCatalog(),
                                          spcLambdaRange,
                                          redshifts );

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", solveResult );
    }

    return true;
}

Bool CProcessFlow::LineModelTplshapeSolve( CProcessFlowContext& ctx, const std::string& CategoryFilter )
{
    TFloat64Range lambdaRange;
    TFloat64Range redshiftRange;
    Float64       redshiftStep;

    ctx.GetParameterStore().Get( "lambdaRange", lambdaRange );
    ctx.GetParameterStore().Get( "redshiftRange", redshiftRange );
    ctx.GetParameterStore().Get( "redshiftStep", redshiftStep );

    TFloat64Range spcLambdaRange;
    ctx.GetSpectrum().GetSpectralAxis().ClampLambdaRange( lambdaRange, spcLambdaRange );

    Log.LogInfo( "Processing LineModel-TplShape for spc:%s (LambdaRange: %f-%f:%f)", ctx.GetSpectrum().GetName().c_str(),
            spcLambdaRange.GetBegin(), spcLambdaRange.GetEnd(), ctx.GetSpectrum().GetResolution());

    // Create redshift initial list by spanning redshift acdross the given range, with the given delta
    TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
    DebugAssert( redshifts.size() > 0 );

    // Remove Star category, and filter the list with regard to input variable CategoryFilter
    TStringList templateCategoryList;
    ctx.GetParameterStore().Get( "templateCategoryList", templateCategoryList );
    TStringList   filteredTemplateCategoryList;
    for( UInt32 i=0; i<templateCategoryList.size(); i++ )
    {
        std::string category = templateCategoryList[i];
        if( category == "star" )
        {
        }
        else if(CategoryFilter == "all" || CategoryFilter == category)
        {
            filteredTemplateCategoryList.push_back( category );
        }
    }

    CLineModelTplshapeSolve Solve;
    std::shared_ptr<const CLineModelTplshapeSolveResult> solveResult = Solve.Compute( ctx.GetDataStore(),
                                          ctx.GetSpectrum(),
                                          ctx.GetSpectrumWithoutContinuum(),
                                          ctx.GetTemplateCatalog(),
                                          filteredTemplateCategoryList,
                                          ctx.GetRayCatalog(),
                                          spcLambdaRange,
                                          redshifts );

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", solveResult );
    }

    return true;
}

Bool CProcessFlow::DecisionalTree7( CProcessFlowContext& ctx )
{
    Log.LogInfo( "Process Decisional Tree 7" );

    TFloat64Range lambdaRange;
    TFloat64Range redshiftRange;
    Float64       redshiftStep;
    TStringList     templateCategoryList;

    ctx.GetParameterStore().Get( "lambdaRange", lambdaRange );
    ctx.GetParameterStore().Get( "redshiftRange", redshiftRange );
    ctx.GetParameterStore().Get( "redshiftStep", redshiftStep );
    ctx.GetParameterStore().Get( "templateCategoryList", templateCategoryList );

    COperatorDTree7Solve Solve;
    std::shared_ptr<CDTree7SolveResult> solveResult = Solve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), templateCategoryList, ctx.GetRayCatalog(),
                                                                        lambdaRange, redshiftRange, redshiftStep);

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", solveResult );
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
    std::shared_ptr<const CDTreeASolveResult> solveResult = Solve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), templateCategoryList, ctx.GetRayCatalog(),
                                                                        lambdaRange, redshiftRange, redshiftStep,
                                                                        correlationExtremumCount , overlapThreshold );

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", solveResult );
    }

    return true;
}


Bool CProcessFlow::DecisionalTreeB( CProcessFlowContext& ctx )
{
    TFloat64Range lambdaRange;
    TFloat64Range redshiftRange;
    Float64       redshiftStep;
    TStringList     templateCategoryList;

    ctx.GetParameterStore().Get( "lambdaRange", lambdaRange );
    ctx.GetParameterStore().Get( "redshiftRange", redshiftRange );
    ctx.GetParameterStore().Get( "redshiftStep", redshiftStep );
    ctx.GetParameterStore().Get( "templateCategoryList", templateCategoryList );

    const CSpectrumSpectralAxis& spcSpectralAxis = ctx.GetSpectrum().GetSpectralAxis();
    TFloat64Range spcLambdaRange;
    spcSpectralAxis.ClampLambdaRange( lambdaRange, spcLambdaRange );

    Log.LogInfo( "Processing dtreeb for spc:%s (LambdaRange: %.2f-%.2f:%.2f)", ctx.GetSpectrum().GetName().c_str(),
            spcLambdaRange.GetBegin(), spcLambdaRange.GetEnd(), ctx.GetSpectrum().GetResolution());

    // Create redshift initial list by spanning redshift acdross the given range, with the given delta
    TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
    DebugAssert( redshifts.size() > 0 );

    COperatorDTreeBSolve Solve;
    std::shared_ptr<const CDTreeBSolveResult> solveResult = Solve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), templateCategoryList, ctx.GetRayCatalog(),
                                                                        spcLambdaRange, redshifts);

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", solveResult );
    }

    return true;
}

Bool CProcessFlow::DecisionalTreeC( CProcessFlowContext& ctx )
{
    TFloat64Range lambdaRange;
    TFloat64Range redshiftRange;
    Float64       redshiftStep;
    TStringList     templateCategoryList;

    ctx.GetParameterStore().Get( "lambdaRange", lambdaRange );
    ctx.GetParameterStore().Get( "redshiftRange", redshiftRange );
    ctx.GetParameterStore().Get( "redshiftStep", redshiftStep );
    ctx.GetParameterStore().Get( "templateCategoryList", templateCategoryList );

    const CSpectrumSpectralAxis& spcSpectralAxis = ctx.GetSpectrum().GetSpectralAxis();
    TFloat64Range spcLambdaRange;
    spcSpectralAxis.ClampLambdaRange( lambdaRange, spcLambdaRange );

    Log.LogInfo( "Processing dtreec for spc:%s (LambdaRange: %.2f-%.2f:%.2f)", ctx.GetSpectrum().GetName().c_str(),
            spcLambdaRange.GetBegin(), spcLambdaRange.GetEnd(), ctx.GetSpectrum().GetResolution());

    // Create redshift initial list by spanning redshift acdross the given range, with the given delta
    TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
    DebugAssert( redshifts.size() > 0 );

    COperatorDTreeCSolve Solve;
    std::shared_ptr<const CDTreeCSolveResult> solveResult = Solve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), templateCategoryList, ctx.GetRayCatalog(),
                                                                        spcLambdaRange, redshifts);

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", solveResult );
    }

    return true;
}
