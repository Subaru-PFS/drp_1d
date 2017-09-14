#include <RedshiftLibrary/processflow/processflow.h>

#include <RedshiftLibrary/continuum/standard.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>

#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/continuum/median.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/operator/peakdetection.h>
#include <RedshiftLibrary/gaussianfit/gaussianfit.h>
#include <RedshiftLibrary/ray/ray.h>
#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/operator/raymatching.h>
#include <RedshiftLibrary/common/median.h>

#include <RedshiftLibrary/operator/correlation.h>
#include <RedshiftLibrary/operator/chicorr.h>
#include <RedshiftLibrary/method/blindsolveresult.h>
#include <RedshiftLibrary/operator/raydetectionresult.h>
#include <RedshiftLibrary/operator/raydetection.h>
#include <RedshiftLibrary/method/blindsolve.h>
#include <RedshiftLibrary/method/blindsolveresult.h>
#include <RedshiftLibrary/operator/chisquare.h>

#include <RedshiftLibrary/operator/peakdetectionresult.h>
#include <RedshiftLibrary/operator/raydetectionresult.h>
#include <RedshiftLibrary/operator/raymatchingresult.h>
#include <RedshiftLibrary/operator/raymatchingresult.h>

#include <RedshiftLibrary/method/chisquaresolve.h>
#include <RedshiftLibrary/method/chisquare2solve.h>
#include <RedshiftLibrary/method/chisquarelogsolve.h>
#include <RedshiftLibrary/method/fullsolve.h>
#include <RedshiftLibrary/method/correlationsolve.h>
#include <RedshiftLibrary/method/linematchingsolve.h>
#include <RedshiftLibrary/method/dtree7solve.h>
#include <RedshiftLibrary/method/dtree7solveresult.h>
#include <RedshiftLibrary/method/dtreeasolve.h>
#include <RedshiftLibrary/method/dtreeasolveresult.h>
#include <RedshiftLibrary/method/dtreebsolve.h>
#include <RedshiftLibrary/method/dtreebsolveresult.h>
#include <RedshiftLibrary/method/dtreecsolve.h>
#include <RedshiftLibrary/method/dtreecsolveresult.h>
#include <RedshiftLibrary/method/linematching2solve.h>
#include <RedshiftLibrary/method/linemodelsolve.h>
#include <RedshiftLibrary/method/linemodelsolveresult.h>
#include <RedshiftLibrary/method/linemodeltplshapesolve.h>
#include <RedshiftLibrary/method/linemodeltplshapesolveresult.h>

#include <RedshiftLibrary/reliability/zqual.h>

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
    Log.LogInfo("<proc-spc><%s>", ctx.GetSpectrum().GetName().c_str());

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

    if(methodName  == "chisquarelogsolve" )
    {
        return ChisquareLog( ctx );
    }

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
    }else{
        return false;
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
    }else{
        return false;
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
    std::string opt_dustFit;
    ctx.GetDataStore().GetScopedParam( "chisquare2solve.dustfit", opt_dustFit, "no" );

    Log.LogInfo( "Process Chisquare using overlapThreshold: %.3f", overlapThreshold);
    Log.LogInfo( "Process Chisquare using component: %s", opt_spcComponent.c_str());
    Log.LogInfo( "Process Chisquare using extinction: %s", opt_extinction.c_str());
    Log.LogInfo( "Process Chisquare using dust-fit: %s", opt_dustFit.c_str());
    Log.LogInfo( "..." );

    // prepare the unused masks
    std::vector<CMask> maskList;
    //retrieve the calibration dir path
    std::string calibrationDirPath;
    ctx.GetParameterStore().Get( "calibrationDir", calibrationDirPath );
    CMethodChisquare2Solve solve(calibrationDirPath);
    std::shared_ptr< const CChisquare2SolveResult> solveResult = solve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), filteredTemplateCategoryList,
                                                                        spcLambdaRange, redshifts, overlapThreshold, maskList, opt_spcComponent, opt_interp, opt_extinction, opt_dustFit);

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", solveResult );
    }else{
        return false;
    }
    return true;
}

Bool CProcessFlow::ChisquareLog( CProcessFlowContext& ctx, const std::string& CategoryFilter)
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

    Log.LogInfo( "Process ChisquareLog for spc:%s (LambdaRange: %f-%f:%f)", ctx.GetSpectrum().GetName().c_str(),
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
    ctx.GetParameterStore().Get( "chisquarelogsolve.overlapThreshold", overlapThreshold, 1.0);
    std::string opt_spcComponent;
    ctx.GetDataStore().GetScopedParam( "chisquarelogsolve.spectrum.component", opt_spcComponent, "raw" );
    std::string opt_interp;
    ctx.GetDataStore().GetScopedParam( "chisquarelogsolve.interpolation", opt_interp, "precomputedfinegrid" );
    std::string opt_extinction;
    ctx.GetDataStore().GetScopedParam( "chisquarelogsolve.extinction", opt_extinction, "no" );
    std::string opt_dustFit;
    ctx.GetDataStore().GetScopedParam( "chisquarelogsolve.dustfit", opt_dustFit, "no" );

    Log.LogInfo( "Process Chisquare using overlapThreshold: %.3f", overlapThreshold);
    Log.LogInfo( "Process Chisquare using component: %s", opt_spcComponent.c_str());
    Log.LogInfo( "Process Chisquare using extinction: %s", opt_extinction.c_str());
    Log.LogInfo( "Process Chisquare using dust-fit: %s", opt_dustFit.c_str());
    Log.LogInfo( "..." );

    // prepare the unused masks
    std::vector<CMask> maskList;
    //retrieve the calibration dir path
    std::string calibrationDirPath;
    ctx.GetParameterStore().Get( "calibrationDir", calibrationDirPath );
    CMethodChisquareLogSolve solve(calibrationDirPath);
    std::shared_ptr< const CChisquareLogSolveResult> solveResult = solve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), filteredTemplateCategoryList,
                                                                        spcLambdaRange, redshifts, overlapThreshold, maskList, opt_spcComponent, opt_interp, opt_extinction, opt_dustFit);

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", solveResult );
    }else{
        return false;
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
    }else{
        return false;
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
    }else{
        return false;
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
    }else{
        return false;
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


    //retrieve the calibration dir path
    std::string calibrationDirPath;
    ctx.GetParameterStore().Get( "calibrationDir", calibrationDirPath );
    CLineModelSolve Solve(calibrationDirPath);
    std::shared_ptr<CLineModelSolveResult> solveResult = Solve.Compute( ctx.GetDataStore(),
                                          ctx.GetSpectrum(),
                                          ctx.GetSpectrumWithoutContinuum(),
                                          ctx.GetTemplateCatalog(),
                                          filteredTemplateCategoryList,
                                          ctx.GetRayCatalog(),
                                          spcLambdaRange,
                                          redshifts );

    //*
    // todo: call the qualz object as done in zbayes branch from S. Jamal.
    bool enableQualz = true;
    if(!solveResult && enableQualz)
    {
        Log.LogWarning( "Reliability skipped - no redshift results found");
        enableQualz = false;
    }
    if ( enableQualz && isPdfValid(ctx) )
    {
        CClassifierStore classifStore = ctx.GetClassifierStore();
        if(!classifStore.m_isInitialized)
        {
            Log.LogWarning( "Reliability not initialized. Skipped.");
        }else
        {
            Log.LogInfo( "Processing reliability for Line Model method");
            CQualz solve2;
            std::shared_ptr<const CQualzResult> solve2Result = solve2.Compute( ctx.GetDataStore(), classifStore, redshiftRange, redshiftStep );

            if(solve2Result)
            {
                std::string predLabel="";
                bool retPredLabel = solve2Result->GetPredictedLabel( ctx.GetDataStore(), predLabel );
                solveResult->SetReliabilityLabel(predLabel);
            }
        }
    }
    //*/


    //finally save the linemodel results with (optionally) the zqual label
    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", solveResult );
    }else{
        return false;
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
    }else{
        return false;
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

    //retrieve the calibration dir path
    std::string calibrationDirPath;
    ctx.GetParameterStore().Get( "calibrationDir", calibrationDirPath );
    COperatorDTree7Solve Solve(calibrationDirPath);
    std::shared_ptr<CDTree7SolveResult> solveResult = Solve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), templateCategoryList, ctx.GetRayCatalog(),
                                                                        lambdaRange, redshiftRange, redshiftStep);

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", solveResult );
    }else{
        return false;
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

    //retrieve the calibration dir path
    std::string calibrationDirPath;
    ctx.GetParameterStore().Get( "calibrationDir", calibrationDirPath );
    COperatorDTreeASolve Solve(calibrationDirPath);
    std::shared_ptr<const CDTreeASolveResult> solveResult = Solve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), templateCategoryList, ctx.GetRayCatalog(),
                                                                        lambdaRange, redshiftRange, redshiftStep,
                                                                        correlationExtremumCount , overlapThreshold );

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", solveResult );
    }else{
        return false;
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

    //retrieve the calibration dir path
    std::string calibrationDirPath;
    ctx.GetParameterStore().Get( "calibrationDir", calibrationDirPath );
    COperatorDTreeBSolve Solve(calibrationDirPath);
    std::shared_ptr<const CDTreeBSolveResult> solveResult = Solve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), templateCategoryList, ctx.GetRayCatalog(),
                                                                        spcLambdaRange, redshifts);

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", solveResult );
    }else{
        return false;
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

    //retrieve the calibration dir path
    std::string calibrationDirPath;
    ctx.GetParameterStore().Get( "calibrationDir", calibrationDirPath );
    COperatorDTreeCSolve Solve(calibrationDirPath);
    std::shared_ptr<const CDTreeCSolveResult> solveResult = Solve.Compute( ctx.GetDataStore(), ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                                                                        ctx.GetTemplateCatalog(), templateCategoryList, ctx.GetRayCatalog(),
                                                                        spcLambdaRange, redshifts);

    if( solveResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", solveResult );
    }else{
        return false;
    }
    return true;
}


/**
 * @brief isPdfValid
 * @return
 */
Bool CProcessFlow::isPdfValid(CProcessFlowContext& ctx) const
{
    std::string scope_res = "zPDF/logposterior.logMargP_Z_data";
    auto results_pdf =  ctx.GetDataStore().GetGlobalResult( scope_res.c_str() );
    auto logzpdf1d = std::dynamic_pointer_cast<const CPdfMargZLogResult>( results_pdf.lock() );

    if(!logzpdf1d)
    {
        return false;
    }

    if(logzpdf1d->Redshifts.size()<2)
    {
        return false;
    }

    return true;
}
