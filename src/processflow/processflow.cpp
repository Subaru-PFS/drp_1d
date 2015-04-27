#include <epic/redshift/processflow/processflow.h>

#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/continuum/standard.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/blindsolveresult.h>
#include <epic/redshift/operator/raydetectionresult.h>
#include <epic/redshift/operator/raydetection.h>
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

    if(ctx.GetParams().method  == CProcessFlowContext::nMethod_LineMatching)
        return ProcessWithEL( ctx );

    if(ctx.GetParams().method  == CProcessFlowContext::nMethod_BlindSolve)
        return ProcessWithoutEL( ctx );


    return false;
}


Bool CProcessFlow::ProcessWithoutEL( CProcessFlowContext& ctx )
{
    const CTemplateCatalog& templateCatalog = ctx.GetTemplateCatalog();

    for( UInt32 i=0; i<ctx.GetParams().templateCategoryList.size(); i++ )
    {
        CTemplate::ECategory category = ctx.GetParams().templateCategoryList[i];
        if( category == CTemplate::nCategory_Star )
        {
        }
        else
        {
            for( UInt32 j=0; j<templateCatalog.GetTemplateCount( category ); j++ )
            {
                const CTemplate& tpl = templateCatalog.GetTemplate( category, j );
                const CTemplate& tplWithoutCont = templateCatalog.GetTemplateWithoutContinuum( category, j );


                BlindSolve( ctx, tpl, tplWithoutCont );
            }
        }
    }

    return true;
}

#include <epic/redshift/operator/peakdetectionresult.h>
#include <epic/redshift/operator/raydetectionresult.h>
#include <epic/redshift/operator/raymatchingresult.h>

Bool CProcessFlow::ProcessWithEL( CProcessFlowContext& ctx )
{
    Log.LogInfo( "Process spectrum with EL (LambdaRange: %f-%f:%f)",
            ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());

    // --- EZ: EL Search
    // detect possible peaks
    Float64 winsize = 250.0;
    Float64 cut = 5.0;

    CPeakDetection peakDetection;
    CConstRef<CPeakDetectionResult> peakDetectionResult = peakDetection.Compute( ctx.GetSpectrum(), ctx.GetSpectrum().GetLambdaRange(), winsize, cut);
    ctx.StoreGlobalResult( "peakdetection", *peakDetectionResult );

    CRayDetection rayDetection;
    CConstRef<CRayDetectionResult> rayDetectionResult = rayDetection.Compute( ctx.GetSpectrum(), ctx.GetSpectrum().GetLambdaRange(), peakDetectionResult->PeakList, peakDetectionResult->EnlargedPeakList );
    ctx.StoreGlobalResult( "raycatalog", *rayDetectionResult );

    if(rayDetectionResult->RayCatalog.GetList().size()<2){
        return false;
        //return ProcessWithoutEL( ctx );
    }

    // --- EZ: EL Match
    CRayMatching rayMatching;
    CRef<CRayMatchingResult> rayMatchingResult = rayMatching.Compute(rayDetectionResult->RayCatalog, ctx.GetRayCatalog(), ctx.GetParams().redshiftRange, 2, 0.002 );

    // Store matching results
    ctx.StoreGlobalResult( "raymatching", *rayMatchingResult );


    return true;
/*
    Int32 maxMatchingNum = rayMatching.GetMaxMatchingNumber();
    if(maxMatchingNum>2){ //ez equivalent to SolveDecisionalTree2:three_lines_match()
        TFloat64List selectedRedshift;
        TRedshiftSolutionSetList selectedResults = rayMatching.GetSolutionsListOverNumber(2);
        for( UInt32 j=0; j<selectedResults.size(); j++ )
        {
            Float64 z = rayMatching.GetMeanRedshiftSolution(selectedResults[j]);
            selectedRedshift.push_back(z);
        }

        // --- EZ: solve_basic_nocorrelation
        //retVal = ComputeMerits( ctx, selectedRedshift);
    }

*/
    return true;
}



bool CProcessFlow::ComputeMerits( CProcessFlowContext& ctx, const TFloat64List& redshifts)
{
    const CSpectrum& spc = ctx.GetSpectrum();
    const TFloat64Range& lambdaRange = ctx.GetParams().lambdaRange;

    const CTemplateCatalog& templateCatalog = ctx.GetTemplateCatalog();
    const CProcessFlowContext::TTemplateCategoryList& templateCategotyList = ctx.GetParams().templateCategoryList;


    Log.LogInfo( "Process spectrum for Merit (LambdaRange: %f-%f:%f)",
            ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());

    for( UInt32 i=0; i<templateCategotyList.size(); i++ )
    {
        Log.LogInfo( "Processing merits for template category: %s", CTemplate::GetCategoryName( templateCategotyList[i] ) );
        Log.Indent();

        if( ctx.GetParams().templateCategoryList[i] == CTemplate::nCategory_Star )
        {
        }
        else
        {
            for( UInt32 j=0; j<templateCatalog.GetTemplateCount( (CTemplate::ECategory) templateCategotyList[i] ); j++ )
            {
                const CTemplate& tpl = templateCatalog.GetTemplate( (CTemplate::ECategory) templateCategotyList[i], j );

                //const CTemplate& tplWithoutCont = templateCatalog.GetTemplateWithoutContinuum( (CTemplate::ECategory) templateCategotyList[i], j );


                // Compute merit function
                COperatorChiSquare meritChiSquare;
                CRef<CChisquareResult>  chisquareResults = (CChisquareResult*)meritChiSquare.Compute( spc, tpl, lambdaRange, redshifts, ctx.GetParams().overlapThreshold );
                if( !chisquareResults )
                {
                    Log.LogInfo( "Failed to compute chi square value");
                    return false;
                }

                // Store results
                {
                    ctx.StorePerTemplateResult( tpl, "merit.chisquare", *chisquareResults );

                    Log.LogInfo( "- Template: %s (LambdaRange: %f-%f:%f)", tpl.GetName().c_str(), tpl.GetLambdaRange().GetBegin(), tpl.GetLambdaRange().GetEnd(), tpl.GetResolution() );

                    Float64 merit = chisquareResults->ChiSquare[0];
                    if( merit < 0.00001 )
                        Log.LogInfo( "|- Redshift: %f Merit: %e", redshifts[0], merit );
                    else
                        Log.LogInfo( "|- Redshift: %f Merit: %f", redshifts[0], merit );

                }
            }
        }
    }

    return true;
}


Bool CProcessFlow::BlindSolve( CProcessFlowContext& ctx, const CTemplate& tpl, const CTemplate& tplWithoutCont )
{
    const CSpectrum& spc = ctx.GetSpectrum();
    const CSpectrum& spcWithoutCont = ctx.GetSpectrumWithoutContinuum();

    CSpectrum s = spc;
    s.GetSpectralAxis().ConvertToLogScale();

    // Create redshift initial list by spanning redshift acdross the given range, with the given delta
    TFloat64List redshifts = ctx.GetParams().redshiftRange.SpreadOver( ctx.GetParams().redshiftStep );
    DebugAssert( redshifts.size() > 0 );

    // Compute correlation factor at each of those redshifts
    COperatorCorrelation correlation;
    CRef<CCorrelationResult> correlationResult = (CCorrelationResult*) correlation.Compute( spcWithoutCont, tplWithoutCont, ctx.GetParams().lambdaRange, redshifts, ctx.GetParams().overlapThreshold );

    // Find redshifts extremum
    Int32 nExtremums = 5;
    TPointList extremumList;
    CExtremum extremum( ctx.GetParams().redshiftRange , nExtremums);
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
    CRef<CCorrelationResult> chisquareResult = (CCorrelationResult*)meritChiSquare.Compute( spc, tpl, ctx.GetParams().lambdaRange, extremumRedshifts, ctx.GetParams().overlapThreshold );
    if( !chisquareResult )
    {
        Log.LogInfo( "Failed to compute chi square value");
        return false;
    }

    // Store results
    ctx.StorePerTemplateResult( tpl, "blindsolve.correlation", *correlationResult );
    ctx.StorePerTemplateResult( tpl, "blindsolve.merit", *chisquareResult );

    CRef<CBlindSolveResult>  chisquareResults = new CBlindSolveResult();
    ctx.StoreGlobalResult( "blindsolve", *chisquareResults );
    return true;
}


