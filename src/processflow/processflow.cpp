#include <epic/redshift/processflow/processflow.h>

#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/continuum/standard.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/chicorr.h>
#include <epic/redshift/operator/blindsolveresult.h>
#include <epic/redshift/operator/raydetectionresult.h>
#include <epic/redshift/operator/raydetection.h>
#include <epic/redshift/operator/blindsolve.h>
#include <epic/redshift/operator/blindsolveresult.h>
#include <epic/redshift/operator/blindsolve.h>
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

#include <epic/redshift/operator/peakdetectionresult.h>
#include <epic/redshift/operator/raydetectionresult.h>
#include <epic/redshift/operator/raymatchingresult.h>

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
        return ProcessLineMatching( ctx );

    if(ctx.GetParams().method  == CProcessFlowContext::nMethod_BlindSolve || ctx.GetParams().method  == CProcessFlowContext::nMethod_FullSolve)
        return ProcessWithoutEL( ctx );

    if(ctx.GetParams().method  == CProcessFlowContext::nMethod_DecisionalTree7)
        return ProcessDecisionalTree7( ctx );


    return false;
}


Bool CProcessFlow::ProcessWithoutEL( CProcessFlowContext& ctx )
{



    // Remove Star category

    if(ctx.GetParams().method  == CProcessFlowContext::nMethod_BlindSolve ){

        COperatorBlindSolve blindSolve;
        CConstRef<CBlindSolveResult> blindsolveResult = blindSolve.Compute( ctx, ctx.GetSpectrum(), ctx.GetSpectrumWithoutContinuum(),
                            ctx.GetTemplateCatalog(), ctx.GetParams().templateCategoryList,
                            ctx.GetParams().lambdaRange, ctx.GetParams().redshiftRange, ctx.GetParams().redshiftStep,
                            ctx.GetParams().correlationExtremumCount, ctx.GetParams().overlapThreshold );

        if( blindsolveResult ) {
            ctx.StoreGlobalResult( "blindsolve", *blindsolveResult );
        }


    }else if(ctx.GetParams().method  == CProcessFlowContext::nMethod_FullSolve){

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

                    FullSolve( ctx, tpl, tplWithoutCont );

                }
            }
        }

    }



    return true;
}

Bool CProcessFlow::ProcessDecisionalTree7( CProcessFlowContext& ctx )
{
    Log.LogInfo( "Process Decisional Tree 7(LambdaRange: %f-%f:%f)",
            ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());


    // Peak Detection
    Float64 winsize = 250.0;
    Float64 cut = 5.0;
    CPeakDetection peakDetection;
    CConstRef<CPeakDetectionResult> peakDetectionResult = peakDetection.Compute( ctx.GetSpectrum(), ctx.GetSpectrum().GetLambdaRange(), winsize, cut);
    ctx.StoreGlobalResult( "peakdetection", *peakDetectionResult );
    Log.LogInfo( "Peak Detection output: %d peaks found", peakDetectionResult->PeakList.size());

    // check Peak Detection results
    if(peakDetectionResult->PeakList.size()<1){
        Log.LogInfo( "No Peak found, switching to ProcessWithoutEL");
        return ProcessWithoutEL( ctx );
    }

    // Ray Detection
    CRayDetection rayDetection;
    CConstRef<CRayDetectionResult> rayDetectionResult = rayDetection.Compute( ctx.GetSpectrum(), ctx.GetSpectrum().GetLambdaRange(), peakDetectionResult->PeakList, peakDetectionResult->EnlargedPeakList );
    ctx.StoreGlobalResult( "raycatalog", *rayDetectionResult );
    Log.LogInfo( "Ray Detection output: %d ray(s) found", rayDetectionResult->RayCatalog.GetList().size());

    // check Ray Detection results
    Int32 nRaysDetected = rayDetectionResult->RayCatalog.GetList().size();
    if( nRaysDetected < 1){
        Log.LogInfo( "Not ray found, switching to ProcessWithoutEL");
        return ProcessWithoutEL( ctx );
    }

    // Ray Match
    CRayMatching rayMatching;
    Int32 MinMatchNum = 1;
    Float64 tol = 0.002;
    CRef<CRayMatchingResult> rayMatchingResult = rayMatching.Compute(rayDetectionResult->RayCatalog, ctx.GetRayCatalog(), ctx.GetParams().redshiftRange, MinMatchNum, tol );
    // Store matching results
    ctx.StoreGlobalResult( "raymatching", *rayMatchingResult );

    //check ray matching results
    if(rayMatchingResult->GetSolutionsListOverNumber(0).size()<1){
        Log.LogInfo( "Not match found, switching to ProcessWithoutEL");
        return ProcessWithoutEL( ctx );
    }

    Int32 matchNum = rayMatchingResult->GetMaxMatchingNumber();
    UInt32 nStrongPeaks = rayDetectionResult->RayCatalog.GetFilteredList(CRay::nType_Emission, CRay::nForce_Strong).size();

    // match num >= 3
    if(matchNum >= 3 || (nStrongPeaks < 1 && matchNum >= 2)){
        if(matchNum >= 3){
            Log.LogInfo( "match num >= 3");
        }
        if(nStrongPeaks < 1 && matchNum >= 2){
            Log.LogInfo( "match num <1, MatchNum>=2");
        }
        Log.LogInfo( "compute merits on redshift candidates from ray matching");
        TFloat64List selectedRedshift;
        CRayMatchingResult::TSolutionSetList selectedResults = rayMatchingResult->GetSolutionsListOverNumber(matchNum-1);
        for( UInt32 j=0; j<selectedResults.size(); j++ )
        {
            Float64 z = rayMatchingResult->GetMeanRedshiftSolution(selectedResults[j]);
            selectedRedshift.push_back(z);
        }
        return ComputeMerits( ctx, selectedRedshift);
    }

//    // 3 lines 1 match or 4 lines 2 matches
//    if(nRaysDetected - matchNum >= 2){
//        Log.LogInfo( "3 lines 1 match or 4 lines 2 matches...");
//    }

    if(nStrongPeaks > 0){
        Log.LogInfo( "Ray Matching with %d strong peaks", nStrongPeaks);
        CRayMatching rayMatchingStrong;
        Int32 MinMatchNum = 1;
        Float64 tol = 0.002;
        CRef<CRayMatchingResult> rayMatchingStrongResult = rayMatchingStrong.Compute(rayDetectionResult->RayCatalog, ctx.GetRayCatalog(), ctx.GetParams().redshiftRange, MinMatchNum, tol, CRay::nType_Emission, CRay::nForce_Strong );
        Int32 matchNumStrong = rayMatchingStrongResult->GetMaxMatchingNumber();

        if(matchNumStrong>1){
            Log.LogInfo( "match num strong >= 2, compute merits on redshift candidates from strong ray matching");
            TFloat64List selectedRedshift;
            CRayMatchingResult::TSolutionSetList selectedResults = rayMatchingStrongResult->GetSolutionsListOverNumber(matchNumStrong-1);
            for( UInt32 j=0; j<selectedResults.size(); j++ )
            {
                Float64 z = rayMatchingStrongResult->GetMeanRedshiftSolution(selectedResults[j]);
                selectedRedshift.push_back(z);
            }
            return ComputeMerits( ctx, selectedRedshift);
        }else{
            Log.LogInfo( "Not match found with strong lines, switching to ProcessWithoutEL");
            return ProcessWithoutEL( ctx );
        }
    }

    Log.LogInfo( "Switching to ProcessWithoutEL");
    return ProcessWithoutEL( ctx );
}

Bool CProcessFlow::ProcessLineMatching( CProcessFlowContext& ctx )
{
    Log.LogInfo( "Process spectrum with EL (LambdaRange: %f-%f:%f)",
            ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());

    // --- EZ: EL Search
    // detect possible peaks
    Float64 winsize = 250.0;
    Float64 cut = 5.0;

    CPeakDetection peakDetection;
    CConstRef<CPeakDetectionResult> peakDetectionResult = peakDetection.Compute( ctx.GetSpectrum(), ctx.GetSpectrum().GetLambdaRange(), winsize, cut);
    if( peakDetectionResult )
        ctx.StoreGlobalResult( "peakdetection", *peakDetectionResult );

    CRayDetection rayDetection;
    CConstRef<CRayDetectionResult> rayDetectionResult = rayDetection.Compute( ctx.GetSpectrum(), ctx.GetSpectrum().GetLambdaRange(), peakDetectionResult->PeakList, peakDetectionResult->EnlargedPeakList );

    if( rayDetectionResult ) {
        ctx.StoreGlobalResult( "raycatalog", *rayDetectionResult );

        if(rayDetectionResult->RayCatalog.GetList().size()<1){
            return false;
        }
    }

    // --- EZ: EL Match
    CRayMatching rayMatching;
    Int32 MinMatchNum = 1;
    Float64 tol = 0.002;
    CRef<CRayMatchingResult> rayMatchingResult = rayMatching.Compute(rayDetectionResult->RayCatalog, ctx.GetRayCatalog(), ctx.GetParams().redshiftRange, MinMatchNum, tol );

    // Store matching results
    if( rayMatchingResult )
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
    const TTemplateCategoryList& templateCategotyList = ctx.GetParams().templateCategoryList;


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

Bool CProcessFlow::FullSolve( CProcessFlowContext& ctx, const CTemplate& tpl, const CTemplate& tplWithoutCont )
{
    const CSpectrum& spc = ctx.GetSpectrum();
    const CSpectrum& spcWithoutCont = ctx.GetSpectrumWithoutContinuum();

    CSpectrum s = spc;
    s.GetSpectralAxis().ConvertToLogScale();

    // Create redshift initial list by spanning redshift acdross the given range, with the given delta
    TFloat64List redshifts = ctx.GetParams().redshiftRange.SpreadOver( ctx.GetParams().redshiftStep );
    DebugAssert( redshifts.size() > 0 );

    // Compute correlation factor at each of those redshifts
    COperatorChicorr chicorr;
    CChisquareResult* chisquareResult = new CChisquareResult();
    CCorrelationResult* correlationResult = new CCorrelationResult();
    chicorr.Compute( spc, spcWithoutCont, tpl, tplWithoutCont, ctx.GetParams().lambdaRange, redshifts, ctx.GetParams().overlapThreshold, correlationResult, chisquareResult);


    // Store results
    ctx.StorePerTemplateResult( tpl, "blindsolve.correlation", *correlationResult );
    ctx.StorePerTemplateResult( tpl, "blindsolve.merit", *chisquareResult );

    return true;
}

Bool CProcessFlow::FullSolveBrute( CProcessFlowContext& ctx, const CTemplate& tpl, const CTemplate& tplWithoutCont )
{
    const CSpectrum& spc = ctx.GetSpectrum();
    const CSpectrum& spcWithoutCont = ctx.GetSpectrumWithoutContinuum();

    CSpectrum s = spc;
    s.GetSpectralAxis().ConvertToLogScale();

    // Create redshift initial list by spanning redshift acdross the given range, with the given delta
    TFloat64List redshifts = ctx.GetParams().redshiftRange.SpreadOver( ctx.GetParams().redshiftStep );
    DebugAssert( redshifts.size() > 0 );

    // Compute correlation factor at each of those redshifts
    // Compute correlation factor at each of those redshifts
    COperatorCorrelation correlation;
    CRef<CCorrelationResult> correlationResult = (CCorrelationResult*) correlation.Compute( spcWithoutCont, tplWithoutCont, ctx.GetParams().lambdaRange, redshifts, ctx.GetParams().overlapThreshold );

    COperatorChiSquare meritChiSquare;
    CRef<CCorrelationResult> chisquareResult = (CCorrelationResult*)meritChiSquare.Compute( spc, tpl, ctx.GetParams().lambdaRange, redshifts, ctx.GetParams().overlapThreshold );
    if( !chisquareResult )
    {
        Log.LogInfo( "Failed to compute chi square value");
        return false;
    }

    //COperatorChicorr chicorr;
    //CChisquareResult* chisquareResult = new CChisquareResult();
    //CCorrelationResult* correlationResult = new CCorrelationResult();
    //chicorr.Compute( spc, spcWithoutCont, tpl, tplWithoutCont, ctx.GetParams().lambdaRange, redshifts, ctx.GetParams().overlapThreshold, correlationResult, chisquareResult);


    // Store results
    ctx.StorePerTemplateResult( tpl, "blindsolve.correlation", *correlationResult );
    ctx.StorePerTemplateResult( tpl, "blindsolve.merit", *chisquareResult );

    return true;
}
