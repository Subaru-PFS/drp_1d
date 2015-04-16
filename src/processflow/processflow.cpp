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
#include <epic/redshift/peak/detection.h>
#include <epic/redshift/gaussianfit/gaussianfit.h>
#include <epic/redshift/ray/ray.h>
#include <epic/redshift/ray/catalog.h>
#include <epic/redshift/ray/matching.h>
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
    return ProcessWithEL( ctx );

    const CProcessFlowContext::TTemplateCategoryList& tplCategoryList = ctx.GetParams().templateCategoryList;
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
    const CProcessFlowContext::TTemplateCategoryList& templateCategotyList = ctx.GetParams().templateCategoryList;

    for( UInt32 i=0; i<templateCategotyList.size(); i++ )
    {
        if( ctx.GetParams().templateCategoryList[i] == CTemplate::nCategory_Star )
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
    //const CTemplateCatalog& templateCatalog = ctx.GetTemplateCatalog();
    //const CProcessFlowContext::TTemplateCategoryList& templateCategoryList = ctx.GetTemplateCategoryList();


    Log.LogInfo( "Process spectrum with EL (LambdaRange: %f-%f:%f)",
            ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());

    // --- EZ: EL Search
    bool retVal = ELSearch(ctx);
    if(ctx.GetDetectedRayCatalog().GetList().size()<2){
        return false;
        //return ProcessWithoutEL( ctx );
    }

    // --- EZ: EL Match
    CRayMatching rayMatching;
    retVal = rayMatching.Compute(ctx.GetDetectedRayCatalog(), ctx.GetRayCatalog(), ctx.GetParams().redshiftRange, 2, 0.002 );
    // Store matching results
    {
        Float64 bestz=-1.0;
        Int32 bestzMatchingNum = -1;
        rayMatching.GetBestRedshift(bestz, bestzMatchingNum);
        ctx.SetRayMatchingResult(rayMatching.GetResults(), bestz, bestzMatchingNum);
    }
    return true;

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


    return true;
}

bool CProcessFlow::ELSearch( CProcessFlowContext& ctx )
{
    const CSpectrum& spc = ctx.GetSpectrum();
    const TFloat64Range& lambdaRange = ctx.GetParams().lambdaRange;

    // detect possible peaks
    Float64 winsize = 250.0;
    Float64 cut = 5;
    CPeakDetection detection;
    Float64 minsize = 3;
    Float64 maxsize = 90;
    bool retVal = detection.Compute( spc, lambdaRange, winsize, cut);
    const TInt32RangeList& resPeaks = detection.GetResults();
    UInt32 nPeaks = resPeaks.size();

    // filter the peaks with gaussian fit and create the detected rays catalog
    CRef<CRayCatalog> detectedRayCatalog = new CRayCatalog();
    for( UInt32 j=0; j<nPeaks; j++ )
    {
        bool toAdd = true;
        //find gaussian fit
        CGaussianFit fitter;
        CGaussianFit::EStatus status = fitter.Compute( spc, TInt32Range( resPeaks[j].GetBegin(), resPeaks[j].GetEnd() ) );
        if(status!=NSEpic::CGaussianFit::nStatus_Success){
            continue;
        }

        Float64 gaussAmp;
        Float64 gaussPos;
        Float64 gaussWidth;
        fitter.GetResults( gaussAmp, gaussPos, gaussWidth );

        /* check pos
        if(){
            toAdd = false;
        }
        //*/
        // check amp
        if(gaussAmp<0){
            toAdd = false;
        }
        // check width
        Float64 FWHM_FACTOR=2.35;
        if(gaussWidth<0){
            toAdd = false;
        }else{
            Float64 fwhm = FWHM_FACTOR*gaussWidth;
            if(fwhm<minsize){
                toAdd = false;
            }
            if(fwhm>maxsize){
                toAdd = false;
            }
        }

        if(toAdd){
            // check type weak or strong
            Int32 type = 2;
            // strong/weak test to do
            //...
            char buffer [64];
            sprintf(buffer,"detected_peak_%d",j);
            std::string peakName = buffer;
            detectedRayCatalog->Add( CRay( peakName, gaussPos, type ) );
        }
    }

    ctx.SetRayDetectionResult(*detectedRayCatalog);

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
    ctx.StorePerTemplateResult( tpl, "blindsolve.chisquare", *chisquareResult );

    return true;
}


