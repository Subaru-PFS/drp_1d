#include <epic/redshift/processflow/processflow.h>

#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/continuum/standard.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/blindsolveresult.h>
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

#include <epic/redshift/operator/raydetectionresult.h>
#include <epic/redshift/operator/raymatchingresult.h>

Bool CProcessFlow::ProcessWithEL( CProcessFlowContext& ctx )
{
    Log.LogInfo( "Process spectrum with EL (LambdaRange: %f-%f:%f)",
            ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());

    // --- EZ: EL Search
    bool retVal = ELSearch(ctx);

    const CRayDetectionResult* rayDetectionResult = (const CRayDetectionResult*)ctx.GetGlobalResult( "raycatalog" );
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

#include <epic/redshift/operator/raydetectionresult.h>

bool CProcessFlow::ELSearch( CProcessFlowContext& ctx )
{
    const CSpectrum& spc = ctx.GetSpectrum();
    const CSpectrumFluxAxis fluxAxis = spc.GetFluxAxis();
    const TFloat64Range& lambdaRange = ctx.GetParams().lambdaRange;


    // detect possible peaks
    Float64 winsize = 250.0;
    Float64 cut = 5.0;
    Float64 strongcut = 2.0;
    CPeakDetection detection;
    Float64 minsize = 3;
    Float64 maxsize = 90;
    bool retVal = detection.Compute( spc, lambdaRange, winsize, cut);
    const TInt32RangeList& resPeaks = detection.GetResults();
    const TInt32RangeList& resPeaksEnlarged = detection.GetResultsEnlarged();

    UInt32 nPeaks = resPeaks.size();

    CRef<CRayDetectionResult> result = new CRayDetectionResult();

    // filter the peaks with gaussian fit and create the detected rays catalog
    for( UInt32 j=0; j<nPeaks; j++ )
    {
        bool toAdd = true;
        //find gaussian fit
        CGaussianFit fitter;
        CGaussianFit::EStatus status = fitter.Compute( spc, TInt32Range( resPeaksEnlarged[j].GetBegin(), resPeaksEnlarged[j].GetEnd() ) );
        if(status!=NSEpic::CGaussianFit::nStatus_Success){
            continue;
        }

        Float64 gaussAmp;
        Float64 gaussPos;
        Float64 gaussWidth;
        fitter.GetResults( gaussAmp, gaussPos, gaussWidth );
        Float64 gaussCont;
        fitter.GetResultsPolyCoeff0( gaussCont );

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
        if(toAdd){
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
        }

        //find max value and pos
        Float64 max_value = DBL_MIN;
        Int32 max_index = -1;
        for( Int32 k=resPeaks[j].GetBegin(); k<resPeaks[j].GetEnd()+1; k++ )
        {
            if(max_value < fluxAxis[k]){
                max_value = fluxAxis[k];
                max_index = k;
            }
        }

        // Check if gaussian fit is very different from peak itself
        if(toAdd){
            // check max_gauss vs max_spectrum
            Float64 gaussAmp_with_cont = gaussAmp + gaussCont;
            if(gaussAmp_with_cont/max_value <= 0.65 || gaussAmp_with_cont/max_value >= 1.35){
                toAdd = false;
            }
            if(abs(gaussPos-spc.GetSpectralAxis()[max_index])>3.*spc.GetResolution()){
                toAdd = false;
            }

        }

        // check type weak or strong
        Int32 force = 1; //weak by default
        if(toAdd){
            // strong/weak test to do
            const Float64* fluxData = fluxAxis.GetSamples();
            Int32 windowSampleCount = winsize / spc.GetResolution();
            CMedian<Float64> medianProcessor;
            int left = max(0, (Int32)(max_index-windowSampleCount/2.0+0.5) ) ;
            int right = min((Int32)fluxAxis.GetSamplesCount()-1, (Int32)(max_index + windowSampleCount/2.0) )+1;
            int size_odd = right - left +1;
            if( (int)(size_odd/2) == size_odd/2 ){
                size_odd -= 1;
            }

            Float64 med = medianProcessor.Find( fluxData + left, size_odd);
            Float64 xmad = XMadFind( fluxData + left, size_odd , med);
            // TODO: check with the noise spectrum
            /*if noise!=None:
            noise_mean=noise[left_index:right_index].mean()
            if noise_mean>xmadm:
                # Use noise file
                xmadm=noise_mean
                self.msg.debug(self.__module__, '18.96', 0, self.spectrum.name, "%4.1f" % pos)
            */
            Float64 max_value_no_continuum = max_value - med;
            Float64 ratioAmp=max_value_no_continuum/xmad;

            if(ratioAmp<cut){
                toAdd = false;
            }else if(ratioAmp>cut*strongcut){
                force = 2; //strong
            }
        }

        if(toAdd){
            char buffer [64];
            sprintf(buffer,"detected_peak_%d",j);
            std::string peakName = buffer;
            result->RayCatalog.Add( CRay( peakName, gaussPos, CRay::nType_Emission, force ) );

        }
    }

    ctx.StoreGlobalResult( "raycatalog", *result );

}

Float64 CProcessFlow::XMadFind( const Float64* x, Int32 n, Float64 median )
{
    std::vector<Float64> xdata;
    Float64 xmadm = 0.0;

    xdata.reserve( n );

    for( Int32 i=0;i<n; i++ )
    {
        xdata[i] = fabs( x[i]-median );
    }

    CQuickSort<Float64> sort;

    sort.Sort( xdata.data(), n);

    if( ((float)n)/2. - int(n/2.) == 0 )
    {
        UInt32 i1 = n/2 -1;
        UInt32 i2 = n/2;
        xmadm = 0.5*(xdata[i1]+xdata[i2]);
    }
    else
    {
        UInt32 i1 = int(n/2);
        xmadm = xdata[i1];
    }

    return xmadm;
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


