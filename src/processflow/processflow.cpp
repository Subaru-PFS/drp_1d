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
#include <epic/redshift/common/median.h>
#include <stdio.h>
#include <math.h>

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
    if(ctx.GetMethod()  == CProcessFlowContext::nMethod_LineMatching)
        return ProcessWithEL( ctx );

    if(ctx.GetMethod()  == CProcessFlowContext::nMethod_BlindSolve)
        return ProcessWithoutEL( ctx );

    /*// Method auto selection
    const CProcessFlowContext::TTemplateCategoryList& tplCategoryList = ctx.GetTemplateCategoryList();
    if( std::find( tplCategoryList.begin(), tplCategoryList.end(), CTemplate::nCategory_Emission )==tplCategoryList.end() )
    {
        return ProcessWithoutEL( ctx );
    }
    else
    {
        return ProcessWithEL( ctx );
    }
    //*/

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
    Log.LogInfo( "Process spectrum with EL (LambdaRange: %f-%f:%f)",
            ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());

    // --- EZ: EL Search
    bool retVal = ELSearch(ctx);
    if(ctx.GetDetectedRayCatalog().GetList().size()<1){
        return false;
        //return ProcessWithoutEL( ctx );
    }

    // --- EZ: EL Match
    CRayMatching rayMatching;
    retVal = rayMatching.Compute(ctx.GetDetectedRayCatalog(), ctx.GetRayCatalog(), ctx.GetRedshiftRange(), 1, 0.002 );
    // Store matching results
    {
        Float64 bestz=-1.0;
        Int32 bestzMatchingNum = -1;
        rayMatching.GetBestRedshift(bestz, bestzMatchingNum);
        std::string strDesc;
        rayMatching.GetDescription(strDesc);
        ctx.SetRayMatchingResult(rayMatching.GetResults(), bestz, bestzMatchingNum, strDesc);
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
    const TFloat64Range& lambdaRange = ctx.GetLambdaRange();
    const CSpectrumFluxAxis fluxAxis = spc.GetFluxAxis();

    // detect possible peaks
    Float64 winsize = 250.0;
    Float64 cut = 5.0;
    Float64 strongcut = 2.0;
    CPeakDetection detection;
    Float64 minsize = 3;
    Float64 maxsize = 70;
    bool retVal = detection.Compute( spc, lambdaRange, winsize, cut);
    const TInt32RangeList& resPeaks = detection.GetResults();
    const TInt32RangeList& resPeaksEnlarged = detection.GetResultsEnlarged();

    UInt32 nPeaks = resPeaks.size();

    // filter the peaks with gaussian fit and create the detected rays catalog
    CRef<CRayCatalog> detectedRayCatalog = new CRayCatalog();
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
            // use noise spectrum
            const Float64* error = fluxAxis.GetError();
            if( error!= NULL ){
                // check if noise file has been loaded
                bool isNoiseOnes = true;
                for ( Int32 i=left; i<right; i++)
                {
                    if(error[i]!=1.0){
                        isNoiseOnes = false;
                        break;
                    }
                }

                if(!isNoiseOnes){
                    Float64 mean_noise = 0.0;
                    Int32 n_mean_noise = 0;
                    for ( Int32 i=left; i<right; i++)
                    {
                        mean_noise += error[i];
                        n_mean_noise ++;
                    }
                    if(n_mean_noise>0){
                        mean_noise /= n_mean_noise;
                    }
                    // choose between noise mean or xmad
                    if(mean_noise>xmad){
                        xmad = mean_noise;
                    }
                }

            }
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
            detectedRayCatalog->Add( CRay( peakName, gaussPos, 2, force ) );
        }
    }

    ctx.SetRayDetectionResult(*detectedRayCatalog);

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
    const TFloat64Range& lambdaRange = ctx.GetLambdaRange();

    const CTemplateCatalog& templateCatalog = ctx.GetTemplateCatalog();
    const CProcessFlowContext::TTemplateCategoryList& templateCategotyList = ctx.GetTemplateCategoryList();


    Log.LogInfo( "Process spectrum for Merit (LambdaRange: %f-%f:%f)",
            ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());

    for( UInt32 i=0; i<templateCategotyList.size(); i++ )
    {
        Log.LogInfo( "Processing merits for template category: %s", CTemplate::GetCategoryName( templateCategotyList[i] ) );
        Log.Indent();

        if( ctx.GetTemplateCategoryList()[i] == CTemplate::nCategory_Star )
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
                Int32 retVal = meritChiSquare.Compute( spc, tpl, lambdaRange, redshifts, ctx.GetOverlapThreshold() );
                if( !retVal )
                {
                    Log.LogInfo( "Failed to compute chi square value");
                    return false;
                }

                // Store results
                {
                    ctx.AddMeritResults( tpl,
                                      redshifts,
                                      meritChiSquare.GetResults(), meritChiSquare.GetStatus(),
                                      redshifts );
                    Log.LogInfo( "- Template: %s (LambdaRange: %f-%f:%f)", tpl.GetName().c_str(), tpl.GetLambdaRange().GetBegin(), tpl.GetLambdaRange().GetEnd(), tpl.GetResolution() );

                    Float64 merit = meritChiSquare.GetResults()[0];
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


