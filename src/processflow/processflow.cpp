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


namespace NSEpic
{

class CBlindSolveTask
{
public:

    CProcessFlowContext&    m_Ctx;
    const CTemplate&        m_Tpl;
    const CSpectrum&        m_Spc;
    const CTemplate&        m_TplWithoutCont;
    const CSpectrum&        m_SpcWithoutCont;
    const TFloat64Range&    m_LambdaRanges;
    boost::mutex&           m_Mutex;

    CBlindSolveTask( CProcessFlowContext& ctx, const CTemplate& tpl, const CTemplate& tplWithoutCont, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const TFloat64Range& lambdaRange, boost::mutex& mutex ) :
        m_Ctx( ctx ),
        m_Tpl( tpl ),
        m_Spc( spc ),
        m_TplWithoutCont( tplWithoutCont ),
        m_SpcWithoutCont( spcWithoutCont ),
        m_LambdaRanges( lambdaRange ),
        m_Mutex( mutex )
    {

    }

    void operator()()
    {
        Bool retVal = true;

        CSpectrum s = m_Spc;
        s.GetSpectralAxis().ConvertToLogScale();

        // Create redshift initial list by spanning redshift acdross the given range, with the given delta
        CRedshifts redshifts( m_Ctx.GetRedshiftRange(), m_Ctx.GetRedshiftStep() );
        DebugAssert( redshifts.GetRedshiftsCount() > 0 );

        // Compute correlation factor at each of those redshifts
        COperatorCorrelation correlation;
        correlation.Compute( m_SpcWithoutCont, m_TplWithoutCont, m_LambdaRanges, redshifts, m_Ctx.GetOverlapThreshold() );

        // Find redshifts extremum
        Int32 nExtremums = 5;
        TPointList extremumList;
        CExtremum extremum( m_Ctx.GetRedshiftRange() , nExtremums);
        extremum.Find( redshifts.GetRedshifts(), correlation.GetResults().data(), redshifts.GetRedshiftsCount(), extremumList );

        if( extremumList.size() == 0 )
        {
            Log.LogInfo( "- Template: %s (LambdaRange: %f-%f:%f)", m_Tpl.GetName().c_str(), m_Tpl.GetLambdaRange().GetBegin(), m_Tpl.GetLambdaRange().GetEnd(), m_Tpl.GetResolution()  );
            Log.LogInfo( "|- No Redshift found" );
            return;
        }

        // Compute merit function
        CRedshifts newRedshifts( extremumList );
        COperatorChiSquare meritChiSquare;
        retVal = meritChiSquare.Compute( m_Spc, m_Tpl, m_LambdaRanges, newRedshifts, m_Ctx.GetOverlapThreshold() );
        if( !retVal )
        {
            Log.LogInfo( "Failed to compute chi square value");
            return;
        }

        // Store results
        {
            boost::lock_guard<boost::mutex> lock( m_Mutex );
            m_Ctx.AddResults( m_Tpl,
                              newRedshifts,
                              meritChiSquare.GetResults(), meritChiSquare.GetStatus(),
                              redshifts, correlation.GetResults() );
            Log.LogInfo( "- Template: %s (LambdaRange: %f-%f:%f)", m_Tpl.GetName().c_str(), m_Tpl.GetLambdaRange().GetBegin(), m_Tpl.GetLambdaRange().GetEnd(), m_Tpl.GetResolution() );

            Float64 merit = meritChiSquare.GetResults()[0];
            if( merit < 0.00001 )
                Log.LogInfo( "|- Redshift: %f Merit: %e", newRedshifts[0], merit );
            else
                Log.LogInfo( "|- Redshift: %f Merit: %f", newRedshifts[0], merit );

        }

    }

};

}

CProcessFlow::CProcessFlow( Int32 nbThread ) :
    m_ThreadPool( nbThread )
{

}

CProcessFlow::~CProcessFlow()
{

}

bool CProcessFlow::Process( CProcessFlowContext& ctx )
{
    return ProcessWithEL( ctx );

    const CProcessFlowContext::TTemplateCategoryList& tplCategoryList = ctx.GetTemplateCategoryList();
    if( std::find( tplCategoryList.begin(), tplCategoryList.end(), CTemplate::nCategory_Emission )==tplCategoryList.end() )
    {
        return ProcessWithoutEL( ctx );
    }
    else
    {
        return ProcessWithEL( ctx );
    }
}


bool CProcessFlow::ProcessWithoutEL( CProcessFlowContext& ctx )
{
    const CTemplateCatalog& templateCatalog = ctx.GetTemplateCatalog();
    const CProcessFlowContext::TTemplateCategoryList& templateCategotyList = ctx.GetTemplateCategoryList();

    m_ThreadPool.Reset();

    Log.LogInfo( "Process spectrum without EL (LambdaRange: %f-%f:%f)",
            ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());

    for( UInt32 i=0; i<templateCategotyList.size(); i++ )
    {
        Log.LogInfo( "Processing template category: %s", CTemplate::GetCategoryName( templateCategotyList[i] ) );
        Log.Indent();

        if( ctx.GetTemplateCategoryList()[i] == CTemplate::nCategory_Star )
        {
        }
        else
        {
            for( UInt32 j=0; j<templateCatalog.GetTemplateCount( (CTemplate::ECategory) templateCategotyList[i] ); j++ )
            {
                const CTemplate& tpl = templateCatalog.GetTemplate( (CTemplate::ECategory) templateCategotyList[i], j );
                const CTemplate& tplWithoutCont = templateCatalog.GetTemplateWithoutContinuum( (CTemplate::ECategory) templateCategotyList[i], j );
                const CSpectrum& spc = ctx.GetSpectrum();
                const CSpectrum& spcWithoutCont = ctx.GetSpectrumWithoutContinuum();
                const TFloat64Range& lambdaRange = ctx.GetLambdaRange();

                m_ThreadPool.AddTask( CBlindSolveTask( ctx, tpl, tplWithoutCont, spc, spcWithoutCont, lambdaRange, m_SyncMutex ) );
            }

        }

        m_ThreadPool.WaitForAllTaskToFinish();

        Log.UnIndent();
    }

    return true;
}


bool CProcessFlow::ProcessWithEL( CProcessFlowContext& ctx )
{
    //const CTemplateCatalog& templateCatalog = ctx.GetTemplateCatalog();
    //const CProcessFlowContext::TTemplateCategoryList& templateCategoryList = ctx.GetTemplateCategoryList();


    Log.LogInfo( "Process spectrum with EL (LambdaRange: %f-%f:%f)",
            ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());

    // --- EZ: EL Search
    bool retVal = ELSearch(ctx);
    if(ctx.GetDetectedRayCatalog().GetList().size()<2){
        return ProcessWithoutEL( ctx );
    }

    // --- EZ: EL Match
    CRayMatching rayMatching;
    retVal = rayMatching.Compute(ctx.GetDetectedRayCatalog(), ctx.GetRayCatalog(), ctx.GetRedshiftRange(), 2, 0.001 );
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
        retVal = ComputeMerits( ctx, selectedRedshift);
    }


    return true;
}

bool CProcessFlow::ELSearch( CProcessFlowContext& ctx )
{
    const CSpectrum& spc = ctx.GetSpectrum();
    const TFloat64Range& lambdaRange = ctx.GetLambdaRange();

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
    const TFloat64Range& lambdaRange = ctx.GetLambdaRange();

    const CTemplateCatalog& templateCatalog = ctx.GetTemplateCatalog();
    const CProcessFlowContext::TTemplateCategoryList& templateCategotyList = ctx.GetTemplateCategoryList();


    Log.LogInfo( "Process spectrum without EL (LambdaRange: %f-%f:%f)",
            ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution());

    for( UInt32 i=0; i<templateCategotyList.size(); i++ )
    {
        Log.LogInfo( "Processing template category: %s", CTemplate::GetCategoryName( templateCategotyList[i] ) );
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
                CRedshifts newRedshifts( redshifts );
                COperatorChiSquare meritChiSquare;
                Int32 retVal = meritChiSquare.Compute( spc, tpl, lambdaRange, newRedshifts, ctx.GetOverlapThreshold() );
                if( !retVal )
                {
                    Log.LogInfo( "Failed to compute chi square value");
                    return false;
                }

                // Store results
                {
                    ctx.AddMeritResults( tpl,
                                      newRedshifts,
                                      meritChiSquare.GetResults(), meritChiSquare.GetStatus(),
                                      redshifts );
                    Log.LogInfo( "- Template: %s (LambdaRange: %f-%f:%f)", tpl.GetName().c_str(), tpl.GetLambdaRange().GetBegin(), tpl.GetLambdaRange().GetEnd(), tpl.GetResolution() );

                    Float64 merit = meritChiSquare.GetResults()[0];
                    if( merit < 0.00001 )
                        Log.LogInfo( "|- Redshift: %f Merit: %e", newRedshifts[0], merit );
                    else
                        Log.LogInfo( "|- Redshift: %f Merit: %f", newRedshifts[0], merit );

                }
            }

        }


        Log.UnIndent();
    }

    return true;
}


