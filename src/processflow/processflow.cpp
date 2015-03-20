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
        TPointList extremumList;
        CExtremum extremum( m_Ctx.GetRedshiftRange() );
        extremum.Find( redshifts.GetRedshifts(), correlation.GetResults().data(), redshifts.GetRedshiftsCount(), extremumList );

        if( extremumList.size() == 0 )
        {
            Log.LogInfo( "- Template: %s (LambdaRange: %f-%f:%f)", m_Tpl.GetName().c_str(), m_Tpl.GetLambdaRange().GetBegin(), m_Tpl.GetLambdaRange().GetEnd(), m_Tpl.GetResolution()  );
            Log.LogInfo( "|- No Redshift found" );
            return;
        }

        // Compute merit function
        CRedshifts newRedshifts( extremumList[0].X );
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
            m_Ctx.AddResults( m_Tpl, newRedshifts, meritChiSquare.GetResults(), redshifts, correlation.GetResults() );
            Log.LogInfo( "- Template: %s (LambdaRange: %f-%f:%f)", m_Tpl.GetName().c_str(), m_Tpl.GetLambdaRange().GetBegin(), m_Tpl.GetLambdaRange().GetEnd(), m_Tpl.GetResolution() );
            Log.LogInfo( "|- Redshift: %f Merit: %f", newRedshifts[0], meritChiSquare.GetResults()[0] );
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
    const CTemplateCatalog& templateCatalog = ctx.GetTemplateCatalog();
    const CProcessFlowContext::TTemplateCategoryList& templateCategotyList = ctx.GetTemplateCategoryList();

    m_ThreadPool.Reset();

    for( UInt32 i=0; i<templateCategotyList.size(); i++ )
    {
        Log.LogInfo( "Processing template category: %s", CTemplate::GetCategoryName( templateCategotyList[i] ) );
        Log.Indent();

        if( ctx.GetTemplateCategoryList()[i] == CTemplate::nCategory_Star )
        {
            Log.LogError( "Processing not yet implemented for Star template" );
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

        Log.UnIndent();
    }

    m_ThreadPool.WaitForAllTaskToFinish();


    return true;
}
