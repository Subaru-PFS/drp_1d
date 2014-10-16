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


using namespace __NS__;

IMPLEMENT_MANAGED_OBJECT( CProcessFlow )

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

namespace __NS__
{

class CBlindSolveTask
{
public:

    CProcessFlowContext&    m_Ctx;
    const CTemplate&        m_Tpl;
    const CSpectrum&        m_Spc;
    const TFloat64Range&    m_LambdaRanges;
    boost::mutex&           m_Mutex;

    CBlindSolveTask( CProcessFlowContext& ctx, const CTemplate& tpl, const CSpectrum& spc, const TFloat64Range& lambdaRange, boost::mutex& mutex ) :
        m_Ctx( ctx ),
        m_Tpl( tpl ),
        m_Spc( spc ),
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

        // Compute "Prepared" spctrum, i.e: spectrum with log scale and continuum removed
        CSpectrum preparedSpc = m_Spc;
        CTemplate preparedTpl = m_Tpl;

        preparedSpc.RemoveContinuum<CContinuumMedian>();
        preparedSpc.ConvertToLogScale();

        preparedTpl.RemoveContinuum<CContinuumMedian>();
        preparedTpl.ConvertToLogScale();

        // Compute correlation factor at each of those redshifts
        COperatorCorrelation correlation;
        correlation.Compute( preparedSpc, preparedTpl, m_LambdaRanges, redshifts, m_Ctx.GetOverlapThreshold() );

        // Find redshifts extremum
        TPointList extremumList;
        CExtremum extremum( m_Ctx.GetRedshiftRange(), m_Ctx.GetMaxCorrelationExtremumCount() );
        extremum.Find( redshifts.GetRedshifts(), correlation.GetResults().data(), redshifts.GetRedshiftsCount(), extremumList );


        // Compute fine grained correlation over each extremum
        TFloat64List newRedshiftsValues;
        for( UInt32 k=0; k<extremumList.size(); k++ )
        {
            SPoint p = extremumList[k];
            TFloat64Range newRange;

            if( p.X > redshifts.GetRange().GetEnd() || p.X < redshifts.GetRange().GetBegin() )
            {
                continue;
            }

            newRange.SetBegin( std::max(redshifts.GetRange().GetBegin(), p.X - m_Ctx.GetFineGrainedCorrelationRadius() ) );
            newRange.SetEnd( std::min(redshifts.GetRange().GetEnd(),     p.X + m_Ctx.GetFineGrainedCorrelationRadius() ) );

            CRedshifts   tmpRedshifts( newRange, m_Ctx.GetRedshiftStep()  );
            DebugAssert( tmpRedshifts.GetRedshiftsCount() > 0 );

            COperatorCorrelation tmpCorrelation;
            retVal = tmpCorrelation.Compute( preparedSpc, preparedTpl, m_LambdaRanges, tmpRedshifts, m_Ctx.GetOverlapThreshold()  );

            TPointList tmpExtremumList;
            CExtremum tmpExtremum;
            tmpExtremum.Find( tmpRedshifts.GetRedshifts(), tmpCorrelation.GetResults().data(), tmpRedshifts.GetRedshiftsCount(), tmpExtremumList );

            if( tmpExtremumList.size() )
            {
                newRedshiftsValues.push_back( tmpExtremumList[0].X );
            }
        }

        if( newRedshiftsValues.size( ) == 0 )
        {
            Log.LogInfo( "Template: %s (LambdaRange: %f-%f:%f:%f)", m_Tpl.GetName().c_str(), m_Tpl.GetLambdaRange().GetBegin(), m_Tpl.GetLambdaRange().GetEnd(), m_Tpl.GetResolution(), m_Tpl.GetLambdaRange().GetLength() );
            Log.LogInfo( "No Redshift found" );
            return;
        }

        preparedSpc = m_Spc;
        preparedTpl = m_Tpl;

        preparedSpc.ConvertToLinearScale();
        preparedTpl.ConvertToLinearScale();

        CRedshifts newRedshifts( newRedshiftsValues.data(), newRedshiftsValues.size() );
        COperatorChiSquare meritChiSquare;
        retVal = meritChiSquare.Compute( preparedSpc, preparedTpl, m_LambdaRanges, newRedshifts, m_Ctx.GetOverlapThreshold() );

        {
            boost::lock_guard<boost::mutex> lock( m_Mutex );
            m_Ctx.AddCorrelationResult( m_Tpl, newRedshifts, meritChiSquare.GetResults() );
            Log.LogInfo( "Template: %s (LambdaRange: %f-%f:%f:%f)", m_Tpl.GetName().c_str(), m_Tpl.GetLambdaRange().GetBegin(), m_Tpl.GetLambdaRange().GetEnd(), m_Tpl.GetResolution(), m_Tpl.GetLambdaRange().GetLength() );
            Log.LogInfo( "Redshift: %f Merit: %f", newRedshifts[0], meritChiSquare.GetResults()[0] );
        }

    }

};

}

bool CProcessFlow::ProcessWithoutEL( CProcessFlowContext& ctx )
{
    const CTemplateCatalog& templateCatalog = ctx.GetTemplateCatalog();
    const CProcessFlowContext::TTemplateCategoryList& templateCategotyList = ctx.GetTemplateCategoryList();

    m_ThreadPool.Reset();

    Log.LogInfo( "Process pectrum without EL (LambdaRange: %f-%f:%f:%f)",
            ctx.GetSpectrum().GetLambdaRange().GetBegin(), ctx.GetSpectrum().GetLambdaRange().GetEnd(), ctx.GetSpectrum().GetResolution(), ctx.GetSpectrum().GetLambdaRange().GetLength() );

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
                const CSpectrum& spc = ctx.GetSpectrum();
                const TFloat64Range& lambdaRange = ctx.GetLambdaRange();

                m_ThreadPool.AddTask( CBlindSolveTask( ctx, tpl, spc, lambdaRange, m_SyncMutex ) );

            }

        }

        m_ThreadPool.WaitForAllTaskToFinish();

        Log.UnIndent();
    }


    std::string tplName;
    Float64 redshift = 0.0;
    Float64 merit = 0.0;

    if( !ctx.GetBestCorrelationResult( redshift, merit, tplName ) )
        return false;

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
                const CSpectrum& spc = ctx.GetSpectrum();
                const TFloat64Range& lambdaRange = ctx.GetLambdaRange();

                m_ThreadPool.AddTask( CBlindSolveTask( ctx, tpl, spc, lambdaRange, m_SyncMutex ) );

            }

        }

        Log.UnIndent();
    }

    m_ThreadPool.WaitForAllTaskToFinish();

    std::string tplName;
    Float64 redshift = 0.0;
    Float64 merit = 0.0;

    if( !ctx.GetBestCorrelationResult( redshift, merit, tplName ) )
        return false;

    return true;
}
