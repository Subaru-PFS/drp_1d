#include <epic/redshift/processflow/context.h>

#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/io/genericreader.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/noise/flat.h>
#include <epic/redshift/noise/fromfile.h>
#include <epic/redshift/ray/catalog.h>
#include <epic/redshift/continuum/median.h>
#include <epic/core/log/log.h>

#include <float.h>

#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CProcessFlowContext )

CProcessFlowContext::SParam::SParam()
{
    redshiftRange = TFloat64Range( 0.0, 3.0 );
    redshiftStep = 0.0001;
    lambdaRange = TFloat64Range( 3800.0, 12500.0 );
    smoothWidth = 3;
    overlapThreshold = 1.0;

    //templateCategoryList.push_back( CTemplate::nCategory_Emission );
    templateCategoryList.push_back( CTemplate::nCategory_Galaxy );
    templateCategoryList.push_back( CTemplate::nCategory_Star );
    templateCategoryList.push_back( CTemplate::nCategory_Qso );
}

CProcessFlowContext::CProcessFlowContext() :
    m_FineGrainedCorrelationRadius( 0.001 ),
    m_RedshiftRange( 0.0, 3.0 ),
    m_OverlapThreshold( 1.0 ),
    m_RedshiftStep( 0.0001 ),
    m_MaxCorrelationExtremumCount( 5 )
{

}

CProcessFlowContext::~CProcessFlowContext()
{

}

bool CProcessFlowContext::Init( const char* spectrumPath, const char* noisePath, const char* templateCatalogPath, const char* rayCatalogPath, const SParam& params )
{
    m_Spectrum = new CSpectrum();

    CSpectrumIOGenericReader reader;
    Bool rValue = reader.Read( spectrumPath, *m_Spectrum );
    if( !rValue )
    {
        Log.LogError("Failed to read input spectrum file: (%s)", spectrumPath );
        m_Spectrum = NULL;
        return false;
    }

    // add noise if any or add flat noise
    if( noisePath == NULL )
    {
        CNoiseFlat noise;
        noise.SetStatErrorLevel( 1.0 );
        if (! noise.AddNoise( *m_Spectrum ) )
        {
            Log.LogError( "Failed to apply flat noise" );
            return false;
        }
    }
    else
    {
        CNoiseFromFile noise;
        if( ! noise.SetNoiseFilePath( noisePath ) )
        {
            Log.LogError("Failled to load noise spectrum");
            return false;
        }

        if( ! noise.AddNoise( *m_Spectrum ) )
        {
            Log.LogError( "Failed to apply noise from spectrum: %s", noisePath );
            return false;
        }
    }

    // Smooth flux
    m_Spectrum->GetFluxAxis().ApplyMeanSmooth( params.smoothWidth );

    // Load template catalog
    m_TemplateCatalog = new CTemplateCatalog;
    rValue = m_TemplateCatalog->Load( templateCatalogPath );
    if( !rValue )
    {
        Log.LogError("Failed to load template catalog: (%s)", templateCatalogPath );
        m_TemplateCatalog = NULL;
        return false;
    }

    // Load ray catalog
    m_RayCatalog = new CRayCatalog;
    rValue = m_RayCatalog->Load( rayCatalogPath );
    if( !rValue )
    {
        Log.LogError("Failed to load ray catalog: (%s)", rayCatalogPath );
        m_RayCatalog = NULL;
        return false;
    }

    m_LambdaRanges = params.lambdaRange;
    m_RedshiftRange = params.redshiftRange;
    m_RedshiftStep = params.redshiftStep;
    m_TemplateCategoryList = params.templateCategoryList;
    m_OverlapThreshold = params.overlapThreshold;

    return true;
}

const TFloat64Range& CProcessFlowContext::GetLambdaRange() const
{
    return m_LambdaRanges;
}

const TFloat64Range& CProcessFlowContext::GetRedshiftRange() const
{
    return m_RedshiftRange;
}

CSpectrum& CProcessFlowContext::GetSpectrum()
{
    return *m_Spectrum;
}

CTemplateCatalog& CProcessFlowContext::GetTemplateCatalog()
{
    return *m_TemplateCatalog;
}

CRayCatalog& CProcessFlowContext::GetRayCatalog()
{
    return *m_RayCatalog;
}

Float64 CProcessFlowContext::GetFineGrainedCorrelationRadius() const
{
    return m_FineGrainedCorrelationRadius;
}

Float64 CProcessFlowContext::GetOverlapThreshold() const
{
    return m_OverlapThreshold;
}

Float64 CProcessFlowContext::GetRedshiftStep() const
{
    return m_RedshiftStep;
}


const CProcessFlowContext::TTemplateCategoryList& CProcessFlowContext::GetTemplateCategoryList() const
{
    return m_TemplateCategoryList;
}

Float64 CProcessFlowContext::GetMaxCorrelationExtremumCount() const
{
    return m_MaxCorrelationExtremumCount;
}

Bool CProcessFlowContext::AddCorrelationResult( const CTemplate& tpl, const CRedshifts& redshifts, const TFloat64List& merits )
{
    if( m_CorrelationResult.find( tpl.GetName() ) != m_CorrelationResult.end() )
    {
        return false;
    }

    m_CorrelationResult[ tpl.GetName() ] = SCorrelationResult();
    m_CorrelationResult[ tpl.GetName() ].Merits = merits;
    m_CorrelationResult[ tpl.GetName() ].Redshifts = redshifts;

    return true;
}

Bool CProcessFlowContext::GetBestCorrelationResult( Float64& redshift, Float64& merit, std::string& tplName ) const
{
    Int32 maxIndex = 0;
    TCorrelationResults::const_iterator maxIt = m_CorrelationResult.end();
    TCorrelationResults::const_iterator it = m_CorrelationResult.begin();

    Float64 min = DBL_MAX ;
    for( it = m_CorrelationResult.begin(); it != m_CorrelationResult.end(); it++ )
    {
        const SCorrelationResult& r = (*it).second;
        for( Int32 i=0; i<r.Merits.size(); i++ )
        {
            if( r.Merits[i] < min )
            {
                min = r.Merits[i];
                maxIndex = i;
                maxIt = it;
            }
        }
    }


    if( maxIt != m_CorrelationResult.end() )
    {
        const SCorrelationResult& r = (*maxIt).second;
        redshift = r.Redshifts[maxIndex];
        merit = r.Merits[maxIndex];
        tplName = (*maxIt).first;
        return true;
    }
    return false;

}
