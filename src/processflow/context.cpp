#include <epic/redshift/processflow/context.h>

#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/io/genericreader.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/noise/flat.h>
#include <epic/redshift/noise/fromfile.h>
#include <epic/redshift/ray/ray.h>
#include <epic/redshift/ray/catalog.h>
#include <epic/redshift/continuum/median.h>
#include <epic/core/log/log.h>
#include <epic/core/debug/assert.h>

#include <stdio.h>
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
    smoothWidth = 0;
    overlapThreshold = 0.9;
    correlationExtremumCount = 5;

    //method = nMethod_LineMatching;
    //method = nMethod_BlindSolve;
    //method = nMethod_FullSolve;
    method = nMethod_DecisionalTree7;

    templateCategoryList.push_back( CTemplate::nCategory_Emission );
    templateCategoryList.push_back( CTemplate::nCategory_Galaxy );
    templateCategoryList.push_back( CTemplate::nCategory_Star );
    templateCategoryList.push_back( CTemplate::nCategory_Qso );
}

CProcessFlowContext::CProcessFlowContext()
{

}

CProcessFlowContext::~CProcessFlowContext()
{

}


bool CProcessFlowContext::Init( const char* spectrumPath, const char* noisePath, const CTemplateCatalog& templateCatalog, const CRayCatalog& rayCatalog, const SParam& params  )
{
    m_SpectrumName = bfs::path( spectrumPath ).stem().string();

    m_Spectrum = new CSpectrum();

    m_Params = params;

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

    // This should be moved to CProcessFlow

    // Smooth flux
    if( params.smoothWidth > 0 )
        m_Spectrum->GetFluxAxis().ApplyMeanSmooth( params.smoothWidth );


    // Compute continuum substracted spectrum
    m_SpectrumWithoutContinuum = new CSpectrum();
    *m_SpectrumWithoutContinuum = *m_Spectrum;

    m_SpectrumWithoutContinuum->RemoveContinuum<CContinuumMedian>();
    m_SpectrumWithoutContinuum->ConvertToLogScale();


    m_TemplateCatalog = ( CTemplateCatalog*) &templateCatalog;
    m_RayCatalog = ( CRayCatalog*) &rayCatalog;

    return true;
}

bool CProcessFlowContext::Init( const char* spectrumPath, const char* noisePath, const char* templateCatalogPath, const char* rayCatalogPath, const SParam& params )
{
    CRef<CTemplateCatalog> templateCatalog = new CTemplateCatalog;
    CRef<CRayCatalog> rayCatalog = new CRayCatalog;

    Bool rValue;

    // Load template catalog
    if( templateCatalogPath )
    {
        rValue = templateCatalog->Load( templateCatalogPath );
        if( !rValue )
        {
            Log.LogError("Failed to load template catalog: (%s)", templateCatalogPath );
            m_TemplateCatalog = NULL;
            return false;
        }
    }

    // Load ray catalog
    if( rayCatalogPath )
    {
        rValue = rayCatalog->Load( rayCatalogPath );
        if( !rValue )
        {
            Log.LogError("Failed to load ray catalog: (%s)", rayCatalogPath );
            m_RayCatalog = NULL;
            return false;
        }
    }

    return Init( spectrumPath, noisePath, *templateCatalog, *rayCatalog, params );
}

const CSpectrum& CProcessFlowContext::GetSpectrum() const
{
    return *m_Spectrum;
}

const CSpectrum& CProcessFlowContext::GetSpectrumWithoutContinuum() const
{
    return *m_SpectrumWithoutContinuum;
}

const CTemplateCatalog& CProcessFlowContext::GetTemplateCatalog() const
{
    return *m_TemplateCatalog;
}

const CRayCatalog& CProcessFlowContext::GetRayCatalog() const
{
    return *m_RayCatalog;
}

const CProcessFlowContext::SParam& CProcessFlowContext::GetParams() const
{
    return m_Params;
}

void CProcessFlowContext::StoreResult( TResultsMap& map, const char* name, const COperatorResult& result )
{
    TResultsMap::iterator it = map.find( name );
    if( it != map.end() )
    {
        DebugError( "Resulat already exist");
        return;
    }

    map[ name ] = &result;
}

Void CProcessFlowContext::StorePerTemplateResult( const CTemplate& t, const char* name, const COperatorResult& result )
{
    TPerTemplateResultsMap::iterator it = m_PerTemplateResults.find( t.GetName() );
    if( it == m_PerTemplateResults.end() )
    {
        m_PerTemplateResults[ t.GetName() ] = TResultsMap();
    }

    StoreResult( m_PerTemplateResults[ t.GetName() ], name, result );
}

Void CProcessFlowContext::StoreGlobalResult( const char* name, const COperatorResult& result )
{
    StoreResult( m_GlobalResults, name, result );
}

const COperatorResult* CProcessFlowContext::GetPerTemplateResult( const CTemplate& t, const char* name ) const
{
    TPerTemplateResultsMap::const_iterator it1 = m_PerTemplateResults.find( t.GetName() );
    if( it1 != m_PerTemplateResults.end() )
    {
        const TResultsMap& m = (*it1).second;
        TResultsMap::const_iterator it2 = m.find( name );
        if( it2 != m.end() )
        {
            return (*it2).second;
        }
    }

    return NULL;
}

TOperatorResultMap CProcessFlowContext::GetPerTemplateResult( const char* name ) const
{
    TOperatorResultMap map;
    TPerTemplateResultsMap::const_iterator it;
    for( it = m_PerTemplateResults.begin(); it != m_PerTemplateResults.end(); ++it )
    {
        std::string tplName = (*it).first;

        const TResultsMap& m = (*it).second;
        TResultsMap::const_iterator it2 = m.find( name );
        if( it2 != m.end() )
        {
            map[tplName] = (*it2).second;
        }
    }

    return map;
}

const COperatorResult* CProcessFlowContext::GetGlobalResult( const char* name ) const
{
    TResultsMap::const_iterator it = m_GlobalResults.find( name );
    if( it != m_GlobalResults.end() )
    {
        return (*it).second;
    }

    return NULL;
}
