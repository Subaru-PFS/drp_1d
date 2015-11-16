#include <epic/redshift/processflow/context.h>

#include <epic/redshift/processflow/datastore.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/io/genericreader.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/noise/flat.h>
#include <epic/redshift/noise/fromfile.h>
#include <epic/redshift/ray/ray.h>
#include <epic/redshift/ray/catalog.h>
#include <epic/redshift/continuum/median.h>
#include <epic/redshift/continuum/irregularsamplingmedian.h>

#include <epic/core/log/log.h>
#include <epic/core/debug/assert.h>

#include <epic/redshift/method/blindsolveresult.h>
#include <epic/redshift/operator/raymatchingresult.h>

#include <stdio.h>
#include <float.h>
#include <fstream>
#include <memory>

#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CProcessFlowContext )

CProcessFlowContext::CProcessFlowContext()
{

}

CProcessFlowContext::~CProcessFlowContext()
{

}


bool CProcessFlowContext::Init( const char* spectrumPath, const char* noisePath,
                                const CTemplateCatalog& templateCatalog, const CRayCatalog& rayCatalog,
                                CParameterStore& paramStore  )
{
    m_Spectrum = std::shared_ptr<CSpectrum>( new CSpectrum() );
    m_Spectrum->SetName(bfs::path( spectrumPath ).stem().string().c_str() );

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
    Int64 smoothWidth;
    paramStore.Get( "smoothWidth", smoothWidth, 0 );
    if( smoothWidth > 0 )
        m_Spectrum->GetFluxAxis().ApplyMeanSmooth( smoothWidth );


    // Compute continuum substracted spectrum
    m_SpectrumWithoutContinuum = std::shared_ptr<CSpectrum>( new CSpectrum() );
    *m_SpectrumWithoutContinuum = *m_Spectrum;


    CContinuumIrregularSamplingMedian continuum;

    m_SpectrumWithoutContinuum->RemoveContinuum( continuum );
    m_SpectrumWithoutContinuum->ConvertToLogScale();


    m_TemplateCatalog = ( CTemplateCatalog*) &templateCatalog;
    m_RayCatalog = ( CRayCatalog*) &rayCatalog;
    m_ParameterStore = &paramStore;
    m_ResultStore = new COperatorResultStore;


    m_DataStore = new CDataStore( *m_ResultStore, *m_ParameterStore );
    m_DataStore->SetSpectrumName( bfs::path( spectrumPath ).stem().string() );


    return true;
}

bool CProcessFlowContext::Init( const char* spectrumPath, const char* noisePath,
                                const char* templateCatalogPath, const char* rayCatalogPath,
                                CParameterStore& paramStore )
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

    // Load line catalog
    if( rayCatalogPath )
    {
        rValue = rayCatalog->Load( rayCatalogPath );
        if( !rValue )
        {
            Log.LogError("Failed to load line catalog: (%s)", rayCatalogPath );
            m_RayCatalog = NULL;
            return false;
        }
    }

    return Init( spectrumPath, noisePath, *templateCatalog, *rayCatalog, paramStore );
}

CParameterStore& CProcessFlowContext::GetParameterStore()
{
    return *m_ParameterStore;
}

COperatorResultStore&  CProcessFlowContext::GetResultStore()
{
    return *m_ResultStore;
}


CDataStore& CProcessFlowContext::GetDataStore()
{
    return *m_DataStore;
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


