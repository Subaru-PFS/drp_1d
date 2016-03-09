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


CProcessFlowContext::CProcessFlowContext()
{

}

CProcessFlowContext::~CProcessFlowContext()
{

}

bool CProcessFlowContext::Init( const char* spectrumPath, const char* noisePath,
                                std::shared_ptr<const CTemplateCatalog> templateCatalog,
                                std::shared_ptr<const CRayCatalog> rayCatalog,
                                std::shared_ptr<CParameterStore> paramStore  )
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
    //check if the Spectrum is valid
    if( !m_Spectrum->IsFluxValid() ){
        Log.LogError( "Failed to validate spectrum: %s", spectrumPath );
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
            Log.LogError("Failed to load noise spectrum");
            return false;
        }

        if( ! noise.AddNoise( *m_Spectrum ) )
        {
            Log.LogError( "Failed to apply noise from spectrum: %s", noisePath );
            return false;
        }
    }
    //check if there is any zeroes in the nosie spectrum
    if( !m_Spectrum->IsNoiseValid() ){
        Log.LogError( "Failed to validate noise from spectrum: %s", noisePath );
        return false;
    }



    // This should be moved to CProcessFlow

    // Smooth flux
    Int64 smoothWidth;
    paramStore->Get( "smoothWidth", smoothWidth, 0 );
    if( smoothWidth > 0 )
        m_Spectrum->GetFluxAxis().ApplyMeanSmooth( smoothWidth );


    // Compute continuum substracted spectrum
    m_SpectrumWithoutContinuum = std::shared_ptr<CSpectrum>( new CSpectrum() );
    *m_SpectrumWithoutContinuum = *m_Spectrum;

    std::string medianRemovalMethod;
    paramStore->Get( "continuumRemoval.method", medianRemovalMethod, "IrregularSamplingMedian" );
    //ctx.GetParameterStore().Get( "continuumRemoval.method", medianRemovalMethod, "Median" );
    if( medianRemovalMethod== "IrregularSamplingMedian"){
        CContinuumIrregularSamplingMedian continuum;
        Float64 opt_medianKernelWidth;
        paramStore->Get( "continuumRemoval.medianKernelWidth", opt_medianKernelWidth, 75 );
        continuum.SetMedianKernelWidth(opt_medianKernelWidth);
        m_SpectrumWithoutContinuum->RemoveContinuum( continuum );

    }else{
        CContinuumMedian continuum;
        Float64 opt_medianKernelWidth;
        paramStore->Get( "continuumRemoval.medianKernelWidth", opt_medianKernelWidth, 75 );
        continuum.SetMedianKernelWidth(opt_medianKernelWidth);
        m_SpectrumWithoutContinuum->RemoveContinuum( continuum );
    }

    m_SpectrumWithoutContinuum->ConvertToLogScale();


    m_TemplateCatalog = templateCatalog;
    m_RayCatalog = rayCatalog;
    m_ParameterStore = paramStore;
    m_ResultStore = std::shared_ptr<COperatorResultStore>( new COperatorResultStore );


    m_DataStore = std::shared_ptr<CDataStore>( new CDataStore( *m_ResultStore, *m_ParameterStore ) );
    m_DataStore->SetSpectrumName( bfs::path( spectrumPath ).stem().string() );


    return true;
}

bool CProcessFlowContext::Init( const char* spectrumPath, const char* noisePath,
                                const char* templateCatalogPath, const char* rayCatalogPath,
                                std::shared_ptr<CParameterStore> paramStore )
{
    std::string medianRemovalMethod;
    paramStore->Get( "continuumRemoval.method", medianRemovalMethod, "IrregularSamplingMedian" );
    Float64 opt_medianKernelWidth;
    paramStore->Get( "continuumRemoval.medianKernelWidth", opt_medianKernelWidth, 75 );
    std::shared_ptr<CTemplateCatalog> templateCatalog = std::shared_ptr<CTemplateCatalog>( new CTemplateCatalog( medianRemovalMethod, opt_medianKernelWidth) );
    std::shared_ptr<CRayCatalog> rayCatalog = std::shared_ptr<CRayCatalog>(new CRayCatalog);


    Bool rValue;


    // Load template catalog
    if( templateCatalogPath )
    {
      Log.LogDebug ( "templateCatalogPath exists." );
      rValue = templateCatalog->Load( templateCatalogPath );
      if( !rValue )
        {
	  Log.LogError( "Failed to load template catalog from path: (%s)", templateCatalogPath );
	  m_TemplateCatalog = NULL;
	  return false;
        }
      Log.LogDebug ( "Template catalog loaded." );
    }

    // Load line catalog
    //std::cout << "ctx" << std::endl;
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

    return Init( spectrumPath, noisePath, templateCatalog, rayCatalog, paramStore );

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


