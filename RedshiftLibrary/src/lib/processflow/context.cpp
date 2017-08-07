#include <RedshiftLibrary/processflow/context.h>

#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/io/genericreader.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/noise/flat.h>
#include <RedshiftLibrary/noise/fromfile.h>
#include <RedshiftLibrary/ray/ray.h>
#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/continuum/median.h>
#include <RedshiftLibrary/continuum/waveletsdf.h>
#include <RedshiftLibrary/continuum/irregularsamplingmedian.h>
#include <RedshiftLibrary/continuum/indexes.h>
#include <RedshiftLibrary/continuum/indexesresult.h>


#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/debug/assert.h>

#include <RedshiftLibrary/method/blindsolveresult.h>
#include <RedshiftLibrary/operator/raymatchingresult.h>
#include <RedshiftLibrary/operator/spectraFluxResult.h>

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
                                const std::string processingID,
                                std::shared_ptr<const CTemplateCatalog> templateCatalog,
                                std::shared_ptr<const CRayCatalog> rayCatalog,
                                std::shared_ptr<CParameterStore> paramStore,
                                std::shared_ptr<CClassifierStore> zqualStore  )
{
    m_ClassifierStore=zqualStore;

    m_Spectrum = std::shared_ptr<CSpectrum>( new CSpectrum() );
    m_Spectrum->SetName(bfs::path( spectrumPath ).stem().string().c_str() );
    m_Spectrum->SetFullPath(bfs::path( spectrumPath ).string().c_str() );
    //Log.LogInfo("Setting spectrum name: (%s)", bfs::path( spectrumPath ).stem().string().c_str() );

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
            Log.LogError("Failed to load noise spectrum");
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
    paramStore->Get( "smoothWidth", smoothWidth, 0 );
    if( smoothWidth > 0 )
        m_Spectrum->GetFluxAxis().ApplyMeanSmooth( smoothWidth );


    // Compute continuum substracted spectrum
    m_SpectrumWithoutContinuum = std::shared_ptr<CSpectrum>( new CSpectrum() );
    *m_SpectrumWithoutContinuum = *m_Spectrum;

    std::string medianRemovalMethod;
    paramStore->Get( "continuumRemoval.method", medianRemovalMethod, "IrregularSamplingMedian" );
    //paramStore->Get( "continuumRemoval.method", medianRemovalMethod, "waveletsDF" );

    const char*     nameBaseline;		// baseline filename
    Log.LogInfo( "Continuum estimation: using %s", medianRemovalMethod.c_str() );
    if( medianRemovalMethod== "IrregularSamplingMedian")
    {
        nameBaseline = "preprocess/baselineISMedian";
        CContinuumIrregularSamplingMedian continuum;
        Float64 opt_medianKernelWidth;
        paramStore->Get( "continuumRemoval.medianKernelWidth", opt_medianKernelWidth, 75 );
        continuum.SetMedianKernelWidth(opt_medianKernelWidth);
        continuum.SetMeanKernelWidth(opt_medianKernelWidth);
        m_SpectrumWithoutContinuum->RemoveContinuum( continuum );
        m_SpectrumWithoutContinuum->SetMedianWinsize(opt_medianKernelWidth);
    }else if( medianRemovalMethod== "Median")
    {
        nameBaseline = "preprocess/baselineMedian";
        CContinuumMedian continuum;
        Float64 opt_medianKernelWidth;
        paramStore->Get( "continuumRemoval.medianKernelWidth", opt_medianKernelWidth, 75 );
        continuum.SetMedianKernelWidth(opt_medianKernelWidth);
        m_SpectrumWithoutContinuum->RemoveContinuum( continuum );
        m_SpectrumWithoutContinuum->SetMedianWinsize(opt_medianKernelWidth);

    }else if( medianRemovalMethod== "waveletsDF")
    {
        nameBaseline = "preprocess/baselineDF";
        Int64 nscales;
        paramStore->Get( "continuumRemoval.decompScales", nscales, 6);
        std::string dfBinPath;
        paramStore->Get( "continuumRemoval.binPath", dfBinPath, "absolute_path_to_df_binaries_here");
        CContinuumDF continuum(dfBinPath);
        m_SpectrumWithoutContinuum->SetDecompScales(nscales);
        bool ret = m_SpectrumWithoutContinuum->RemoveContinuum( continuum );
        if( !ret ) //doesn't seem to work. TODO: check that the df errors lead to a ret=false value
        {
            Log.LogError( "Failed to apply continuum substraction for spectrum: %s", spectrumPath );
            return false;
        }
    }else if( medianRemovalMethod== "raw")
    {
        nameBaseline = "preprocess/baselineRAW";
        CSpectrumFluxAxis& spcFluxAxis = m_SpectrumWithoutContinuum->GetFluxAxis();
        spcFluxAxis.SetSize( m_Spectrum->GetSampleCount() );
        CSpectrumSpectralAxis& spcSpectralAxis = m_SpectrumWithoutContinuum->GetSpectralAxis();
        spcSpectralAxis.SetSize( m_Spectrum->GetSampleCount()  );


        for(Int32 k=0; k<m_SpectrumWithoutContinuum->GetSampleCount(); k++)
        {
            spcFluxAxis[k] = 0.0;
        }
    }else if( medianRemovalMethod== "zero")
    {
        nameBaseline = "preprocess/baselineZERO";
        CSpectrumFluxAxis& spcFluxAxis = m_SpectrumWithoutContinuum->GetFluxAxis();
        spcFluxAxis.SetSize( m_Spectrum->GetSampleCount() );
        CSpectrumSpectralAxis& spcSpectralAxis = m_SpectrumWithoutContinuum->GetSpectralAxis();
        spcSpectralAxis.SetSize( m_Spectrum->GetSampleCount()  );


        for(Int32 k=0; k<m_SpectrumWithoutContinuum->GetSampleCount(); k++)
        {
            spcFluxAxis[k] = m_Spectrum->GetFluxAxis()[k];
        }
    }
    Log.LogInfo("===============================================");

    //process continuum relevance
    CContinuumIndexes continuumIndexes;
    CSpectrum _spcContinuum = *m_Spectrum;
    CSpectrumFluxAxis spcfluxAxis = _spcContinuum.GetFluxAxis();
    spcfluxAxis.Subtract( m_SpectrumWithoutContinuum->GetFluxAxis() );
    CSpectrumFluxAxis& sfluxAxisPtr = _spcContinuum.GetFluxAxis();
    sfluxAxisPtr = spcfluxAxis;

    CContinuumIndexes::SContinuumRelevance continuumRelevance = continuumIndexes.getRelevance( *m_Spectrum, _spcContinuum );


    m_SpectrumWithoutContinuum->ConvertToLogScale();


    m_TemplateCatalog = templateCatalog;
    m_RayCatalog = rayCatalog;
    m_ParameterStore = paramStore;
    m_ResultStore = std::shared_ptr<COperatorResultStore>( new COperatorResultStore );


    m_DataStore = std::shared_ptr<CDataStore>( new CDataStore( *m_ResultStore, *m_ParameterStore ) );
    m_DataStore->SetSpectrumName( bfs::path( spectrumPath ).stem().string() );
    m_DataStore->SetProcessingID( processingID );

    // Save the baseline in store
    std::shared_ptr<CSpectraFluxResult> baselineResult = (std::shared_ptr<CSpectraFluxResult>) new CSpectraFluxResult();
    baselineResult->m_optio = 0;
    UInt32 len = m_Spectrum->GetSampleCount();

    baselineResult->fluxes.resize(len);
    baselineResult->wavel.resize(len);
    for( Int32 k=0; k<len; k++ )
    {
        baselineResult->fluxes[k] = (m_Spectrum->GetFluxAxis())[k] - (m_SpectrumWithoutContinuum->GetFluxAxis())[k];
        baselineResult->wavel[k]  = (m_Spectrum->GetSpectralAxis())[k];
    }
    m_DataStore->StoreScopedGlobalResult(nameBaseline, baselineResult);

    //Save the continuum relevance in store
    const char*     nameContinuumIndexesResult;		// continuum indexes filename
    nameContinuumIndexesResult = "preprocess/continuumIndexes";
    std::shared_ptr<CContinuumIndexesResult> continuumIndexesResult = (std::shared_ptr<CContinuumIndexesResult>) new CContinuumIndexesResult();
    continuumIndexesResult->SetValues(continuumRelevance.StdSpectrum, continuumRelevance.StdContinuum);
    m_DataStore->StoreScopedGlobalResult(nameContinuumIndexesResult, continuumIndexesResult);

    return true;
}

bool CProcessFlowContext::Init( const char* spectrumPath, const char* noisePath, std::string processingID,
                                const char* templateCatalogPath, const char* rayCatalogPath,
                                std::shared_ptr<CParameterStore> paramStore,
		   	        std::shared_ptr<CClassifierStore> zqualStore )
{
    std::string medianRemovalMethod;
    paramStore->Get( "continuumRemoval.method", medianRemovalMethod, "IrregularSamplingMedian" );
    //override the continuum removal for the templates :
    //medianRemovalMethod = "noRontinuumRemovalforTemplates";
    //medianRemovalMethod = "raw";

    Float64 opt_medianKernelWidth;
    paramStore->Get( "continuumRemoval.medianKernelWidth", opt_medianKernelWidth, 75 );
    Int64 opt_nscales;
    paramStore->Get( "continuumRemoval.decompScales", opt_nscales, 8);
    std::string dfBinPath;
    paramStore->Get( "continuumRemoval.binPath", dfBinPath, "absolute_path_to_df_binaries_here");

    std::shared_ptr<CTemplateCatalog> templateCatalog = std::shared_ptr<CTemplateCatalog>( new CTemplateCatalog( medianRemovalMethod, opt_medianKernelWidth, opt_nscales, dfBinPath) );
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

    return Init( spectrumPath, noisePath, processingID, templateCatalog, rayCatalog, paramStore, zqualStore );

}

CParameterStore& CProcessFlowContext::GetParameterStore()
{
    return *m_ParameterStore;
}

COperatorResultStore&  CProcessFlowContext::GetResultStore()
{
    return *m_ResultStore;
}

CClassifierStore& CProcessFlowContext::GetClassifierStore()
{
        return *m_ClassifierStore;
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


