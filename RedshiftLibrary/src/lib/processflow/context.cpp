#include <RedshiftLibrary/processflow/context.h>

#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/io/reader.h>
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

static void NewHandler(const char* reason,
                       const char* file,
                       int line,
                       int gsl_errno){

    Log.LogError(" gsl: %s:%d: ERROR: %s (Errtype: %s)",file, line, reason, gsl_strerror(gsl_errno));
    throw std::runtime_error("GSL Error");
    return ;
}

CProcessFlowContext::CProcessFlowContext()
{
    gsl_set_error_handler(NewHandler);
}

CProcessFlowContext::~CProcessFlowContext()
{

}

void CProcessFlowContext::SmoothFlux(std::shared_ptr<CParameterStore> paramStore)
{
    Int64 smoothWidth;
    paramStore->Get( "smoothWidth", smoothWidth, 0 );
    if( smoothWidth > 0 )
      m_Spectrum->GetFluxAxis().ApplyMeanSmooth( smoothWidth );
}

bool CProcessFlowContext::Init( std::shared_ptr<CSpectrum> spectrum,
                                const std::string processingID,
                                std::shared_ptr<const CTemplateCatalog> templateCatalog,
                                std::shared_ptr<const CRayCatalog> rayCatalog,
                                std::shared_ptr<CParameterStore> paramStore,
                                std::shared_ptr<CClassifierStore> zqualStore  )
{
    Log.LogInfo("Processing context initialization");

    m_ClassifierStore = zqualStore;

    m_Spectrum = spectrum;
    m_SpectrumWithoutContinuum = std::shared_ptr<CSpectrum>( new CSpectrum() );
    *m_SpectrumWithoutContinuum = *spectrum;

    m_TemplateCatalog = templateCatalog;
    m_RayCatalog = rayCatalog;
    m_ParameterStore = paramStore;
    m_ResultStore = std::shared_ptr<COperatorResultStore>( new COperatorResultStore );

    // DataStore initialization
    m_DataStore = std::shared_ptr<CDataStore>( new CDataStore( *m_ResultStore, *m_ParameterStore ) );
    m_DataStore->SetSpectrumName( m_Spectrum->GetName() );
    m_DataStore->SetProcessingID( processingID );

    // Smooth flux
    SmoothFlux(paramStore);

    return true;
}

bool CProcessFlowContext::Init( std::shared_ptr<CSpectrum> spectrum,
				std::string processingID,
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

    // Load template catalog
    if( templateCatalogPath )
    {
      templateCatalog->Load( templateCatalogPath );
      Log.LogDebug ( "Template catalog loaded." );
    }

    // Load line catalog
    //std::cout << "ctx" << std::endl;
    if( rayCatalogPath )
    {
        rayCatalog->Load( rayCatalogPath );
    }

    return Init( spectrum, processingID, templateCatalog, rayCatalog, paramStore, zqualStore );

}

CParameterStore& CProcessFlowContext::GetParameterStore()
{
    return *m_ParameterStore;
}

COperatorResultStore& CProcessFlowContext::GetResultStore()
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

bool CProcessFlowContext::correctSpectrum(Float64 LambdaMin,  Float64 LambdaMax)
{
    return m_Spectrum->correctSpectrum( LambdaMin, LambdaMax ) ;
}

CSpectrum& CProcessFlowContext::GetSpectrumWithoutContinuum()
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
