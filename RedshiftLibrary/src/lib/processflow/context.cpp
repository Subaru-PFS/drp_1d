#include <RedshiftLibrary/processflow/context.h>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/debug/assert.h>

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

bool CProcessFlowContext::Init( std::shared_ptr<CSpectrum> spectrum,
                                const std::string processingID,
                                std::shared_ptr<const CTemplateCatalog> templateCatalog,
                                std::shared_ptr<const CRayCatalog> rayCatalog,
                                std::shared_ptr<CParameterStore> paramStore,
                                std::shared_ptr<CClassifierStore> zqualStore )
{
    Log.LogInfo("Processing context initialization");

    m_ClassifierStore = zqualStore;

    m_TemplateCatalog = templateCatalog;
    m_RayCatalog = rayCatalog;
    m_ParameterStore = paramStore;
    m_ResultStore = std::shared_ptr<COperatorResultStore>( new COperatorResultStore );

    // Spectrum initialization
    m_Spectrum = spectrum;
    InitSpectrum();

    // calzetti ISM & Meiksin IGM initialiation
    std::string calibrationPath;
    m_ParameterStore->Get( "calibrationDir", calibrationPath );
    InitIsmIgm(calibrationPath);

    // DataStore initialization
    m_DataStore = std::shared_ptr<CDataStore>( new CDataStore( *m_ResultStore, *m_ParameterStore ) );
    m_DataStore->SetSpectrumName( m_Spectrum->GetName() );
    m_DataStore->SetProcessingID( processingID );

    Log.LogInfo("Processing context is ready");

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

void CProcessFlowContext::InitSpectrum()
{
    // Smooth flux
    Int64 smoothWidth;
    m_ParameterStore->Get( "smoothWidth", smoothWidth, 0 );
    if( smoothWidth > 0 )
        m_Spectrum->GetFluxAxis().ApplyMeanSmooth(smoothWidth);

    // Continuum removal params
    std::string medianRemovalMethod;
    m_ParameterStore->Get( "continuumRemoval.method", medianRemovalMethod, "IrregularSamplingMedian" );
    m_Spectrum->SetContinuumEstimationMethod(medianRemovalMethod);

    Float64 medianKernelWidth;
    m_ParameterStore->Get( "continuumRemoval.medianKernelWidth", medianKernelWidth, 75.0 );
    m_Spectrum->SetMedianWinsize(medianKernelWidth);

    Float64 nscales;
    m_ParameterStore->Get( "continuumRemoval.decompScales", nscales, 6.0 );
    m_Spectrum->SetDecompScales((Int32)nscales);

    std::string dfBinPath;
    m_ParameterStore->Get( "continuumRemoval.binPath", dfBinPath, "absolute_path_to_df_binaries_here" );
    m_Spectrum->SetWaveletsDFBinPath(dfBinPath);
}


void CProcessFlowContext::InitIsmIgm(const std::string & calibrationPath)
{
    //ISM
    auto ismCorrectionCalzetti = std::make_shared<CSpectrumFluxCorrectionCalzetti>();
    ismCorrectionCalzetti->Init(calibrationPath, 0.0, 0.1, 10);
    //IGM
    auto igmCorrectionMeiksin = std::make_shared<CSpectrumFluxCorrectionMeiksin>();
    igmCorrectionMeiksin->Init(calibrationPath);

    //push in all galaxy templates
    TTemplateRefList  TplList = m_TemplateCatalog->GetTemplate(TStringList{"galaxy"});
    for (auto tpl : TplList)
    {
        tpl->m_ismCorrectionCalzetti = ismCorrectionCalzetti;
        tpl->m_igmCorrectionMeiksin = igmCorrectionMeiksin;
    }   
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

bool CProcessFlowContext::correctSpectrum(Float64 LambdaMin,  Float64 LambdaMax)
{
    return m_Spectrum->correctSpectrum( LambdaMin, LambdaMax );
}

const CSpectrum& CProcessFlowContext::GetSpectrum() const
{
    return *m_Spectrum;
}

const CTemplateCatalog& CProcessFlowContext::GetTemplateCatalog() const
{
    return *m_TemplateCatalog;
}

const CRayCatalog& CProcessFlowContext::GetRayCatalog() const
{
    return *m_RayCatalog;
}
