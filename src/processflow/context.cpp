#include <epic/redshift/processflow/context.h>

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

#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CProcessFlowContext )

CProcessFlowContext::SParam::SParam()
{
    redshiftRange = TFloat64Range( 0.0, 5.0 );
    redshiftStep = 0.0001;
    lambdaRange = TFloat64Range( 3000.0, 9600.0 );
    smoothWidth = 0;
    overlapThreshold = 1.0;
    correlationExtremumCount = 5;

    //method = nMethod_Correlation;
    //method = nMethod_Chisquare;
    //method = nMethod_LineMatching;
    //method = nMethod_LineMatching2;
    //method = nMethod_BlindSolve;
    //method = nMethod_FullSolve;
    //method = nMethod_DecisionalTree7;
    //method = nMethod_DecisionalTreeA;
    method = nMethod_DecisionalTreeB;

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
    SetSpectrumName( bfs::path( spectrumPath ).stem().string().c_str() );

    m_Spectrum = new CSpectrum();
    m_Spectrum->SetName(bfs::path( spectrumPath ).stem().string().c_str() );

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

    m_SpectrumWithoutContinuum->RemoveContinuum<CContinuumIrregularSamplingMedian>();
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

std::string CProcessFlowContext::GetMethodName( EMethod method )
{
    std::string methodStr = "Invalid method name";

    if(method== CProcessFlowContext::nMethod_BlindSolve){
        methodStr = "BlindSolve";
    } else if (method == CProcessFlowContext::nMethod_Correlation){
        methodStr = "CorrelationSolve";
    } else if (method == CProcessFlowContext::nMethod_Chisquare){
        methodStr = "ChisquareSolve";
    } else if (method == CProcessFlowContext::nMethod_LineMatching){
        methodStr = "LineMatching";
    } else if (method == CProcessFlowContext::nMethod_LineMatching2){
        methodStr = "LineMatching2";
    } else if (method == CProcessFlowContext::nMethod_LineModel){
        methodStr = "LineModel";
    } else if (method == CProcessFlowContext::nMethod_DecisionalTree7){
        methodStr = "DecisionalTree7";
    } else if (method == CProcessFlowContext::nMethod_DecisionalTreeA){
        methodStr = "DecisionalTreeA";
    } else if (method == CProcessFlowContext::nMethod_DecisionalTreeB){
        methodStr = "DecisionalTreeB";
    } else if (method == CProcessFlowContext::nMethod_FullSolve){
        methodStr = "FullSolve";
    }
    return methodStr;
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
