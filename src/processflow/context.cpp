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

#include <epic/redshift/operator/blindsolveresult.h>
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
    SetSpectrumName( bfs::path( spectrumPath ).stem().string().c_str() );

    m_Spectrum = new CSpectrum();

    m_Params = params;

    m_dtreepath = nDtreePath_None;
    m_dtreepathnum = -1.0;

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

std::string CProcessFlowContext::GetMethodName( EMethod method )
{
    std::string methodStr = "Invalid method name";

    if(method== CProcessFlowContext::nMethod_BlindSolve){
        methodStr = "BlindSolve";
    } else if (method == CProcessFlowContext::nMethod_LineMatching){
        methodStr = "LineMatching";
    } else if (method == CProcessFlowContext::nMethod_DecisionalTree7){
        methodStr = "DecisionalTree7";
    } else if (method == CProcessFlowContext::nMethod_FullSolve){
        methodStr = "FullSolve";
    }
    return methodStr;
}

void CProcessFlowContext::SaveRedshift( const char* dir )
{
    // Append best redshift, best merit and selected template
    {
        std::fstream outputStream;
        // Save result at root of output directory
        CreateResultStorage( outputStream, bfs::path( "redshift.csv" ), bfs::path( dir ) );

        // Dump output of blindsolve
        if(m_Params.method == CProcessFlowContext::nMethod_BlindSolve)
        {
            Float64 redshift;
            Float64 merit;
            std::string tplName;

            const CBlindSolveResult* blindSolveResult = (CBlindSolveResult*)GetGlobalResult( "blindsolve" );
            blindSolveResult->GetBestFitResult( *this, redshift, merit, tplName );


            outputStream  << m_SpectrumName  << "\t"
                        << redshift << "\t"
                        << merit << "\t"
                        << tplName << std::endl;
        }
        // Dump output of full solve
        else if(m_Params.method == CProcessFlowContext::nMethod_FullSolve)
        {
            Float64 redshift;
            Float64 merit;
            std::string tplName;

            const CBlindSolveResult* blindSolveResult = (CBlindSolveResult*)GetGlobalResult( "blindsolve" );
            blindSolveResult->GetBestCorrelationPeakResult( *this, redshift, merit, tplName );


            outputStream  << m_SpectrumName  << "\t"
                        << redshift << "\t"
                        << merit << "\t"
                        << tplName << std::endl;
        }// Dump output of DecisionalTree7
        else if(m_Params.method == CProcessFlowContext::nMethod_DecisionalTree7)
        {
            Float64 redshift;
            Float64 merit;
            std::string tplName;

            const CBlindSolveResult* blindSolveResult = (CBlindSolveResult*)GetGlobalResult( "blindsolve" );
            if(m_dtreepath == CProcessFlowContext::nDtreePath_BlindSolve || m_dtreepath == CProcessFlowContext::nDtreePath_OnlyFit ){
                blindSolveResult->GetBestFitResult( *this, redshift, merit, tplName );
            }else if(m_dtreepath == CProcessFlowContext::nDtreePath_OnlyCorrelation){
                blindSolveResult->GetBestCorrelationPeakResult( *this, redshift, merit, tplName );
            }else{
                redshift = -1.0;
                merit = -1.0;
                tplName = "None";
            }

            outputStream  << m_SpectrumName  << "\t"
                        << redshift << "\t"
                        << merit << "\t"
                        << tplName << "\t"
                        << m_dtreepathnum  << "\t"
                        << std::endl;
        }
        // Dump output of raymatching
        else if(m_Params.method  == CProcessFlowContext::nMethod_LineMatching)
        {
            Float64 redshift = -1;
            Int32 matchingNum = -1;

            const CRayMatchingResult* rayMatchingResult = (CRayMatchingResult*)GetGlobalResult( "raymatching" );
            if(rayMatchingResult != NULL){
                rayMatchingResult->GetBestRedshift( redshift, matchingNum );
            }

            outputStream.precision(6);
            outputStream  << m_SpectrumName  << "\t"
                        << redshift << "\t"
                        << matchingNum << "\t"
                        << "Ray Matching" << std::endl;
        }

    }
    return;
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
