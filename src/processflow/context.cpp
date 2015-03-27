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
    smoothWidth = 0;
    overlapThreshold = 0.9;

    //templateCategoryList.push_back( CTemplate::nCategory_Emission );
    templateCategoryList.push_back( CTemplate::nCategory_Galaxy );
    templateCategoryList.push_back( CTemplate::nCategory_Star );
    templateCategoryList.push_back( CTemplate::nCategory_Qso );
}

CProcessFlowContext::CProcessFlowContext() :
    m_RedshiftRange( 0.0, 3.0 ),
    m_OverlapThreshold( 1.0 ),
    m_RedshiftStep( 0.0001 )
{

}

CProcessFlowContext::~CProcessFlowContext()
{

}


bool CProcessFlowContext::Init( const char* spectrumPath, const char* noisePath, CTemplateCatalog& templateCatalog, CRayCatalog& rayCatalog, const SParam& params  )
{
    m_SpectrumName = bfs::path( spectrumPath ).stem().string();

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
    if( params.smoothWidth > 0 )
        m_Spectrum->GetFluxAxis().ApplyMeanSmooth( params.smoothWidth );


    // Compute continuum substracted spectrum
    m_SpectrumWithoutContinuum = new CSpectrum();
    *m_SpectrumWithoutContinuum = *m_Spectrum;

    m_SpectrumWithoutContinuum->RemoveContinuum<CContinuumMedian>();
    m_SpectrumWithoutContinuum->ConvertToLogScale();


    m_TemplateCatalog = &templateCatalog;
    m_RayCatalog = &rayCatalog;


    m_LambdaRanges = params.lambdaRange;
    m_RedshiftRange = params.redshiftRange;
    m_RedshiftStep = params.redshiftStep;
    m_TemplateCategoryList = params.templateCategoryList;
    m_OverlapThreshold = params.overlapThreshold;

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

CSpectrum& CProcessFlowContext::GetSpectrumWithoutContinuum()
{
    return *m_SpectrumWithoutContinuum;
}

CTemplateCatalog& CProcessFlowContext::GetTemplateCatalog()
{
    return *m_TemplateCatalog;
}

CRayCatalog& CProcessFlowContext::GetRayCatalog()
{
    return *m_RayCatalog;
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

Bool CProcessFlowContext::AddResults( const CTemplate& tpl,
                                      const CRedshifts& selectedRedshifts,
                                      const TFloat64List& selectedMerits, const COperator::TStatusList& selectedMeritsStatus,
                                      const CRedshifts& redshifts, const TFloat64List& correlationValues  )
{
    if( m_ProcessResult.find( tpl.GetName() ) != m_ProcessResult.end() )
    {
        return false;
    }

    m_ProcessResult[ tpl.GetName() ] = SCorrelationResult();
    m_ProcessResult[ tpl.GetName() ].SelectedMerits = selectedMerits;
    m_ProcessResult[ tpl.GetName() ].SelectedMeritsStatus = selectedMeritsStatus;
    m_ProcessResult[ tpl.GetName() ].SelectedRedshifts = selectedRedshifts;
    m_ProcessResult[ tpl.GetName() ].Redshifts = redshifts;
    m_ProcessResult[ tpl.GetName() ].CorrelationValues = correlationValues;

    return true;
}

Bool CProcessFlowContext::DumpCorrelationResultsToCSV( const char* dir ) const
{
    char outputDir[256];

    sprintf(outputDir, "%s/%s", dir, m_SpectrumName.c_str() );

    if( bfs::exists( outputDir ) )
    {
        if( bfs::is_directory( outputDir ) )
        {
            // Remove everything by default
            bfs::remove_all( outputDir );
        }
        else
        {
            Log.LogError("%s exists, but is not a directory", outputDir );
            return false;
        }
    }


    if( ! bfs::create_directories( outputDir ) )
    {
        Log.LogError("Failed to create directory: %s", outputDir );
        return false;
    }



    TCorrelationResults::const_iterator it;

    for( it = m_ProcessResult.begin(); it != m_ProcessResult.end(); ++it )
    {
        char outputFileName[256];
        sprintf(outputFileName, "%s/%s.txt", outputDir, it->first.c_str() );


        FILE* f = fopen( outputFileName, "w+" );
        if( f == NULL )
        {
            Log.LogError("Failed to create file: %s", outputFileName );
            return false;
        }

        const SCorrelationResult& r = it->second;

        for( int i=0;i<r.Redshifts.GetRedshiftsCount();i++)
        {
            if( r.CorrelationValues[i] < 0.0001 ){
                fprintf( f, "%f %e\n", r.Redshifts[i], r.CorrelationValues[i]);
            }else{
                fprintf( f, "%f %f\n", r.Redshifts[i], r.CorrelationValues[i]);
            }
        }

        fclose( f );
    }


    return true;
}

Bool CProcessFlowContext::GetBestCorrelationResult( Float64& redshift, Float64& merit, std::string& tplName ) const
{
    Int32 maxIndex = 0;
    TCorrelationResults::const_iterator maxIt = m_ProcessResult.end();
    TCorrelationResults::const_iterator it = m_ProcessResult.begin();


    Float64 min = DBL_MAX ;
    for( it = m_ProcessResult.begin(); it != m_ProcessResult.end(); it++ )
    {
        const SCorrelationResult& r = (*it).second;
        for( Int32 i=0; i<r.SelectedMerits.size(); i++ )
        {
            if( r.SelectedMerits[i] < min && r.SelectedMeritsStatus[i] == COperator::nStatus_OK )
            {
                min = r.SelectedMerits[i];
                maxIndex = i;
                maxIt = it;
            }
        }
    }


    if( maxIt != m_ProcessResult.end() )
    {
        const SCorrelationResult& r = (*maxIt).second;
        redshift = r.SelectedRedshifts[maxIndex];
        merit = r.SelectedMerits[maxIndex];
        tplName = (*maxIt).first;
        return true;
    }

    return false;

}
