#include <epic/redshift/noise/fromfile.h>

#include <epic/redshift/spectrum/axis.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/io/genericreader.h>

using namespace NSEpic;
using namespace std;

IMPLEMENT_MANAGED_OBJECT( CNoiseFromFile )

CNoiseFromFile::CNoiseFromFile()
{
}

CNoiseFromFile::~CNoiseFromFile()
{
}

Bool CNoiseFromFile::SetNoiseFilePath( const char* filePath )
{
    m_NoiseSpectrum = std::shared_ptr<CSpectrum>( new CSpectrum() );
    CSpectrumIOGenericReader reader;

    if( reader.Read( filePath, *m_NoiseSpectrum ) )
    {
        return true;
    }

    return false;
}

Bool CNoiseFromFile::AddNoise( CSpectrum& s1 ) const
{
    if( m_NoiseSpectrum == NULL )
        return false;

    if( s1.GetSampleCount() != m_NoiseSpectrum->GetSampleCount() )
        return false;

    Float64* dstError = s1.GetFluxAxis().GetError();
    Float64* srcError = m_NoiseSpectrum->GetFluxAxis().GetSamples();

    for( UInt32 i=0; i<s1.GetFluxAxis().GetSamplesCount(); i++ )
    {
        dstError[i] = srcError[i];
    }

    return true;
}
