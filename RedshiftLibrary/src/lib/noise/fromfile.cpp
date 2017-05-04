#include <RedshiftLibrary/noise/fromfile.h>

#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/io/genericreader.h>

using namespace NSEpic;
using namespace std;


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
