#include <RedshiftLibrary/noise/fromfile.h>

#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/io/reader.h>
#include <RedshiftLibrary/spectrum/io/genericreader.h>

using namespace NSEpic;
using namespace std;


CNoiseFromFile::CNoiseFromFile()
{
}

CNoiseFromFile::~CNoiseFromFile()
{
}

Void CNoiseFromFile::SetNoiseFilePath( const char* filePath,
				       std::shared_ptr<CSpectrumIOReader> noise_reader )
{
    m_NoiseSpectrum = std::shared_ptr<CSpectrum>( new CSpectrum() );

    noise_reader->Read( filePath, m_NoiseSpectrum );
}

Void CNoiseFromFile::AddNoise( CSpectrum& s1 ) const
{
    if( m_NoiseSpectrum == NULL )
      throw string("Noise wasn't initialized");

    if( s1.GetSampleCount() != m_NoiseSpectrum->GetSampleCount() )
      throw string("Sample counts don't match");

    Float64* dstError = s1.GetFluxAxis().GetError();
    Float64* srcError = m_NoiseSpectrum->GetFluxAxis().GetSamples();

    for( UInt32 i=0; i<s1.GetFluxAxis().GetSamplesCount(); i++ )
    {
        dstError[i] = srcError[i];
    }
}
