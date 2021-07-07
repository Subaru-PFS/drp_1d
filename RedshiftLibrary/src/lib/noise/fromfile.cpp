#include "RedshiftLibrary/noise/fromfile.h"

#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/io/reader.h"
#include "RedshiftLibrary/spectrum/io/genericreader.h"

using namespace NSEpic;
using namespace std;


CNoiseFromFile::CNoiseFromFile():
  m_initialized(false)
{
}

CNoiseFromFile::~CNoiseFromFile()
{
}

void CNoiseFromFile::SetNoiseFilePath( const char* filePath,
				       CSpectrumIOReader& noise_reader )
{
  m_initialized = true;
  noise_reader.Read( filePath, m_NoiseSpectrum );
}

void CNoiseFromFile::AddNoise( CSpectrum& s1 ) const
{
    if( !m_initialized )
    {
        throw runtime_error("Noise wasn't initialized");
    }

    if( s1.GetSampleCount() != m_NoiseSpectrum.GetSampleCount() )
    {
        throw runtime_error("Sample counts don't match");
    }

    CSpectrumNoiseAxis dstError = CSpectrumNoiseAxis(m_NoiseSpectrum.GetFluxAxis().GetSamplesVector());
    s1.SetErrorAxis(std::move(dstError));
}
