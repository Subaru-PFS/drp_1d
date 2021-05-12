#ifndef _REDSHIFT_NOISE_FROMFILE_
#define _REDSHIFT_NOISE_FROMFILE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/noise/noise.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/io/reader.h>

#include <memory>

namespace NSEpic
{

class CNoiseFromFile : public CNoise
{

public:

    CNoiseFromFile( );
    ~CNoiseFromFile();

    void SetNoiseFilePath( const char* filePath, CSpectrumIOReader& noise_reader);

    void AddNoise( CSpectrum& s1 ) const;

private:

    CSpectrum     m_NoiseSpectrum;//a weird spectra with no real flux, only error but considered as a flux!
    Bool m_initialized;
};


}

#endif
