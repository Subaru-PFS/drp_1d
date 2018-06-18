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

    Void SetNoiseFilePath( const char* filePath, CSpectrumIOReader& noise_reader);

    Void AddNoise( CSpectrum& s1 ) const;

private:

    CSpectrum     m_NoiseSpectrum;
    Bool m_initialized;
};


}

#endif
