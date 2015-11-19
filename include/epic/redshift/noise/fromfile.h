#ifndef _REDSHIFT_NOISE_FROMFILE_
#define _REDSHIFT_NOISE_FROMFILE_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/noise/noise.h>

#include <memory>

namespace NSEpic
{

class CSpectrum;

class CNoiseFromFile : public CNoise
{

public:

    CNoiseFromFile( );
    ~CNoiseFromFile();

    Bool SetNoiseFilePath( const char* filePath );

    Bool AddNoise( CSpectrum& s1 ) const; 

private:

    std::shared_ptr<CSpectrum>     m_NoiseSpectrum;

};


}

#endif
