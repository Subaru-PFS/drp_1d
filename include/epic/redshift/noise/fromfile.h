#ifndef _REDSHIFT_NOISE_FROMFILE_
#define _REDSHIFT_NOISE_FROMFILE_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/ref.h>
#include <epic/redshift/noise/noise.h>

namespace __NS__
{

class CSpectrum;

class CNoiseFromFile : public CNoise
{

    DEFINE_MANAGED_OBJECT( CNoiseFromFile )

public:

    CNoiseFromFile( );
    ~CNoiseFromFile();

    Bool SetNoiseFilePath( const char* filePath );

    Bool AddNoise( CSpectrum& s1 ) const; 

private:

    CRef<CSpectrum>     m_NoiseSpectrum;

};


}

#endif
