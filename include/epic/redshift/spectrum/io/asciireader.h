#ifndef _REDSHIFT_SPECTRUM_IO_ASCIIREADER_
#define _REDSHIFT_SPECTRUM_IO_ASCIIREADER_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/io/reader.h>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <vector>

namespace __NS__
{

class CSpectrumIOAsciiReader : public CSpectrumIOReader
{

public:

    CSpectrumIOAsciiReader();
    ~CSpectrumIOAsciiReader();

    Bool Read( const char* filePath, CSpectrum& s );

private:

    Bool    IsAsciiDataFile( boost::filesystem::ifstream& file );
    Int32   GetAsciiDataLength( boost::filesystem::ifstream& file );
};


}

#endif
