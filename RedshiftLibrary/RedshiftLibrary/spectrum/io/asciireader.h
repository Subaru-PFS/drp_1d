#ifndef _REDSHIFT_SPECTRUM_IO_ASCIIREADER_
#define _REDSHIFT_SPECTRUM_IO_ASCIIREADER_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/io/reader.h>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <vector>

namespace NSEpic
{

class CSpectrumIOAsciiReader : public CSpectrumIOReader
{

public:

    CSpectrumIOAsciiReader();
    ~CSpectrumIOAsciiReader();

    virtual Void Read( const char* filePath, std::shared_ptr<CSpectrum> s );

private:

    Bool    IsAsciiDataFile( boost::filesystem::ifstream& file );
    Int32   GetAsciiDataLength( boost::filesystem::ifstream& file );
};


}

#endif
