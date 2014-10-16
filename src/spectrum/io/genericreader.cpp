#include <epic/redshift/spectrum/io/genericreader.h>

#include <epic/redshift/spectrum/io/asciireader.h>
#include <epic/redshift/spectrum/io/fitsreader.h>

#include <boost/filesystem.hpp>

using namespace __NS__;
using namespace std;

namespace bfs = boost::filesystem;

CSpectrumIOGenericReader::CSpectrumIOGenericReader()
{

}

CSpectrumIOGenericReader::~CSpectrumIOGenericReader()
{

}

Bool CSpectrumIOGenericReader::CanRead( const char* filePath )
{
    bfs::path path( filePath );

    string fileExtension = path.extension().string();

    if( fileExtension == ".fits" )
    {
        return true;
    }
    else if( fileExtension == ".txt" )
    {
        return true;
    }

    return false;
}

Bool CSpectrumIOGenericReader::Read( const char* filePath, CSpectrum& spectrum )
{
    bfs::path path( filePath );

    string fileExtension = path.extension().string();

    if( !bfs::exists( path ) )
    {
        return false;
    }

    if( fileExtension == ".fits" )
    {
        CSpectrumIOFitsReader reader;

        return reader.Read( filePath, spectrum );
    }
    else if( fileExtension == ".txt" )
    {
        CSpectrumIOAsciiReader reader;

        return reader.Read( filePath, spectrum );
    }

    return false;
}
