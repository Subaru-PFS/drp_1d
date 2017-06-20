#include <RedshiftLibrary/spectrum/io/genericreader.h>

#include <RedshiftLibrary/spectrum/io/asciireader.h>
#include <RedshiftLibrary/spectrum/io/fitsreader.h>

#include <boost/filesystem.hpp>

using namespace NSEpic;
using namespace std;

namespace bfs = boost::filesystem;

/**
 * Empty constructor.
 */
CSpectrumIOGenericReader::CSpectrumIOGenericReader()
{

}

/**
 * Empty destructor.
 */
CSpectrumIOGenericReader::~CSpectrumIOGenericReader()
{

}

/**
 * Will return true if the file extension is either fits or txt. Returns false otherwise.
 */
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

/**
 * Will return a call to reader.Read if the file extension is fits, txt or dat. Returns false otherwise.
 */
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
    else if( fileExtension == ".txt" || fileExtension == ".dat"  )
    {
        CSpectrumIOAsciiReader reader;

        return reader.Read( filePath, spectrum );
    }

    return false;
}
