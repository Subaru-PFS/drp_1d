#include <RedshiftLibrary/spectrum/io/genericreader.h>
#include <RedshiftLibrary/spectrum/io/asciireader.h>
#include <RedshiftLibrary/spectrum/io/fitsreader.h>
#include <RedshiftLibrary/noise/flat.h>
#include <RedshiftLibrary/noise/fromfile.h>

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
 * Will return a call to reader.Read if the file extension is fits, txt or dat. Throw otherwise.
 */
Void CSpectrumIOGenericReader::Read( const char* filePath, std::shared_ptr<CSpectrum> spectrum )
{
    bfs::path path( filePath );

    string fileExtension = path.extension().string();

    if( !bfs::exists( path ) )
    {
      throw string("File doesn't exist :") + filePath;
    }

    if( fileExtension == ".fits" )
    {
        CSpectrumIOFitsReader reader;
        reader.Read( filePath, spectrum );
    }
    else if( fileExtension == ".txt" || fileExtension == ".dat"  )
    {
        CSpectrumIOAsciiReader reader;
        reader.Read( filePath, spectrum );
    } else {
      throw string(filePath) + string(" : Unknown file type");
    }
}
