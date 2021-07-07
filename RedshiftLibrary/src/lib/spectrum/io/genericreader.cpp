#include "RedshiftLibrary/spectrum/io/genericreader.h"
#include "RedshiftLibrary/spectrum/io/asciireader.h"
#include "RedshiftLibrary/spectrum/io/fitsreader.h"
#include "RedshiftLibrary/noise/flat.h"
#include "RedshiftLibrary/noise/fromfile.h"
#include "RedshiftLibrary/log/log.h"

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
void CSpectrumIOGenericReader::Read( const char* filePath, CSpectrum& spectrum )
{
    bfs::path path( filePath );

    string fileExtension = path.extension().string();

    if( !bfs::exists( path ) )
    {
      Log.LogError("File doesn't exist : %s", filePath);
      throw runtime_error("File doesn't exist");
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
      Log.LogError("%s : Unknown file type", filePath);
      throw runtime_error("Unknown file type");
    }
}
