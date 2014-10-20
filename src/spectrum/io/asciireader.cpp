#include <epic/redshift/spectrum/io/asciireader.h>

#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/axis.h>

#include <boost/algorithm/string/predicate.hpp>
#include <sstream>

using namespace NSEpic;
using namespace std;
namespace bfs = boost::filesystem;


CSpectrumIOAsciiReader::CSpectrumIOAsciiReader()
{

}

CSpectrumIOAsciiReader::~CSpectrumIOAsciiReader()
{

}

Bool CSpectrumIOAsciiReader::Read( const char* filePath, CSpectrum& spectrum )
{
    if( !bfs::exists( filePath ) )
        return false;

    bfs::ifstream file;
    file.open( filePath );

    if( !IsAsciiDataFile( file ) )
        return false;

    Int32 length = GetAsciiDataLength( file );
    if( length == -1 )
        return false;

    CSpectrumAxis& spcFluxAxis = spectrum.GetFluxAxis();
    spcFluxAxis.SetSize( length );

    CSpectrumAxis& spcSpectralAxis = spectrum.GetSpectralAxis();
    spcSpectralAxis.SetSize( length );

    Int32 i = 0;
    file.clear();
    file.seekg(0);
    int l = file.tellg();
    for( std::string line; std::getline( file, line ); )
    {
        if( ! boost::starts_with( line, "#" ) )
        {
            std::istringstream iss( line );
            Float64 x, y;
            iss >> x >> y;
            spcSpectralAxis[i] = x;
            spcFluxAxis[i] = y;
            i++;
        }
    }

    return true;
}

Bool CSpectrumIOAsciiReader::IsAsciiDataFile( bfs::ifstream& file  )
{
    return true;
}



Int32 CSpectrumIOAsciiReader::GetAsciiDataLength( bfs::ifstream& file )
{
    Int32 len = 0;
    for( std::string line; std::getline( file, line ); )
    {
        if( ! boost::starts_with( line, "#" ) )
        {
            std::istringstream iss( line );
            Float64 x, y;
            iss >> x >> y;
            if( iss.rdstate() & std::ifstream::failbit )
            {
                file.clear();
                file.seekg ( 0 );
                return -1;
            }
            len++;
        }
    }

    file.clear();
    file.seekg ( 0 );
    if( file.rdstate() & std::ifstream::failbit )
    {
        return -1;
    }
    return len;
}


