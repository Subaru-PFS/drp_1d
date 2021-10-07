// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"

#include "RedshiftLibrary/spectrum/io/asciireader.h"

#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/axis.h"

#include <boost/algorithm/string/predicate.hpp>
#include <sstream>

using namespace NSEpic;
using namespace std;
namespace bfs = boost::filesystem;

/**
 *
 */
CSpectrumIOAsciiReader::CSpectrumIOAsciiReader()
{

}

/**
 *
 */
CSpectrumIOAsciiReader::~CSpectrumIOAsciiReader()
{

}

/**
 *
 */
void CSpectrumIOAsciiReader::Read( const char* filePath, CSpectrum& spectrum )
{
  //Uncomment below when --verbose works properly.
  Log.LogDebug ( "Parsing ASCII file %s.", filePath );
  if( !bfs::exists( filePath ) )
    {
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"Read: Path for spectrum file does not exist. : "<< filePath);
    }

  bfs::ifstream file;
  file.open( filePath );

  if( !IsAsciiDataFile( file ) )
    {
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"Read: file is not ASCII : "<<filePath);
    }

  Int32 length = GetAsciiDataLength( file );
  if( length == -1 )
    {
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"Read: file length == -1 : "<< filePath);
    }

  CSpectrumFluxAxis spcFluxAxis(length);
  CSpectrumSpectralAxis spcSpectralAxis(length);

  Int32 i = 0;
  file.clear();
  file.seekg( 0 );

  for( std::string line; std::getline( file, line ); )
    {
      if( !boost::starts_with( line, "#" ) )
        {
	  std::istringstream iss( line );
	  Float64 x, y;
	  iss >> x >> y;
	  spcSpectralAxis[i] = x;
	  spcFluxAxis[i] = y;
	  i++;
        }
    }
  file.close();
  spectrum.SetSpectralAndFluxAxes(std::move(spcSpectralAxis),std::move(spcFluxAxis));
  
  Log.LogDebug ( "File contents read as ASCII characters." );
}

/**
 *
 */
Bool CSpectrumIOAsciiReader::IsAsciiDataFile( bfs::ifstream& file  )
{
    return true;
}


/**
 * Calculates the number of non-commented-out ASCII characters in the file.
 */
Int32 CSpectrumIOAsciiReader::GetAsciiDataLength( bfs::ifstream& file )
{
  Log.LogDebug ( "GetAsciiDataLength: started." );
    Int32 len = 0;
    for( std::string line; std::getline( file, line ); )
    {
        if( !boost::starts_with( line, "#" ) )
        {
            std::istringstream iss( line );
            Float64 x, y;
            iss >> x >> y;
            if( iss.rdstate() & std::ifstream::failbit )
            {
                file.clear();
                file.seekg ( 0 );
		Log.LogError( "GetAsciiDataLength: iss failbit was set." );
                return -1;
            }
            len++;
        }
    }

    file.clear();
    file.seekg ( 0 );
    if( file.rdstate() & std::ifstream::failbit )
      {
	Log.LogError( "GetAsciiDataLength: file failbit was set." );
        return -1;
      }
    Log.LogDebug( "GetAsciiDataLength: return %d.", len );
    return len;
}
