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
#include "RedshiftLibrary/spectrum/io/genericreader.h"
#include "RedshiftLibrary/spectrum/io/asciireader.h"
#include "RedshiftLibrary/spectrum/io/fitsreader.h"
#include "RedshiftLibrary/noise/flat.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"

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
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"File doesn't exist : "<< filePath);
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
      throw GlobalException(INTERNAL_ERROR,Formatter()<<filePath<<": Unknown file type");
    }
}
