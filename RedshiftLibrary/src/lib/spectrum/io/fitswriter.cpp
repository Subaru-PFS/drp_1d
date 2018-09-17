#include <RedshiftLibrary/spectrum/io/fitswriter.h>

#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/axis.h>


#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;

using namespace NSEpic;
using namespace std;

/**
 * Empty constructor.
 */
CSpectrumIOFitsWriter::CSpectrumIOFitsWriter()
{

}

/**
 * Empty destructor.
 */
CSpectrumIOFitsWriter::~CSpectrumIOFitsWriter()
{

}

/**
 * 
 */
Bool CSpectrumIOFitsWriter::Write( const char* filePath, CSpectrum& spectrum )
{
    fitsfile *fptr = NULL;
    Int32 status = 0;

    Bool retv = true;

    if( bfs::exists( filePath ) )
	{
    	bfs::remove( filePath );

	}

	if( fits_create_file(&fptr, filePath, &status) )
	{
		return false;
	}


    const char* ttype[2] = {"WAVE","FLUX"};
    const char* tform[2] = {"E","E"};
    if( fits_create_tbl(fptr, BINARY_TBL, spectrum.GetSampleCount(), 2,(char**)ttype, (char**)tform, NULL, NULL, &status) )
    {
    	return false;
    }

    if( fits_write_col(fptr, TDOUBLE, 1, 1, 1, spectrum.GetSampleCount(), spectrum.GetSpectralAxis().GetSamples(), &status ) )
    {
    	return false;
    }

    if( fits_write_col(fptr, TDOUBLE, 2, 1, 1, spectrum.GetSampleCount(), spectrum.GetFluxAxis().GetSamples(), &status ) )
    {
    	return false;
    }

    if( fits_close_file(fptr, &status) )
    {
    	return false;
    }


    return retv;
}
