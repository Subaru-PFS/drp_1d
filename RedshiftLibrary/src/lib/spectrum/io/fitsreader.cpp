#include <RedshiftLibrary/spectrum/io/fitsreader.h>
#include <RedshiftLibrary/log/log.h>

#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/axis.h>
#include <limits>


using namespace NSEpic;
using namespace std;

/**
 * Empty constructor.
 */
CSpectrumIOFitsReader::CSpectrumIOFitsReader()
{

}

/**
 * Empty destructor.
 */
CSpectrumIOFitsReader::~CSpectrumIOFitsReader()
{

}

/**
 * If the fits file has 2 HDUs, if the 2nd HDU is binary table, if this table only has 2 columns, then continue - otherwise return false.
 * Call the spectrum methods to get the flux and spectral axii, if the column data can be read as samples, return true, otherwise return false.
 */
Bool CSpectrumIOFitsReader::Read2( fitsfile* fptr, CSpectrum& spectrum )
{
    Int32 status = 0;
    Int32 length = 0;
    Int32 hdunum = 2;
    Int32 hdutype = 0;
    Float64 nullval = std::numeric_limits<Float64>::quiet_NaN(); //see https://heasarc.gsfc.nasa.gov/fitsio/c/c_user/node80.html
    Int32 anynul = 0;

    // Move to second hdu
    if( fits_movabs_hdu( fptr, hdunum, &hdutype, &status ) )
        return false;

    // Check type
    if ( hdutype!=BINARY_TBL )
        return false;

    Int32 nbCols = 0;
    long nbRows = 0;
    if( fits_get_num_cols( fptr, &nbCols, &status ) )
        return false;

    if( fits_get_num_rows( fptr, &nbRows, &status ) )
        return false;

    // only two columns must be present: WAVE and FLUX
    if( nbCols != 2 )
        return false;

    length = (Int32) nbRows;
    CSpectrumAxis& spcFluxAxis = spectrum.GetFluxAxis();
    CSpectrumAxis& spcSpectralAxis = spectrum.GetSpectralAxis();

    spcFluxAxis.SetSize( length );
    spcSpectralAxis.SetSize( length );

    if( fits_read_col( fptr, TDOUBLE, 1, 1, 1, length, &nullval, spcSpectralAxis.GetSamples(), &anynul, &status ) )
        return false;

    if( fits_read_col( fptr, TDOUBLE, 2, 1, 1, length, &nullval, spcFluxAxis.GetSamples(), &anynul, &status ) )
        return false;

    return true;
}

/**
 * Read FITS file formated with the following rules:
 *
 * - One Primary HDU of Image type with the following keys in it's header:
 *
 * - NAXIS: 2
 * - NAXIS1: Number of pixel in the flux axis
 * - NAXIS2: 1
 *
 * - CRPIX1: -1
 * - CRVAL1: Start lambda range
 * - CDELT1: lambda range between two consecutives samples
 */
Bool CSpectrumIOFitsReader::Read1( fitsfile* fptr, CSpectrum& spectrum )
{
    Int32 status = 0;
    long naxiss[2];
    Int32 naxis;
    Int32 length = 0;
    Int32 nfound;
    Int32 hdunum = 0;
    Int32 hdutype = 0;

    // Move to first hdu
    if( fits_movabs_hdu( fptr, 1, &hdutype, &status ) )
        return false;

    // Check type
    if ( hdutype!=IMAGE_HDU )
        return false;

    // read NAXIS k/w
    if( fits_read_key( fptr, TINT, "NAXIS", &naxis, NULL, &status ) )
        return false;

    // check axis number
    if( naxis>2 )
        return false;

    // read NAXIS lengths
    if( fits_read_keys_lng( fptr, "NAXIS", 1, naxis, naxiss, &nfound, &status ) )
        return -1;

    // spectrum length found
    length = naxiss[0];

    // Read data
    CSpectrumAxis& spcFluxAxis = spectrum.GetFluxAxis();

    Float64 nullval = std::numeric_limits<Float64>::quiet_NaN();
    Int32 anynul = 0;
    spcFluxAxis.SetSize( length );
    if( fits_read_img( fptr, TDOUBLE, 1, length, &nullval, spcFluxAxis.GetSamples(), &anynul, &status ) )
        return false;

    Log.LogDebug("    CSpectrumIOFitsReader: loaded flux values n samples=%d", length);
    Log.LogDebug("    CSpectrumIOFitsReader: loaded flux values first sample=%f", spcFluxAxis.GetSamples()[0]);

    // read keywords
    float crpix1, crval1, cdelt1;
    if( fits_read_key( fptr, TFLOAT, "CRPIX1", &crpix1, NULL, &status ) )
        return false;

    if( fits_read_key( fptr, TFLOAT, "CRVAL1", &crval1, NULL, &status ) )
        return false;

    if( fits_read_key( fptr, TFLOAT, "CDELT1", &cdelt1, NULL, &status ) )
        return false;

    Log.LogDebug("    CSpectrumIOFitsReader: loaded CRPIX1=%f", crpix1);
    Log.LogDebug("    CSpectrumIOFitsReader: loaded CRVAL1=%f", crval1);
    Log.LogDebug("    CSpectrumIOFitsReader: loaded CDELT1=%f", cdelt1);

    // wavelength array
    CSpectrumAxis& spcSpectralAxis = spectrum.GetSpectralAxis();

    spcSpectralAxis.SetSize( length );
    double wave_value = crval1 - cdelt1 * (crpix1-1);
    //double wave_value = crval1 + cdelt1 + crpix1*cdelt1; //Warning, modified from: wave_value = crval1 - cdelt1 * (crpix1-1);, to be further checked...
    spcSpectralAxis[0] = wave_value;

    // set wavelength
    for(Int32 i=1; i<length; i++)
    {
        wave_value += cdelt1;
        spcSpectralAxis[i]=wave_value;
    }

    Log.LogDebug("    CSpectrumIOFitsReader: loaded spectral values first sample=%f", spcSpectralAxis.GetSamples()[0]);

    return true;
}

/**
 * Attempts to read the file as a fits file containing 2 HDUs, using the Read1 and Read2 methods to read each HDU, respectively. If all calls return true, this method returns true as well - and returns false otherwise.
 */
void CSpectrumIOFitsReader::Read( const char* filePath, CSpectrum& spectrum )
{
    fitsfile *fptr = NULL;
    Int32 status = 0;
    Int32 hdunum=0;

    // open the fits file
    if( !fits_open_file( &fptr, filePath, READONLY, &status ) )
    {
        // recover hdu number
        if( !fits_get_num_hdus( fptr, &hdunum, &status ) )
        {
            // check hdus
            if ( hdunum==1 )
            {
                Log.LogDebug("    CSpectrumIOFitsReader: Read1");
                if( !Read1( fptr, spectrum ) )
		  {
		    fits_close_file( fptr, &status );
		    Log.LogError("error in Read1 : %s", filePath);
		    throw runtime_error("error in Read1");
		  }
            }
            else if( hdunum == 2 )
            {
	      Log.LogDebug("    CSpectrumIOFitsReader: Read2");
                if( !Read2( fptr, spectrum ) )
		  {
		    fits_close_file( fptr, &status );
		    Log.LogError("error in Read2 : %s", filePath);
		    throw runtime_error("error in Read2 ");
		  }
            }
            else
            {
	      fits_close_file( fptr, &status );
	      Log.LogError("bad hdu count : %s", filePath);
	      throw runtime_error("bad hdu count");
            }
        } else {
	  fits_close_file( fptr, &status );
	  Log.LogError("bad hdu count : %s", filePath);
	  throw runtime_error("bad hdu count");
	}
    }
    else
    {
      fits_close_file( fptr, &status );
      Log.LogError("error opening fits file : %s", filePath);
      throw runtime_error("error opening fits file");
    }
    fits_close_file( fptr, &status );
}
