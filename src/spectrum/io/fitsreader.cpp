#include <epic/redshift/spectrum/io/fitsreader.h>

#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/axis.h>


using namespace NSEpic;
using namespace std;

CSpectrumIOFitsReader::CSpectrumIOFitsReader()
{

}

CSpectrumIOFitsReader::~CSpectrumIOFitsReader()
{

}

Bool CSpectrumIOFitsReader::Read2( fitsfile* fptr, CSpectrum& spectrum )
{
    Int32 status = 0;
    Int32 length=0;
    Int32 hdunum=0;
    Int32 hdutype=0;
    Int32 nullval = 0;
    Int32 anynul = 0;

    // Move to second hdu
    if(fits_movabs_hdu(fptr, 2, &hdutype, &status))
        return false;

    // Check type
    if (hdutype!=BINARY_TBL)
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
    Int32 length=0;
    Int32 nfound;
    Int32 hdunum=0;
    Int32 hdutype=0;

    // Move to first hdu
    if(fits_movabs_hdu(fptr, 1, &hdutype, &status))
        return false;

    // Check type
    if (hdutype!=IMAGE_HDU)
        return false;

    // read NAXIS k/w
    if(fits_read_key(fptr, TINT, "NAXIS", &naxis, NULL, &status))
        return false;

    // check axis number
    if( naxis > 2)
        return false;

    // read NAXIS lengths
    if(fits_read_keys_lng(fptr, "NAXIS", 1, naxis, naxiss, &nfound, &status))
        return -1;

    // spectrum length found
    length=naxiss[0];

    // Read data
    CSpectrumAxis& spcFluxAxis = spectrum.GetFluxAxis();

    Int32 nullval = 0;
    Int32 anynul = 0;
    spcFluxAxis.SetSize( length );
    if(fits_read_img(fptr, TDOUBLE, 1, length, &nullval, spcFluxAxis.GetSamples(), &anynul, &status))
        return false;

    // read keywords
    float crpix1, crval1, cdelt1;
    if(fits_read_key(fptr, TFLOAT, "CRPIX1", &crpix1, NULL, &status))
        return false;

    if(fits_read_key(fptr, TFLOAT, "CRVAL1", &crval1, NULL, &status))
        return false;

    if(fits_read_key(fptr, TFLOAT, "CDELT1", &cdelt1, NULL, &status))
        return false;

    // wavelength array
    CSpectrumAxis& spcSpectralAxis = spectrum.GetSpectralAxis();

    spcSpectralAxis.SetSize( length );
    double wave_value = crval1 - cdelt1 * (crpix1-1);
    spcSpectralAxis[0] = wave_value;

    // set wavelength
    for(Int32 i=1; i<length; i++)
    {
        wave_value += cdelt1;
        spcSpectralAxis[i]=wave_value;
    }

    fits_close_file(fptr, &status);
    return true;
}

Bool CSpectrumIOFitsReader::Read( const char* filePath, CSpectrum& spectrum )
{
    fitsfile *fptr = NULL;
    Int32 status = 0;
    Int32 hdunum=0;

    Bool retv = true;
    // open the fits file
    if( ! fits_open_file(&fptr, filePath, READONLY, &status))
    {
        // recover hdu number
        if(!fits_get_num_hdus(fptr, &hdunum, &status))
        {
            // check hdus
            if (hdunum==1)
            {
                if( !Read1( fptr, spectrum ) )
                    retv = false;
            }
            else if(hdunum == 2 )
            {
                if( !Read2( fptr, spectrum ) )
                    retv = false;
            }
            else
            {
                retv = false;
            }
        }
    }
    else
    {
        return false;
    }

    fits_close_file(fptr, &status);

    return retv;
}
