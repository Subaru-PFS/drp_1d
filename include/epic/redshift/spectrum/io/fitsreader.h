#ifndef _REDSHIFT_SPECTRUM_IO_FITSREADER_
#define _REDSHIFT_SPECTRUM_IO_FITSREADER_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/io/reader.h>

#include <vector>
#include <fitsio.h>

namespace NSEpic
{

class CSpectrumIOFitsReader : public CSpectrumIOReader
{

public:

    CSpectrumIOFitsReader();
    ~CSpectrumIOFitsReader();

    Bool Read( const char* filePath, CSpectrum& s );

private:

    Bool Read1( fitsfile* fptr, CSpectrum& spectrum );
    Bool Read2( fitsfile* fptr, CSpectrum& spectrum );
};


}

#endif
