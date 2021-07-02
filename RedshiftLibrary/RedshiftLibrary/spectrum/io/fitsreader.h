#ifndef _REDSHIFT_SPECTRUM_IO_FITSREADER_
#define _REDSHIFT_SPECTRUM_IO_FITSREADER_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/spectrum/io/reader.h"

#include <vector>
#include <fitsio.h>

namespace NSEpic
{

class CSpectrumIOFitsReader : public CSpectrumIOReader
{

public:

    CSpectrumIOFitsReader();
    ~CSpectrumIOFitsReader();

    virtual void Read( const char* filePath, CSpectrum& s );

private:
    Bool Read1( fitsfile* fptr, CSpectrum& spectrum );
    Bool Read2( fitsfile* fptr, CSpectrum& spectrum );
};


}

#endif
