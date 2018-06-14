#ifndef _REDSHIFT_SPECTRUM_IO_FITSREADER_
#define _REDSHIFT_SPECTRUM_IO_FITSREADER_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/io/reader.h>

#include <vector>
#include <fitsio.h>

namespace NSEpic
{

class CSpectrumIOFitsReader : public CSpectrumIOReader
{

public:

    CSpectrumIOFitsReader();
    ~CSpectrumIOFitsReader();

    virtual Void Read( const char* filePath, std::shared_ptr<CSpectrum> s );

private:
    Bool Read1( fitsfile* fptr, std::shared_ptr<CSpectrum> spectrum );
    Bool Read2( fitsfile* fptr, std::shared_ptr<CSpectrum> spectrum );
};


}

#endif
