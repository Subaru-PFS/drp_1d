#ifndef _REDSHIFT_SPECTRUM_IO_GENERICREADER_
#define _REDSHIFT_SPECTRUM_IO_GENERICREADER_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/io/reader.h>


namespace NSEpic
{

class CSpectrumIOGenericReader : public CSpectrumIOReader
{

public:

    CSpectrumIOGenericReader();
    ~CSpectrumIOGenericReader();

    static Bool CanRead( const char* filePath );
    Bool Read( const char* filePath, CSpectrum& s );

private:

};


}

#endif
