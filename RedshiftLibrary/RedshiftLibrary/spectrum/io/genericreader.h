#ifndef _REDSHIFT_SPECTRUM_IO_GENERICREADER_
#define _REDSHIFT_SPECTRUM_IO_GENERICREADER_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/io/reader.h>


namespace NSEpic
{

class CSpectrumIOGenericReader : public CSpectrumIOReader
{

public:

    CSpectrumIOGenericReader();
    virtual ~CSpectrumIOGenericReader();

    static Bool CanRead( const char* filePath );
    virtual Bool Read( const char* filePath, std::shared_ptr<CSpectrum> spectrum );

private:

};


}

#endif
