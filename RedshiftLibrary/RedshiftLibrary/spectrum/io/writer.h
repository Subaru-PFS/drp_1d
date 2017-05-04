#ifndef _REDSHIFT_SPECTRUM_IO_WRITER_
#define _REDSHIFT_SPECTRUM_IO_WRITER_

#include <RedshiftLibrary/common/datatypes.h>

#include <vector>

namespace NSEpic
{

class CSpectrum;

class CSpectrumIOWriter
{

public:

    CSpectrumIOWriter();
    virtual ~CSpectrumIOWriter();

    virtual Bool Write( const char* path, CSpectrum& s ) = 0;

private:

};


}

#endif
