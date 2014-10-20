#ifndef _REDSHIFT_SPECTRUM_IO_READER_
#define _REDSHIFT_SPECTRUM_IO_READER_

#include <epic/core/common/datatypes.h>

#include <vector>

namespace NSEpic
{

class CSpectrum;

class CSpectrumIOReader
{

public:

    CSpectrumIOReader();
    virtual ~CSpectrumIOReader();

    virtual Bool Read( const char* fluxPath, CSpectrum& s ) = 0;

private:

};


}

#endif
