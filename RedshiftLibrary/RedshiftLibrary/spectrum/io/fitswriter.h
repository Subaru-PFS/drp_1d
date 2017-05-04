#ifndef _REDSHIFT_SPECTRUM_IO_FITSWRITER_
#define _REDSHIFT_SPECTRUM_IO_FITSWRITER_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/io/writer.h>

#include <vector>
#include <fitsio.h>

namespace NSEpic
{

class CSpectrumIOFitsWriter : public CSpectrumIOWriter
{

public:

	CSpectrumIOFitsWriter();
    ~CSpectrumIOFitsWriter();

    Bool Write( const char* filePath, CSpectrum& s );

private:

};


}

#endif
