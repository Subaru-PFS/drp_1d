#ifndef _REDSHIFT_RAY_PDFZ_
#define _REDSHIFT_RAY_PDFZ_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>


#include <string>

namespace NSEpic
{

/**
 * \ingroup Redshift
 * Pdfz
 */
class CPdfz
{

public:


    CPdfz();
    ~CPdfz();

    Int32 Compute(TFloat64List merits, TFloat64List redshifts, Float64 cstLog, TFloat64List &logPdf);

private:


};


}

#endif

