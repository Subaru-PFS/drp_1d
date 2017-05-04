#ifndef _REDSHIFT_RAY_DELTAZ_
#define _REDSHIFT_RAY_DELTAZ_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>


#include <string>

namespace NSEpic
{

/**
 * \ingroup Redshift
 * Deltaz
 */
class CDeltaz
{

public:


    CDeltaz();
    ~CDeltaz();

    Int32 Compute(TFloat64List merits, TFloat64List redshifts, Float64 redshift, TFloat64Range redshiftRange, Float64 &deltaz);
    Int32 Compute3ddl(TFloat64List merits, TFloat64List redshifts, Float64 redshift, TFloat64Range redshiftRange, Float64 &deltaz);

private:


};


}

#endif

