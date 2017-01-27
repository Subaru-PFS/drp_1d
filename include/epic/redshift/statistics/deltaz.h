#ifndef _REDSHIFT_RAY_DELTAZ_
#define _REDSHIFT_RAY_DELTAZ_

#include <epic/core/common/datatypes.h>

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

    Float64 Compute(TFloat64List merits, TFloat64List redshifts, Float64 redshift);

private:


};


}

#endif

