#ifndef _REDSHIFT_RAY_DELTAZ_
#define _REDSHIFT_RAY_DELTAZ_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>


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
    Float64 GetDeltaz(TFloat64List redshifts, TFloat64List pdf, Float64 z, Int32 gslfit = 0);
    Int32 Compute3ddl(TFloat64List merits, TFloat64List redshifts, Float64 redshift, TFloat64Range redshiftRange, Float64 &deltaz);
    Int32 Compute(TFloat64List merits, TFloat64List redshifts, Float64 redshift, TFloat64Range redshiftRange, Float64& sigma);
private:


};


}

#endif

