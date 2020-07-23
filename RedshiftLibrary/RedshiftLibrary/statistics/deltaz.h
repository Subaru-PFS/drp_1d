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
    Float64 GetDeltaz(const TFloat64List & redshifts, const TFloat64List & pdf, const Float64 z, const Int32 gslfit=0);
    Int32 Compute3ddl(const TFloat64List &merits, const TFloat64List &redshifts, const Float64 redshift, const TFloat64Range & redshiftRange, Float64& deltaz);
    Int32 Compute(const TFloat64List &merits, const TFloat64List &redshifts, const Float64 redshift, const TFloat64Range &redshiftRange, Float64& sigma);
private:


};


}

#endif

