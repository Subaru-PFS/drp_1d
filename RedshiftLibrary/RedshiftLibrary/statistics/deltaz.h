#ifndef _REDSHIFT_STATISTICS_DELTAZ_
#define _REDSHIFT_STATISTICS_DELTAZ_

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
    Int32 GetIndices(const TFloat64List & redshifts, const Float64 redshift, const Int32 HalfNbSamples, 
                          Int32 & iz, Int32 & izmin, Int32 & izmax );
    Int32 GetRangeIndices(const TFloat64List & redshifts, const Float64 redshift, const TFloat64Range & redshiftRange, 
                          Int32 & iz, Int32 & izmin, Int32 & izmax );
    Int32 Compute3ddl(const TFloat64List &merits, const TFloat64List &redshifts, const Int32 iz, const Int32 izmin, const Int32 izmax, Float64& sigma);
    Int32 Compute(const TFloat64List & merits, const TFloat64List & redshifts, const Int32 iz, const Int32 izmin, const Int32 izmax, Float64& sigma);

private:


};


}

#endif
