#ifndef _REDSHIFT_STATISTICS_ZPRIOR_
#define _REDSHIFT_STATISTICS_ZPRIOR_

#include "RedshiftLibrary/common/datatypes.h"

namespace NSEpic
{

/**
 * \ingroup Redshift
 * CZPrior
 */
class CZPrior
{

public:

    CZPrior();
    ~CZPrior();

    TFloat64List GetConstantLogZPrior(UInt32 nredshifts);
    TFloat64List GetStrongLinePresenceLogZPrior(const TBoolList & linePresence, const Float64 penalization_factor);
    TFloat64List GetNLinesSNRAboveCutLogZPrior(const TInt32List & nlinesAboveSNR, const Float64 penalization_factor);
    TFloat64List GetEuclidNhaLogZPrior(const TFloat64List & redshifts, const Float64 aCoeff);
    TFloat64List CombineLogZPrior(const TFloat64List & logprior1, const TFloat64List & logprior2);
};


}

#endif
