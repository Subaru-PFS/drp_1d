#ifndef _REDSHIFT_OPERATOR_PEAKDETECTIONRESULT_
#define _REDSHIFT_OPERATOR_PEAKDETECTIONRESULT_

#include "RedshiftLibrary/processflow/result.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"

#include <vector>

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class CPeakDetectionResult : public COperatorResult
{

public:

    CPeakDetectionResult();
    virtual ~CPeakDetectionResult();

    TInt32RangeList PeakList;
    TInt32RangeList EnlargedPeakList;

};


}

#endif
