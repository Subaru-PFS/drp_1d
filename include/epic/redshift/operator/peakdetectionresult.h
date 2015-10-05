#ifndef _REDSHIFT_OPERATOR_PEAKDETECTIONRESULT_
#define _REDSHIFT_OPERATOR_PEAKDETECTIONRESULT_

#include <epic/redshift/processflow/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>

#include <vector>

namespace NSEpic
{

class CPeakDetectionResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CPeakDetectionResult )

public:

    CPeakDetectionResult();
    virtual ~CPeakDetectionResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;

    TInt32RangeList PeakList;
    TInt32RangeList EnlargedPeakList;

};


}

#endif
