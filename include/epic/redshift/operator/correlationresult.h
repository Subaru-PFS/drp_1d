#ifndef _REDSHIFT_OPERATOR_CORRELATIONRESULT_
#define _REDSHIFT_OPERATOR_CORRELATIONRESULT_

#include <epic/redshift/operator/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/operator/operator.h>

namespace NSEpic
{

class CCorrelationResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CCorrelationResult )

public:

    CCorrelationResult();
    virtual ~CCorrelationResult();

    Void Save( std::ostream& stream ) const;

    TFloat64List    Redshifts;
    TFloat64List    Correlation;
    TFloat64List    Overlap;
    COperator::TStatusList  Status;

};


}

#endif
