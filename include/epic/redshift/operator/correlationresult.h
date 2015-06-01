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

    Void Save( const COperatorResultStore& store, std::ostream& stream ) const;
    Void Load( std::istream& stream );

    TFloat64List    Redshifts;
    TFloat64List    Correlation;
    TFloat64List    Overlap;
    COperator::TStatusList  Status;

};


}

#endif
