#ifndef _REDSHIFT_OPERATOR_CORRELATIONRESULT_
#define _REDSHIFT_OPERATOR_CORRELATIONRESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>

namespace NSEpic
{

class CCorrelationResult : public COperatorResult
{

public:

    CCorrelationResult();
    virtual ~CCorrelationResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    Void Load( std::istream& stream );

    TFloat64List    Redshifts;
    TFloat64List    Correlation;
    TFloat64List    Overlap;
    TFloat64List            Extrema;
    COperator::TStatusList  Status;

};


}

#endif
