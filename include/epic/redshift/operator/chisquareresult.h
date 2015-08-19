#ifndef _REDSHIFT_OPERATOR_CHISQUARERESULT_
#define _REDSHIFT_OPERATOR_CHISQUARERESULT_

#include <epic/redshift/operator/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/operator/operator.h>

namespace NSEpic
{

class CChisquareResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CChisquareResult )

public:

    CChisquareResult();
    virtual ~CChisquareResult();

    Void Save( const COperatorResultStore& store, std::ostream& stream ) const;
    Void SaveLine( const COperatorResultStore& store, std::ostream& stream ) const;
    Void Load( std::istream& stream );

    TFloat64List            Redshifts;
    TFloat64List            ChiSquare;
    TFloat64List            Overlap;
    TFloat64List            Extrema;
    COperator::TStatusList  Status;

};


}

#endif
