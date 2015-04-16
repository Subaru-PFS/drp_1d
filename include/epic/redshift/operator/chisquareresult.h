#ifndef _REDSHIFT_OPERATOR_CHISQUARERESULT_
#define _REDSHIFT_OPERATOR_CHISQUARERESULT_

#include <epic/redshift/operator/result.h>
#include <epic/core/common/datatypes.h>

namespace NSEpic
{

class CChisquareResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CChisquareResult )

public:

    CChisquareResult();
    virtual ~CChisquareResult();

    Void Save( std::ostream& stream ) const;

    TFloat64List    Redshifts;
    TFloat64List    ChiSquare;
    TFloat64List    Overlap;

};


}

#endif
