#ifndef _REDSHIFT_OPERATOR_CORRELATIONSOLVERESULT_
#define _REDSHIFT_OPERATOR_CORRELATIONSOLVERESULT_

#include <epic/redshift/operator/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/catalog.h>

#include <vector>

namespace NSEpic
{

class CProcessFlowContext;

class CCorrelationSolveResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CCorrelationSolveResult )

public:

    CCorrelationSolveResult();
    virtual ~CCorrelationSolveResult();

    Void Save( const COperatorResultStore& store, std::ostream& stream ) const;
    Void SaveLine( const COperatorResultStore& store, std::ostream& stream ) const;
    Bool GetBestCorrelationResult( const COperatorResultStore& store, Float64& redshift, Float64& merit, std::string& tplName ) const;

};


}

#endif
