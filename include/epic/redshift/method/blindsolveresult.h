#ifndef _REDSHIFT_OPERATOR_BLINDSOLVERESULT_
#define _REDSHIFT_OPERATOR_BLINDSOLVERESULT_

#include <epic/redshift/operator/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/catalog.h>

#include <vector>

namespace NSEpic
{

class CProcessFlowContext;

class CBlindSolveResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CBlindSolveResult )

public:

    CBlindSolveResult();
    virtual ~CBlindSolveResult();

    Void Save( const COperatorResultStore& store, std::ostream& stream ) const;
    Void SaveLine( const COperatorResultStore& store, std::ostream& stream ) const;
    Bool GetBestFitResult( const COperatorResultStore& store, Float64& redshift, Float64& merit, std::string& tplName ) const;
};


}

#endif
