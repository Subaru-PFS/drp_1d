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

    Void Save( std::ostream& stream ) const;
    Bool GetBestCorrelationResult( const CProcessFlowContext& ctx, Float64& redshift, Float64& merit, std::string& tplName ) const;
    Bool GetBestCorrelationPeakResult( const CProcessFlowContext& ctx, Float64& redshift, Float64& merit, std::string& tplName ) const;

};


}

#endif