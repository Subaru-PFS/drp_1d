#ifndef _REDSHIFT_OPERATOR_CHISQUARESOLVERESULT_
#define _REDSHIFT_OPERATOR_CHISQUARESOLVERESULT_

#include <epic/redshift/operator/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/catalog.h>

#include <vector>

namespace NSEpic
{

class CProcessFlowContext;

class CChisquareSolveResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CChisquareSolveResult )

public:

    CChisquareSolveResult();
    virtual ~CChisquareSolveResult();

    Void Save( const COperatorResultStore& store, std::ostream& stream ) const;
    Void SaveLine( const COperatorResultStore& store, std::ostream& stream ) const;
    Bool GetBestRedshift( const COperatorResultStore& store, Float64& redshift, Float64& merit, std::string& tplName ) const;

};


}

#endif


