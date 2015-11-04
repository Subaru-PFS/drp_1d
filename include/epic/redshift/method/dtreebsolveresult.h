#ifndef _REDSHIFT_OPERATOR_DTREEBSOLVERESULT_
#define _REDSHIFT_OPERATOR_DTREEBSOLVERESULT_

#include <epic/redshift/operator/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/catalog.h>

#include <vector>

namespace NSEpic
{

class CProcessFlowContext;

class CDTreeBSolveResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CDTreeBSolveResult )

public:

    CDTreeBSolveResult();
    virtual ~CDTreeBSolveResult();

    Void Save( const COperatorResultStore& store, std::ostream& stream ) const;
    Void SaveLine( const COperatorResultStore& store, std::ostream& stream ) const;

    Bool GetBestRedshift( const COperatorResultStore& store, Float64& redshift, Float64& merit ) const;
    Bool GetBestRedshiftChi2( const COperatorResultStore& store, std::string scopeStr, Float64 targetz, Float64& redshift, Float64& merit, std::string& tplName ) const;
};


}

#endif
