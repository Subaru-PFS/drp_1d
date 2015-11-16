#ifndef _REDSHIFT_OPERATOR_FULLSOLVERESULT_
#define _REDSHIFT_OPERATOR_FULLSOLVERESULT_

#include <epic/redshift/processflow/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/catalog.h>

#include <vector>

namespace NSEpic
{

class CProcessFlowContext;

/**
 * \ingroup Redshift
 */
class CFullSolveResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CFullSolveResult )

public:

    CFullSolveResult();
    virtual ~CFullSolveResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    Bool GetBestCorrelationResult( const CDataStore& store, Float64& redshift, Float64& merit, std::string& tplName ) const;

};


}

#endif
