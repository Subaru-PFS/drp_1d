#ifndef _REDSHIFT_OPERATOR_BLINDSOLVERESULT_
#define _REDSHIFT_OPERATOR_BLINDSOLVERESULT_

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
class CBlindSolveResult : public COperatorResult
{

public:

    CBlindSolveResult();
    virtual ~CBlindSolveResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    Bool GetBestFitResult( const CDataStore& store, Float64& redshift, Float64& merit, std::string& tplName ) const;
};


}

#endif
