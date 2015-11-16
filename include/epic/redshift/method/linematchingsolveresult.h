#ifndef _REDSHIFT_OPERATOR_LINEMATCHINGSOLVERESULT_
#define _REDSHIFT_OPERATOR_LINEMATCHINGSOLVERESULT_

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
class CLineMatchingSolveResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CLineMatchingSolveResult )

public:

    CLineMatchingSolveResult();
    virtual ~CLineMatchingSolveResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    Bool GetBestResult( const CDataStore& store, Float64& redshift, Float64& merit ) const;

};


}

#endif
