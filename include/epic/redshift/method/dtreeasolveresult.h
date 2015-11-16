#ifndef _REDSHIFT_OPERATOR_DTREEASOLVERESULT_
#define _REDSHIFT_OPERATOR_DTREEASOLVERESULT_

#include <epic/redshift/processflow/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/catalog.h>

#include <vector>

namespace NSEpic
{

class CProcessFlowContext;
class CDataStore;

/**
 * \ingroup Redshift
 */
class CDTreeASolveResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CDTreeASolveResult )

public:

    CDTreeASolveResult();
    virtual ~CDTreeASolveResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;

};


}

#endif
