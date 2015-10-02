#ifndef _REDSHIFT_OPERATOR_DTREEASOLVERESULT_
#define _REDSHIFT_OPERATOR_DTREEASOLVERESULT_

#include <epic/redshift/processflow/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/catalog.h>

#include <vector>

namespace NSEpic
{

class CProcessFlowContext;

class CDTreeASolveResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CDTreeASolveResult )

public:

    CDTreeASolveResult();
    virtual ~CDTreeASolveResult();

    Void Save( const COperatorResultStore& store, std::ostream& stream ) const;
    Void SaveLine( const COperatorResultStore& store, std::ostream& stream ) const;

};


}

#endif
