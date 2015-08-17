#ifndef _REDSHIFT_OPERATOR_DTREE7SOLVERESULT_
#define _REDSHIFT_OPERATOR_DTREE7SOLVERESULT_

#include <epic/redshift/operator/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/catalog.h>

#include <vector>

namespace NSEpic
{

class CProcessFlowContext;

class CDTree7SolveResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CDTree7SolveResult )

public:

    CDTree7SolveResult();
    virtual ~CDTree7SolveResult();

    Void Save( const COperatorResultStore& store, std::ostream& stream ) const;
    Void SaveLine( const COperatorResultStore& store, std::ostream& stream ) const;

};


}

#endif
