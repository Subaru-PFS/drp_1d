#ifndef _REDSHIFT_OPERATOR_DTREE7SOLVERESULT_
#define _REDSHIFT_OPERATOR_DTREE7SOLVERESULT_

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
class CDTree7SolveResult : public COperatorResult
{

public:

    CDTree7SolveResult();
    virtual ~CDTree7SolveResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;

};


}

#endif
