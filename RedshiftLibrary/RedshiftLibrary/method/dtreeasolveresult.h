#ifndef _REDSHIFT_OPERATOR_DTREEASOLVERESULT_
#define _REDSHIFT_OPERATOR_DTREEASOLVERESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/ray/catalog.h>

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

public:

    CDTreeASolveResult();
    virtual ~CDTreeASolveResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;

};


}

#endif
