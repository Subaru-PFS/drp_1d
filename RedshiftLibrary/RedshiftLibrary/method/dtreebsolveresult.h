#ifndef _REDSHIFT_OPERATOR_DTREEBSOLVERESULT_
#define _REDSHIFT_OPERATOR_DTREEBSOLVERESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/processflow/result.h>

#include <vector>

namespace NSEpic
{

class CProcessFlowContext;

/**
 * \ingroup Redshift
 */
class CDTreeBSolveResult : public COperatorResult
{

public:

    CDTreeBSolveResult();
    virtual ~CDTreeBSolveResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;

    Bool GetBestRedshift( const CDataStore& store, Float64& redshift, Float64& merit, std::string& dtreepath) const;
};


}

#endif
