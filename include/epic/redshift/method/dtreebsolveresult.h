#ifndef _REDSHIFT_OPERATOR_DTREEBSOLVERESULT_
#define _REDSHIFT_OPERATOR_DTREEBSOLVERESULT_

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
class CDTreeBSolveResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CDTreeBSolveResult )

public:

    CDTreeBSolveResult();
    virtual ~CDTreeBSolveResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;

    Bool GetBestRedshift( const CDataStore& store, Float64& redshift, Float64& merit, std::string& dtreepath) const;

    Bool GetBestRedshiftChi2( const CDataStore& store, std::string scopeStr, Float64 targetz, Float64& redshift, Float64& merit, std::string& tplName ) const;
    TFloat64List GetBestRedshiftChi2List(const CDataStore& store, std::string scopeStr, Float64 &minmerit) const;
};


}

#endif
