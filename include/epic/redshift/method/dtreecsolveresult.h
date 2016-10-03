#ifndef _REDSHIFT_OPERATOR_DTREECSOLVERESULT_
#define _REDSHIFT_OPERATOR_DTREECSOLVERESULT_

#include <epic/redshift/processflow/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/catalog.h>
#include <epic/redshift/processflow/result.h>

#include <vector>

namespace NSEpic
{

class CProcessFlowContext;

/**
 * \ingroup Redshift
 */
class CDTreeCSolveResult : public COperatorResult
{

public:

    CDTreeCSolveResult();
    virtual ~CDTreeCSolveResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;

    Bool GetBestRedshift(const CDataStore& store, Float64& redshift, Float64& merit, std::string &tplName, std::string& dtreepath) const;
    std::string GetBestContinuumTplNameAtRedshift( const CDataStore& store, Float64 z) const;


};


}

#endif
