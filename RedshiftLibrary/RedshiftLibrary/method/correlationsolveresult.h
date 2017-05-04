#ifndef _REDSHIFT_OPERATOR_CORRELATIONSOLVERESULT_
#define _REDSHIFT_OPERATOR_CORRELATIONSOLVERESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/ray/catalog.h>

#include <vector>

namespace NSEpic
{

class CProcessFlowContext;

/**
 * \ingroup Redshift
 */
class CCorrelationSolveResult : public COperatorResult
{

public:

    CCorrelationSolveResult();
    virtual ~CCorrelationSolveResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    Bool GetBestCorrelationResult( const CDataStore& store, Float64& redshift, Float64& merit, std::string& tplName ) const;

};


}

#endif
