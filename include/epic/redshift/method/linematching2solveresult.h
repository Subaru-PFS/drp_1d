#ifndef _REDSHIFT_OPERATOR_LINEMATCHINGSOLVE2RESULT_
#define _REDSHIFT_OPERATOR_LINEMATCHINGSOLVE2RESULT_

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
class CLineMatching2SolveResult : public COperatorResult
{

public:

    CLineMatching2SolveResult();
    virtual ~CLineMatching2SolveResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    Bool GetBestResult( const CDataStore& store, Float64& redshift, Float64& merit ) const;

};


}

#endif
