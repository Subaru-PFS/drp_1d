#ifndef _REDSHIFT_METHOD_LINEMATCHINGSOLVERESULT_
#define _REDSHIFT_METHOD_LINEMATCHINGSOLVERESULT_

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
  class CLineMatchingSolveResult : public COperatorResult
{

public:

    CLineMatchingSolveResult();
    virtual ~CLineMatchingSolveResult();

    void Save(std::ostream& stream ) const;
    void SaveLine(std::ostream& stream ) const;

    Bool GetBestResult( const CDataStore& store, Float64& redshift, Float64& merit ) const;

};


}

#endif
