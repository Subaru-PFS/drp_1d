#ifndef _REDSHIFT_METHOD_LINEMATCHINGSOLVE2RESULT_
#define _REDSHIFT_METHOD_LINEMATCHINGSOLVE2RESULT_

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
class CLineMatching2SolveResult : public COperatorResult
{

public:

    CLineMatching2SolveResult();
    virtual ~CLineMatching2SolveResult();

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }

    Bool GetBestResult( const CDataStore& store, Float64& redshift, Float64& merit ) const;

};


}

#endif
