#ifndef _REDSHIFT_METHOD_DTREE7SOLVERESULT_
#define _REDSHIFT_METHOD_DTREE7SOLVERESULT_

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
class CDTree7SolveResult : public COperatorResult
{

public:

    CDTree7SolveResult();
    virtual ~CDTree7SolveResult();

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }

};


}

#endif
