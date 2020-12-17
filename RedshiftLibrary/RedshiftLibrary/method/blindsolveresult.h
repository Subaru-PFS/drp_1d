#ifndef _REDSHIFT_METHOD_BLINDSOLVERESULT_
#define _REDSHIFT_METHOD_BLINDSOLVERESULT_

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
class CBlindSolveResult : public COperatorResult
{

public:

    CBlindSolveResult();
    virtual ~CBlindSolveResult();

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    Bool GetBestFitResult( const CDataStore& store, Float64& redshift, Float64& merit, std::string& tplName ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }
};


}

#endif
