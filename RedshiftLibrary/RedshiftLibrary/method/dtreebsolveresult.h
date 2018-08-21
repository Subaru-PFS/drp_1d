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

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }

    Bool GetBestRedshift( const CDataStore& store, Float64& redshift, Float64& merit, std::string& dtreepath) const;
};


}

#endif
