#ifndef _REDSHIFT_OPERATOR_CHISQUARESOLVERESULT_
#define _REDSHIFT_OPERATOR_CHISQUARESOLVERESULT_

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
class CChisquareSolveResult : public COperatorResult
{

public:

    CChisquareSolveResult();
    virtual ~CChisquareSolveResult();

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }

    Bool GetBestRedshift( const CDataStore& store, Float64& redshift, Float64& merit, std::string& tplName ) const;
    Bool GetBestRedshiftPerTemplateString( const CDataStore& store, std::string& output ) const;


};


}

#endif


