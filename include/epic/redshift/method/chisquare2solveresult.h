#ifndef _REDSHIFT_OPERATOR_CHISQUARE2SOLVERESULT_
#define _REDSHIFT_OPERATOR_CHISQUARE2SOLVERESULT_

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
class CChisquare2SolveResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CChisquare2SolveResult )

public:

    enum EType
    {
             nType_raw = 1,
             nType_continuumOnly = 2,
             nType_noContinuum = 3,
             nType_all = 4,
    };

    CChisquare2SolveResult();
    virtual ~CChisquare2SolveResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    Bool GetBestRedshift( const CDataStore& store, Float64& redshift, Float64& merit, std::string& tplName ) const;
    Bool GetBestRedshiftPerTemplateString( const CDataStore& store, std::string& output ) const;

    Int32 m_type;

};


}

#endif


