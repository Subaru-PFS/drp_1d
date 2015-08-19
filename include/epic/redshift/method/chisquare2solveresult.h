#ifndef _REDSHIFT_OPERATOR_CHISQUARE2SOLVERESULT_
#define _REDSHIFT_OPERATOR_CHISQUARE2SOLVERESULT_

#include <epic/redshift/operator/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/catalog.h>

#include <vector>

namespace NSEpic
{

class CProcessFlowContext;

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

    Void Save( const COperatorResultStore& store, std::ostream& stream ) const;
    Void SaveLine( const COperatorResultStore& store, std::ostream& stream ) const;
    Bool GetBestRedshift( const COperatorResultStore& store, Float64& redshift, Float64& merit, std::string& tplName ) const;
    Bool GetBestRedshiftPerTemplateString( const COperatorResultStore& store, std::string& output ) const;

    Int32 m_type;

};


}

#endif


