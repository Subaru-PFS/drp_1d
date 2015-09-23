#ifndef _REDSHIFT_OPERATOR_LINEMODELSOLVERESULT_
#define _REDSHIFT_OPERATOR_LINEMODELSOLVERESULT_

#include <epic/redshift/operator/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/catalog.h>

#include <vector>

namespace NSEpic
{

class CProcessFlowContext;

class CLineModelSolveResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CLineModelSolveResult )

public:

    enum EType
    {
         nType_raw = 1,
         nType_continuumOnly = 2,
         nType_noContinuum = 3,
         nType_all = 4,
    };

    CLineModelSolveResult();
    virtual ~CLineModelSolveResult();

    Void Save( const COperatorResultStore& store, std::ostream& stream ) const;
    Void SaveLine( const COperatorResultStore& store, std::ostream& stream ) const;
    Bool GetBestRedshift( const COperatorResultStore& store, Float64& redshift, Float64& merit ) const;
    Bool GetBestRedshiftLogArea( const COperatorResultStore& store, Float64& redshift, Float64& merit ) const;

};


}

#endif

