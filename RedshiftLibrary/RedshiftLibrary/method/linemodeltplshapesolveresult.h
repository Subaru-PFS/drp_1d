#ifndef _REDSHIFT_OPERATOR_LINEMODELTPLSHAPESOLVERESULT_
#define _REDSHIFT_OPERATOR_LINEMODELTPLSHAPESOLVERESULT_

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
class CLineModelTplshapeSolveResult : public COperatorResult
{

public:

    enum EType
    {
         nType_raw = 1,
         nType_continuumOnly = 2,
         nType_noContinuum = 3,
         nType_all = 4,
    };

    CLineModelTplshapeSolveResult();
    virtual ~CLineModelTplshapeSolveResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    Bool GetBestRedshift( const CDataStore& store, Float64& redshift, Float64& merit, std::string& tplName) const;
    Bool GetBestRedshiftPerTemplateString( const CDataStore& store, std::string& output ) const;


};


}

#endif

