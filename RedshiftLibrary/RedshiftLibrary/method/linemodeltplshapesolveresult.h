#ifndef _REDSHIFT_OPERATOR_LINEMODELTPLSHAPESOLVERESULT_
#define _REDSHIFT_OPERATOR_LINEMODELTPLSHAPESOLVERESULT_

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

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }

    Bool GetBestRedshift( const CDataStore& store, Float64& redshift, Float64& merit, std::string& tplName) const;
    Bool GetBestRedshiftPerTemplateString( const CDataStore& store, std::string& output ) const;


};


}

#endif

