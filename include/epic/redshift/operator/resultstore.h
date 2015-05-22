#ifndef _REDSHIFT_OPERATOR_RESULTSTORE_
#define _REDSHIFT_OPERATOR_RESULTSTORE_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/managedobject.h>
#include <epic/redshift/operator/result.h>

#include <vector>
#include <ostream>

namespace NSEpic
{

class CTemplate;

class COperatorResultStore : public CManagedObject
{

public:

    typedef std::map< std::string, CConstRef<COperatorResult> > TResultsMap;
    typedef std::map< std::string, TResultsMap>                 TPerTemplateResultsMap;


    COperatorResultStore();
    virtual ~COperatorResultStore();

    Void  StorePerTemplateResult( const CTemplate& t, const char* name, const COperatorResult& result );
    Void  StoreGlobalResult( const char* name, const COperatorResult& result );

    const COperatorResult*  GetPerTemplateResult( const CTemplate& t, const char* name ) const;
    TOperatorResultMap      GetPerTemplateResult( const char* name ) const;
    const COperatorResult*  GetGlobalResult( const char* name ) const;

protected:

    void StoreResult( TResultsMap& map, const char* name, const COperatorResult& result );

    TPerTemplateResultsMap          m_PerTemplateResults;
    TResultsMap                     m_GlobalResults;

};


}

#endif
