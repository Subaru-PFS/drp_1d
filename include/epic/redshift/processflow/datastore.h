#ifndef _REDSHIFT_PROCESSFLOW_DATASTORE_
#define _REDSHIFT_PROCESSFLOW_DATASTORE_

#include <epic/core/common/datatypes.h>

#include <epic/redshift/processflow/parameterstore.h>
#include <epic/core/common/managedobject.h>
#include <epic/redshift/processflow/resultstore.h>

#include <boost/filesystem.hpp>
#include <vector>
#include <ostream>

namespace NSEpic
{

class CDataStore  : public CManagedObject
{

public:

    typedef std::vector<std::string>    TScopeStack;

    class CAutoScope {
        public:
            CAutoScope( CDataStore& store, const char* name );
            ~CAutoScope();
        private:
            CDataStore* m_Store;
    };

    CDataStore();
    virtual ~CDataStore();

    Void                PushScope( const char* name );
    Void                PopScope();

    std::string         GetCurrentScopeName() const;

    std::string         GetScope(CConstRef<COperatorResult>  result) const;

    // Wrapper functions
    Bool                            GetScopedParam( const char* name, TFloat64List& v, const TFloat64List& defaultValue = TFloat64List() ) const;
    Bool                            GetScopedParam( const char* name, TInt64List& v, const TInt64List& defaultValue = TInt64List() ) const;
    Bool                            GetScopedParam( const char* name, TBoolList& v, const TBoolList& defaultValue = TBoolList() ) const;
    Bool                            GetScopedParam( const char* name, Float64& v, Float64 defaultValue  = 0 ) const;
    Bool                            GetScopedParam( const char* name, Int64& v, Int64 defaultValue = 0 ) const;
    Bool                            GetScopedParam( const char* name, Bool& v, Bool defaultValue = true ) const;
    Bool                            GetScopedParam( const char* name, std::string& v, std::string defaultValue = "" ) const;

    Bool                            SetScopedParam( const char* name, const TFloat64List& v );
    Bool                            SetScopedParam( const char* name, const TInt64List& v );
    Bool                            SetScopedParam( const char* name, const TBoolList& v );
    Bool                            SetScopedParam( const char* name, Float64 v );
    Bool                            SetScopedParam( const char* name, Int64 v );
    Bool                            SetScopedParam( const char* name, Bool v );
    Bool                            SetScopedParam( const char* name, const std::string& v );

    Bool                            GetParam( const char* name, TFloat64List& v, const TFloat64List& defaultValue = TFloat64List() ) const;
    Bool                            GetParam( const char* name, TInt64List& v, const TInt64List& defaultValue = TInt64List() ) const;
    Bool                            GetParam( const char* name, TBoolList& v, const TBoolList& defaultValue = TBoolList() ) const;
    Bool                            GetParam( const char* name, Float64& v, Float64 defaultValue  = 0 ) const;
    Bool                            GetParam( const char* name, Int64& v, Int64 defaultValue = 0 ) const;
    Bool                            GetParam( const char* name, Bool& v, Bool defaultValue = true ) const;
    Bool                            GetParam( const char* name, std::string& v, std::string defaultValue = "" ) const;

    Bool                            SetParam( const char* name, const TFloat64List& v );
    Bool                            SetParam( const char* name, const TInt64List& v );
    Bool                            SetParam( const char* name, const TBoolList& v );
    Bool                            SetParam( const char* name, Float64 v );
    Bool                            SetParam( const char* name, Int64 v );
    Bool                            SetParam( const char* name, Bool v );
    Bool                            SetParam( const char* name, const std::string& v );

    Void                            StoreScopedPerTemplateResult( const CTemplate& t, const char* name, const COperatorResult& result );
    Void                            StoreScopedGlobalResult( const char* name, const COperatorResult& result );

    const COperatorResult*          GetPerTemplateResult( const CTemplate& t, const char* name ) const;
    TOperatorResultMap              GetPerTemplateResult( const char* name ) const;
    const COperatorResult*          GetGlobalResult( const char* name ) const;

    Void                            SaveRedshiftResultHeader( const char* dir );
    Void                            SaveRedshiftResult( const char* dir );
    Void                            SaveAllResults( const char* dir ) const;


protected:

    COperatorResultStore            m_ResultStore;

    CParameterStore                 m_ParameterStore;

    TScopeStack                     m_ScopeStack;


};


}

#endif
