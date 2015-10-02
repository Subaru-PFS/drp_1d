#ifndef _REDSHIFT_PROCESSFLOW_DATASTORE_
#define _REDSHIFT_PROCESSFLOW_DATASTORE_

#include <epic/core/common/datatypes.h>

#include <epic/redshift/processflow/parameterstore.h>
#include <epic/redshift/processflow/resultstore.h>

#include <boost/filesystem.hpp>
#include <vector>
#include <ostream>

namespace NSEpic
{

class CDataStore
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

    Void PushScope( const char* name );
    Void PopScope();

    std::string GetCurrentScopeName() const;

    // Wrapper functions
    Bool                            GetScopedParam( const char* name, TFloat64List& v, const TFloat64List& defaultValue = TFloat64List() );
    Bool                            GetScopedParam( const char* name, TInt64List& v, const TInt64List& defaultValue = TInt64List() );
    Bool                            GetScopedParam( const char* name, TBoolList& v, const TBoolList& defaultValue = TBoolList() );
    Bool                            GetScopedParam( const char* name, Float64& v, Float64 defaultValue  = 0 );
    Bool                            GetScopedParam( const char* name, Int64& v, Int64 defaultValue = 0 );
    Bool                            GetScopedParam( const char* name, Bool& v, Bool defaultValue = true );

protected:

    //COperatorResultStore            m_ResultStore;
    CParameterStore                 m_ParameterStore;

    TScopeStack                     m_ScopeStack;

};


}

#endif
