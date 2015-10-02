#include <epic/redshift/processflow/datastore.h>

#include <epic/core/debug/assert.h>

using namespace NSEpic;

CDataStore::CAutoScope::CAutoScope( CDataStore& store, const char* name )
{
    m_Store = & store;
    m_Store->PushScope( name );
}

CDataStore::CAutoScope::~CAutoScope()
{
    m_Store->PopScope();
}

CDataStore::CDataStore()
{
}

CDataStore::~CDataStore()
{

}

std::string CDataStore::GetCurrentScopeName() const
{
    std::string n;

    TScopeStack::const_iterator it;

    if( m_ScopeStack.size() == 0 )
        return n;

    n = m_ScopeStack[0];
    it = m_ScopeStack.begin();
    it++;

    for( ; it != m_ScopeStack.end(); it++ ) {
        n.append(".");
        n.append((*it));
    }

    return n;
}

Void CDataStore::PushScope( const char* name )
{
    m_ScopeStack.push_back( name );
}

Void CDataStore::PopScope()
{
    m_ScopeStack.pop_back();
}

Bool CDataStore::GetScopedParam( const char* name, TFloat64List& v, const TFloat64List& defaultValue )
{
    return m_ParameterStore.Get( GetCurrentScopeName().c_str(), name, v, defaultValue );
}

Bool CDataStore::GetScopedParam( const char* name, TInt64List& v, const TInt64List& defaultValue )
{
    return m_ParameterStore.Get( GetCurrentScopeName().c_str(), name, v, defaultValue );
}

Bool CDataStore::GetScopedParam( const char* name, TBoolList& v, const TBoolList& defaultValue )
{
    return m_ParameterStore.Get( GetCurrentScopeName().c_str(), name, v, defaultValue );
}

Bool CDataStore::GetScopedParam( const char* name, Float64& v, Float64 defaultValue )
{
    return m_ParameterStore.Get( GetCurrentScopeName().c_str(), name, v, defaultValue );
}

Bool CDataStore::GetScopedParam( const char* name, Int64& v, Int64 defaultValue )
{
    return m_ParameterStore.Get( GetCurrentScopeName().c_str(), name, v, defaultValue );
}

Bool CDataStore::GetScopedParam( const char* name, Bool& v, Bool defaultValue )
{
    return m_ParameterStore.Get( GetCurrentScopeName().c_str(), name, v, defaultValue );
}

