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

Void  CDataStore::SaveRedshiftResult( const char* dir )
{
    m_ResultStore.SaveRedshiftResult( *this, dir );
}

Void  CDataStore::SaveAllResults( const char* dir ) const
{
    m_ResultStore.SaveAllResults( *this, dir );
}

Void  CDataStore::StoreScopedPerTemplateResult( const CTemplate& t, const char* name, const COperatorResult& result )
{
    m_ResultStore.StorePerTemplateResult( t, GetCurrentScopeName().c_str(), name, result );
}

Void CDataStore::StoreScopedGlobalResult( const char* name, const COperatorResult& result )
{
    m_ResultStore.StoreGlobalResult( GetCurrentScopeName().c_str(), name, result );
}


const COperatorResult*  CDataStore::GetPerTemplateResult( const CTemplate& t, const char* name ) const
{
    return m_ResultStore.GetPerTemplateResult( t, name );
}

TOperatorResultMap CDataStore::GetPerTemplateResult( const char* name ) const
{
    return m_ResultStore.GetPerTemplateResult( name );
}

const COperatorResult* CDataStore::GetGlobalResult( const char* name ) const
{
    return m_ResultStore.GetGlobalResult( name );
}

Void CDataStore::PushScope( const char* name )
{
    m_ScopeStack.push_back( name );
}

Void CDataStore::PopScope()
{
    m_ScopeStack.pop_back();
}

std::string CDataStore::GetScope(CConstRef<COperatorResult>  result) const
{
    return m_ResultStore.GetScope( result );
}



Bool CDataStore::GetScopedParam( const char* name, TFloat64List& v, const TFloat64List& defaultValue ) const
{
    return m_ParameterStore.Get( GetCurrentScopeName().c_str(), name, v, defaultValue );
}

Bool CDataStore::GetScopedParam( const char* name, TInt64List& v, const TInt64List& defaultValue ) const
{
    return m_ParameterStore.Get( GetCurrentScopeName().c_str(), name, v, defaultValue );
}

Bool CDataStore::GetScopedParam( const char* name, TBoolList& v, const TBoolList& defaultValue ) const
{
    return m_ParameterStore.Get( GetCurrentScopeName().c_str(), name, v, defaultValue );
}

Bool CDataStore::GetScopedParam( const char* name, Float64& v, Float64 defaultValue ) const
{
    return m_ParameterStore.Get( GetCurrentScopeName().c_str(), name, v, defaultValue );
}

Bool CDataStore::GetScopedParam( const char* name, Int64& v, Int64 defaultValue ) const
{
    return m_ParameterStore.Get( GetCurrentScopeName().c_str(), name, v, defaultValue );
}

Bool CDataStore::GetScopedParam( const char* name, Bool& v, Bool defaultValue ) const
{
    return m_ParameterStore.Get( GetCurrentScopeName().c_str(), name, v, defaultValue );
}

Bool CDataStore::GetScopedParam( const char* name, std::string& v, std::string defaultValue ) const
{
    return m_ParameterStore.Get( GetCurrentScopeName().c_str(), name, v, defaultValue );
}

Bool CDataStore::SetScopedParam( const char* name, const TFloat64List& v )
{
    return m_ParameterStore.Set( GetCurrentScopeName().c_str(), name, v );
}

Bool CDataStore::SetScopedParam( const char* name, const TInt64List& v )
{
    return m_ParameterStore.Set( GetCurrentScopeName().c_str(), name, v );
}

Bool CDataStore::SetScopedParam( const char* name, const TBoolList& v )
{
    return m_ParameterStore.Set( GetCurrentScopeName().c_str(), name, v );
}

Bool CDataStore::SetScopedParam( const char* name, Float64 v )
{
    return m_ParameterStore.Set( GetCurrentScopeName().c_str(), name, v );
}

Bool CDataStore::SetScopedParam( const char* name, Int64 v )
{
    return m_ParameterStore.Set( GetCurrentScopeName().c_str(), name, v );
}

Bool CDataStore::SetScopedParam( const char* name, Bool v )
{
    return m_ParameterStore.Set( GetCurrentScopeName().c_str(), name, v );
}

Bool CDataStore::SetScopedParam( const char* name, const std::string& v )
{
    return m_ParameterStore.Set( GetCurrentScopeName().c_str(), name, v );
}





Bool CDataStore::GetParam( const char* name, TFloat64List& v, const TFloat64List& defaultValue ) const
{
    return m_ParameterStore.Get( NULL, name, v, defaultValue );
}

Bool CDataStore::GetParam( const char* name, TInt64List& v, const TInt64List& defaultValue ) const
{
    return m_ParameterStore.Get( NULL, name, v, defaultValue );
}

Bool CDataStore::GetParam( const char* name, TBoolList& v, const TBoolList& defaultValue ) const
{
    return m_ParameterStore.Get( NULL, name, v, defaultValue );
}

Bool CDataStore::GetParam( const char* name, Float64& v, Float64 defaultValue ) const
{
    return m_ParameterStore.Get( NULL, name, v, defaultValue );
}

Bool CDataStore::GetParam( const char* name, Int64& v, Int64 defaultValue ) const
{
    return m_ParameterStore.Get( NULL, name, v, defaultValue );
}

Bool CDataStore::GetParam( const char* name, Bool& v, Bool defaultValue ) const
{
    return m_ParameterStore.Get( NULL, name, v, defaultValue );
}

Bool CDataStore::GetParam( const char* name, std::string& v, std::string defaultValue ) const
{
    return m_ParameterStore.Get( NULL, name, v, defaultValue );
}

Bool CDataStore::SetParam( const char* name, const TFloat64List& v )
{
    return m_ParameterStore.Set( NULL, name, v );
}

Bool CDataStore::SetParam( const char* name, const TInt64List& v )
{
    return m_ParameterStore.Set( NULL, name, v );
}

Bool CDataStore::SetParam( const char* name, const TBoolList& v )
{
    return m_ParameterStore.Set( NULL, name, v );
}

Bool CDataStore::SetParam( const char* name, Float64 v )
{
    return m_ParameterStore.Set( NULL, name, v );
}

Bool CDataStore::SetParam( const char* name, Int64 v )
{
    return m_ParameterStore.Set( NULL, name, v );
}

Bool CDataStore::SetParam( const char* name, Bool v )
{
    return m_ParameterStore.Set( NULL, name, v );
}

Bool CDataStore::SetParam( const char* name, const std::string& v )
{
    return m_ParameterStore.Set( NULL, name, v );
}
