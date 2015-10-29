#include <epic/redshift/processflow/datastore.h>

#include <epic/core/debug/assert.h>

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT_NOT_INSTANCIABLE( CDataStore );

CDataStore::CAutoScope::CAutoScope( CDataStore& store, const std::string& name )
{
    m_Store = & store;
    m_Store->PushScope( name );
}

CDataStore::CAutoScope::~CAutoScope()
{
    m_Store->PopScope();
}

CDataStore::CDataStore( COperatorResultStore& resultStore, CParameterStore& parameStore ) :
   m_ResultStore( resultStore ),
   m_ParameterStore( parameStore )
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

const std::string& CDataStore::GetSpectrumName() const
{
    return m_SpectrumName;
}

Void CDataStore::SetSpectrumName( const std::string& name )
{
    m_SpectrumName = name;
}

Void  CDataStore::SaveRedshiftResult( const boost::filesystem::path& dir )
{
    m_ResultStore.SaveRedshiftResult( *this, dir );
}

Void  CDataStore::SaveAllResults( const boost::filesystem::path& dir ) const
{
    m_ResultStore.SaveAllResults( *this, dir );
}

Void  CDataStore::StoreScopedPerTemplateResult( const CTemplate& t, const std::string& name, const COperatorResult& result )
{
    m_ResultStore.StorePerTemplateResult( t, GetCurrentScopeName(), name, result );
}

Void CDataStore::StoreScopedGlobalResult( const std::string& name, const COperatorResult& result )
{
    m_ResultStore.StoreGlobalResult( GetCurrentScopeName(), name, result );
}


const COperatorResult*  CDataStore::GetPerTemplateResult( const CTemplate& t, const std::string& name ) const
{
    return m_ResultStore.GetPerTemplateResult( t, name );
}

TOperatorResultMap CDataStore::GetPerTemplateResult( const std::string& name ) const
{
    return m_ResultStore.GetPerTemplateResult( name );
}

const COperatorResult* CDataStore::GetGlobalResult( const std::string& name ) const
{
    return m_ResultStore.GetGlobalResult( name );
}

Void CDataStore::PushScope( const std::string& name )
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



std::string  CDataStore::GetScopedName( const std::string& name ) const {

    std::string scopedName = GetCurrentScopeName();

    if( ! scopedName.empty() ) {
        scopedName.append(".");
    }

    scopedName.append( name );

    return scopedName;

}

Bool CDataStore::GetScopedParam( const std::string& name, TFloat64List& v, const TFloat64List& defaultValue ) const
{
    return m_ParameterStore.Get( GetScopedName( name ), v, defaultValue );
}

Bool CDataStore::GetScopedParam( const std::string& name, TInt64List& v, const TInt64List& defaultValue ) const
{
    return m_ParameterStore.Get( GetScopedName( name ), v, defaultValue );
}

Bool CDataStore::GetScopedParam( const std::string& name, TBoolList& v, const TBoolList& defaultValue ) const
{
    return m_ParameterStore.Get( GetScopedName( name ), v, defaultValue );
}

Bool CDataStore::GetScopedParam( const std::string& name, Float64& v, Float64 defaultValue ) const
{
    return m_ParameterStore.Get( GetScopedName( name ), v, defaultValue );
}

Bool CDataStore::GetScopedParam( const std::string& name, Int64& v, Int64 defaultValue ) const
{
    return m_ParameterStore.Get( GetScopedName( name ), v, defaultValue );
}

Bool CDataStore::GetScopedParam( const std::string& name, Bool& v, Bool defaultValue ) const
{
    return m_ParameterStore.Get( GetScopedName( name ), v, defaultValue );
}

Bool CDataStore::GetScopedParam( const std::string& name, std::string& v, std::string defaultValue ) const
{
    return m_ParameterStore.Get( GetScopedName( name ), v, defaultValue );
}

Bool CDataStore::SetScopedParam( const std::string& name, const TFloat64List& v )
{
    return m_ParameterStore.Set( GetScopedName( name ), v );
}

Bool CDataStore::SetScopedParam( const std::string& name, const TInt64List& v )
{
    return m_ParameterStore.Set( GetScopedName( name ), v );
}

Bool CDataStore::SetScopedParam( const std::string& name, const TBoolList& v )
{
    return m_ParameterStore.Set( GetScopedName( name ), v );
}

Bool CDataStore::SetScopedParam( const std::string& name, Float64 v )
{
    return m_ParameterStore.Set( GetScopedName( name ), v );
}

Bool CDataStore::SetScopedParam( const std::string& name, Int64 v )
{
    return m_ParameterStore.Set( GetScopedName( name ), v );
}

Bool CDataStore::SetScopedParam( const std::string& name, Bool v )
{
    return m_ParameterStore.Set( GetScopedName( name ), v );
}

Bool CDataStore::SetScopedParam( const std::string& name, const std::string& v )
{
    return m_ParameterStore.Set( GetScopedName( name ), v );
}





Bool CDataStore::GetParam( const std::string& name, TFloat64List& v, const TFloat64List& defaultValue ) const
{
    return m_ParameterStore.Get( name, v, defaultValue );
}

Bool CDataStore::GetParam( const std::string& name, TInt64List& v, const TInt64List& defaultValue ) const
{
    return m_ParameterStore.Get( name, v, defaultValue );
}

Bool CDataStore::GetParam( const std::string& name, TBoolList& v, const TBoolList& defaultValue ) const
{
    return m_ParameterStore.Get( name, v, defaultValue );
}

Bool CDataStore::GetParam( const std::string& name, Float64& v, Float64 defaultValue ) const
{
    return m_ParameterStore.Get( name, v, defaultValue );
}

Bool CDataStore::GetParam( const std::string& name, Int64& v, Int64 defaultValue ) const
{
    return m_ParameterStore.Get( name, v, defaultValue );
}

Bool CDataStore::GetParam( const std::string& name, Bool& v, Bool defaultValue ) const
{
    return m_ParameterStore.Get( name, v, defaultValue );
}

Bool CDataStore::GetParam( const std::string& name, std::string& v, std::string defaultValue ) const
{
    return m_ParameterStore.Get( name, v, defaultValue );
}

Bool CDataStore::SetParam( const std::string& name, const TFloat64List& v )
{
    return m_ParameterStore.Set( name, v );
}

Bool CDataStore::SetParam( const std::string& name, const TInt64List& v )
{
    return m_ParameterStore.Set( name, v );
}

Bool CDataStore::SetParam( const std::string& name, const TBoolList& v )
{
    return m_ParameterStore.Set( name, v );
}

Bool CDataStore::SetParam( const std::string& name, Float64 v )
{
    return m_ParameterStore.Set( name, v );
}

Bool CDataStore::SetParam( const std::string& name, Int64 v )
{
    return m_ParameterStore.Set( name, v );
}

Bool CDataStore::SetParam( const std::string& name, Bool v )
{
    return m_ParameterStore.Set( name, v );
}

Bool CDataStore::SetParam( const std::string& name, const std::string& v )
{
    return m_ParameterStore.Set( name, v );
}
