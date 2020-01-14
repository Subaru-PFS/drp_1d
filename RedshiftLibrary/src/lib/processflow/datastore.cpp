#include <RedshiftLibrary/processflow/datastore.h>

#include <RedshiftLibrary/debug/assert.h>

using namespace NSEpic;


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

void CDataStore::SetSpectrumName( const std::string& name )
{
    m_SpectrumName = name;
}

const std::string& CDataStore::GetProcessingID() const
{
    return m_ProcessingID;
}

void CDataStore::SetProcessingID( const std::string& valStr )
{
    m_ProcessingID = valStr;
}

void  CDataStore::SaveRedshiftResult( const boost::filesystem::path& dir )
{
    m_ResultStore.SaveRedshiftResult( *this, dir );
}

void  CDataStore::SaveStellarResult( const boost::filesystem::path& dir )
{
    m_ResultStore.SaveStellarResult( *this, dir );
}

void  CDataStore::SaveQsoResult( const boost::filesystem::path& dir )
{
    m_ResultStore.SaveQsoResult( *this, dir );
}

void  CDataStore::SaveClassificationResult( const boost::filesystem::path& dir )
{
    m_ResultStore.SaveClassificationResult( *this, dir );
}

void  CDataStore::SaveCandidatesResult( const boost::filesystem::path& dir )
{
    m_ResultStore.SaveCandidatesResult( *this, dir );
}


void CDataStore::SaveReliabilityResult( const boost::filesystem::path& dir )
{
	m_ResultStore.SaveReliabilityResult( *this, dir );
}

void  CDataStore::SaveAllResults( const boost::filesystem::path& dir, const std::string opt ) const
{
    m_ResultStore.SaveAllResults( *this, dir, opt );
}

void  CDataStore::StoreScopedPerTemplateResult( const CTemplate& t, const std::string& name, std::shared_ptr<const COperatorResult> result )
{
    m_ResultStore.StorePerTemplateResult( t, GetCurrentScopeName(), name, result );
}

void CDataStore::StoreScopedGlobalResult( const std::string& name, std::shared_ptr<const COperatorResult> result )
{
    m_ResultStore.StoreGlobalResult( GetCurrentScopeName(), name, result );
}

void CDataStore::DeleteScopedGlobalResult( const std::string& name )
{
    m_ResultStore.DeleteGlobalResult(GetCurrentScopeName(), name);//""
    
}
void CDataStore::ChangeScopedGlobalResult( const std::string& oldkey, const std::string& newkey )
{
    
    auto  result = GetGlobalResult( oldkey ).lock();
    StoreScopedGlobalResult(newkey, result);
    DeleteScopedGlobalResult(oldkey);
}

void CDataStore::StoreGlobalResult( const std::string& name, std::shared_ptr<const COperatorResult> result )
{
    m_ResultStore.StoreGlobalResult( "", name, result );
}

std::weak_ptr<const COperatorResult>  CDataStore::GetPerTemplateResult( const CTemplate& t, const std::string& name ) const
{
    return m_ResultStore.GetPerTemplateResult( t, name );
}

TOperatorResultMap CDataStore::GetPerTemplateResult( const std::string& name ) const
{
    return m_ResultStore.GetPerTemplateResult( name );
}

std::weak_ptr<const COperatorResult> CDataStore::GetGlobalResult( const std::string& name ) const
{
    return m_ResultStore.GetGlobalResult( name );
}

void CDataStore::PushScope( const std::string& name )
{
    m_ScopeStack.push_back( name );
}

void CDataStore::PopScope()
{
    m_ScopeStack.pop_back();
}

std::string CDataStore::GetScope( const COperatorResult&  result ) const
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

void CDataStore::GetScopedParam( const std::string& name, TFloat64List& v, const TFloat64List& defaultValue ) const
{
    return m_ParameterStore.Get( GetScopedName( name ), v, defaultValue );
}

void CDataStore::GetScopedParam( const std::string& name, TInt64List& v, const TInt64List& defaultValue ) const
{
    return m_ParameterStore.Get( GetScopedName( name ), v, defaultValue );
}

void CDataStore::GetScopedParam( const std::string& name, TBoolList& v, const TBoolList& defaultValue ) const
{
    return m_ParameterStore.Get( GetScopedName( name ), v, defaultValue );
}

void CDataStore::GetScopedParam( const std::string& name, Float64& v, Float64 defaultValue ) const
{
    return m_ParameterStore.Get( GetScopedName( name ), v, defaultValue );
}

void CDataStore::GetScopedParam( const std::string& name, Int64& v, Int64 defaultValue ) const
{
    return m_ParameterStore.Get( GetScopedName( name ), v, defaultValue );
}

void CDataStore::GetScopedParam( const std::string& name, Bool& v, Bool defaultValue ) const
{
    return m_ParameterStore.Get( GetScopedName( name ), v, defaultValue );
}

void CDataStore::GetScopedParam( const std::string& name, std::string& v, std::string defaultValue ) const
{
    return m_ParameterStore.Get( GetScopedName( name ), v, defaultValue );
}

void CDataStore::SetScopedParam( const std::string& name, const TFloat64List& v )
{
    return m_ParameterStore.Set( GetScopedName( name ), v );
}

void CDataStore::SetScopedParam( const std::string& name, const TInt64List& v )
{
    return m_ParameterStore.Set( GetScopedName( name ), v );
}

void CDataStore::SetScopedParam( const std::string& name, const TBoolList& v )
{
    return m_ParameterStore.Set( GetScopedName( name ), v );
}

void CDataStore::SetScopedParam( const std::string& name, Float64 v )
{
    return m_ParameterStore.Set( GetScopedName( name ), v );
}

void CDataStore::SetScopedParam( const std::string& name, Int64 v )
{
    return m_ParameterStore.Set( GetScopedName( name ), v );
}

void CDataStore::SetScopedParam( const std::string& name, Bool v )
{
    return m_ParameterStore.Set( GetScopedName( name ), v );
}

void CDataStore::SetScopedParam( const std::string& name, const std::string& v )
{
    return m_ParameterStore.Set( GetScopedName( name ), v );
}





void CDataStore::GetParam( const std::string& name, TFloat64List& v, const TFloat64List& defaultValue ) const
{
    return m_ParameterStore.Get( name, v, defaultValue );
}

void CDataStore::GetParam( const std::string& name, TInt64List& v, const TInt64List& defaultValue ) const
{
    return m_ParameterStore.Get( name, v, defaultValue );
}

void CDataStore::GetParam( const std::string& name, TBoolList& v, const TBoolList& defaultValue ) const
{
    return m_ParameterStore.Get( name, v, defaultValue );
}

void CDataStore::GetParam( const std::string& name, Float64& v, Float64 defaultValue ) const
{
    return m_ParameterStore.Get( name, v, defaultValue );
}

void CDataStore::GetParam( const std::string& name, Int64& v, Int64 defaultValue ) const
{
    return m_ParameterStore.Get( name, v, defaultValue );
}

void CDataStore::GetParam( const std::string& name, Bool& v, Bool defaultValue ) const
{
    return m_ParameterStore.Get( name, v, defaultValue );
}

void CDataStore::GetParam( const std::string& name, std::string& v, std::string defaultValue ) const
{
    return m_ParameterStore.Get( name, v, defaultValue );
}

void CDataStore::SetParam( const std::string& name, const TFloat64List& v )
{
    return m_ParameterStore.Set( name, v );
}

void CDataStore::SetParam( const std::string& name, const TInt64List& v )
{
    return m_ParameterStore.Set( name, v );
}

void CDataStore::SetParam( const std::string& name, const TBoolList& v )
{
    return m_ParameterStore.Set( name, v );
}

void CDataStore::SetParam( const std::string& name, Float64 v )
{
    return m_ParameterStore.Set( name, v );
}

void CDataStore::SetParam( const std::string& name, Int64 v )
{
    return m_ParameterStore.Set( name, v );
}

void CDataStore::SetParam( const std::string& name, Bool v )
{
    return m_ParameterStore.Set( name, v );
}

void CDataStore::SetParam( const std::string& name, const std::string& v )
{
    return m_ParameterStore.Set( name, v );
}
