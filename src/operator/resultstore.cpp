#include <epic/redshift/operator/resultstore.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/template/template.h>

using namespace NSEpic;

COperatorResultStore::COperatorResultStore()
{

}

COperatorResultStore::~COperatorResultStore()
{

}

void COperatorResultStore::StoreResult( TResultsMap& map, const char* name, const COperatorResult& result )
{
    TResultsMap::iterator it = map.find( name );
    if( it != map.end() )
    {
        DebugError( "Resulat already exist");
        return;
    }

    map[ name ] = &result;
}

Void COperatorResultStore::StorePerTemplateResult( const CTemplate& t, const char* name, const COperatorResult& result )
{
    TPerTemplateResultsMap::iterator it = m_PerTemplateResults.find( t.GetName() );
    if( it == m_PerTemplateResults.end() )
    {
        m_PerTemplateResults[ t.GetName() ] = TResultsMap();
    }

    StoreResult( m_PerTemplateResults[ t.GetName() ], name, result );
}

Void COperatorResultStore::StoreGlobalResult( const char* name, const COperatorResult& result )
{
    StoreResult( m_GlobalResults, name, result );
}

const COperatorResult* COperatorResultStore::GetPerTemplateResult( const CTemplate& t, const char* name ) const
{
    TPerTemplateResultsMap::const_iterator it1 = m_PerTemplateResults.find( t.GetName() );
    if( it1 != m_PerTemplateResults.end() )
    {
        const TResultsMap& m = (*it1).second;
        TResultsMap::const_iterator it2 = m.find( name );
        if( it2 != m.end() )
        {
            return (*it2).second;
        }
    }

    return NULL;
}

TOperatorResultMap COperatorResultStore::GetPerTemplateResult( const char* name ) const
{
    TOperatorResultMap map;
    TPerTemplateResultsMap::const_iterator it;
    for( it = m_PerTemplateResults.begin(); it != m_PerTemplateResults.end(); ++it )
    {
        std::string tplName = (*it).first;

        const TResultsMap& m = (*it).second;
        TResultsMap::const_iterator it2 = m.find( name );
        if( it2 != m.end() )
        {
            map[tplName] = (*it2).second;
        }
    }

    return map;
}

const COperatorResult* COperatorResultStore::GetGlobalResult( const char* name ) const
{
    TResultsMap::const_iterator it = m_GlobalResults.find( name );
    if( it != m_GlobalResults.end() )
    {
        return (*it).second;
    }

    return NULL;
}
