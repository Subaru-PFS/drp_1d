#include <epic/core/common/typeinforegistry.h>


#include <cstring>
#include <typeinfo>
#include <iostream>

using namespace NSEpic;

CTypeInfo*  CTypeInfoRegistry::m_RootEntry = NULL;

CTypeInfoRegistry::CTypeInfoRegistry()
{
    CTypeInfo* n = m_RootEntry;

    while( n )
    {
        m_TypeInfoMap[ n->Name ] = n;
        n = n->NextTypeInfo;
    }
}

CTypeInfoRegistry::~CTypeInfoRegistry()
{
}

Void CTypeInfoRegistry::RegisterTypeInfo( CTypeInfo& objType )
{
    if( m_RootEntry == NULL )
    {
        m_RootEntry = &objType;
        return;
    }
    else
    {
        CTypeInfo* tmp = m_RootEntry;
        m_RootEntry = &objType;
        m_RootEntry->NextTypeInfo = tmp;
    }
}

/**
 * Retrieve a specific managed object type info according to it's class name.
 */
const CTypeInfo* CTypeInfoRegistry::GetTypeInfo( String typeName ) const
{
    TypeInfoMap::const_iterator it = m_TypeInfoMap.find( typeName );
    if( it == m_TypeInfoMap.end() )
        return NULL;

    CTypeInfo* icb = it->second;

    return icb;
}

