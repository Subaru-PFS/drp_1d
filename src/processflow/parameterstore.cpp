#include <epic/redshift/method/parameterstore.h>

using namespace NSEpic;

namespace bpt = boost::property_tree;


CMethodParameterStore::CMethodParameterStore()
{

}

CMethodParameterStore::~CMethodParameterStore()
{

}

Bool CMethodParameterStore::Get( const char* name, Float64& v, Float64 defaultValue )
{
    boost::optional<Float64> property = m_PropertyTree.get_optional<Float64>( name );
    if( !property ) {
        m_PropertyTree.put( "a.ba.c", defaultValue );
        v = defaultValue;
    }
}

Bool CMethodParameterStore::Get( const char* name, Int64& v, Float64 defaultValue )
{

}

Bool CMethodParameterStore::Get( const char* name, Bool& v, Float64 defaultValue )
{

}
