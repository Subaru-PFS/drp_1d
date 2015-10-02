#include <epic/redshift/processflow/parameterstore.h>

#include <boost/property_tree/json_parser.hpp>

using namespace NSEpic;

namespace bpt = boost::property_tree;


CMethodParameterStore::CMethodParameterStore()
{

}

CMethodParameterStore::~CMethodParameterStore()
{

}

Bool CMethodParameterStore::Get( const char* path, const char* name, TBoolList& v, const TBoolList& defaultValue )
{
    std::string finalPath = GetQualifiedName( path, name );

    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( finalPath );

    // If property does not exist, add it
    if( !property ) {
        bpt::ptree array;

        for( Int32 i=0; i<defaultValue.size(); i++ ){
            bpt::ptree item;
            item.put( "", defaultValue[i]);

            array.push_back( std::make_pair( "", item) );
        }

        m_PropertyTree.add_child( finalPath, array );
        v = defaultValue;
    } else {
        v.resize( (*property).size() );

        bpt::ptree::const_iterator it;
        Int32 i=0;
        for( it=property->begin(); it != property->end(); it++ ){
            v[i++] = it->second.get_value<Bool>();
        }
    }

    return true;
}

Bool CMethodParameterStore::Get( const char* path, const char* name, TInt64List& v, const TInt64List& defaultValue )
{
    std::string finalPath = GetQualifiedName( path, name );

    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( finalPath );

    // If property does not exist, add it
    if( !property ) {
        bpt::ptree array;

        for( Int32 i=0; i<defaultValue.size(); i++ ){
            bpt::ptree item;
            item.put( "", defaultValue[i]);

            array.push_back( std::make_pair( "", item) );
        }

        m_PropertyTree.add_child( finalPath, array );
        v = defaultValue;
    } else {
        v.resize( (*property).size() );

        bpt::ptree::const_iterator it;
        Int32 i=0;
        for( it=property->begin(); it != property->end(); it++ ){
            v[i++] = it->second.get_value<Int64>();
        }
    }

    return true;
}

Bool CMethodParameterStore::Get( const char* path, const char* name, TFloat64List& v, const TFloat64List& defaultValue )
{
    std::string finalPath = GetQualifiedName( path, name );

    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( finalPath );

    // If property does not exist, add it
    if( !property ) {
        bpt::ptree array;

        for( Int32 i=0; i<defaultValue.size(); i++ ){
            bpt::ptree item;
            item.put( "", defaultValue[i]);

            array.push_back( std::make_pair( "", item) );
        }

        m_PropertyTree.add_child( finalPath, array );
        v = defaultValue;
    } else {
        v.resize( (*property).size() );

        bpt::ptree::const_iterator it;
        Int32 i=0;
        for( it=property->begin(); it != property->end(); it++ ){
            v[i++] = it->second.get_value<Float64>();
        }
    }

    return true;
}

Bool CMethodParameterStore::Get( const char* path, const char* name, Float64& v, Float64 defaultValue )
{
    std::string finalPath = GetQualifiedName( path, name );

    boost::optional<Float64> property = m_PropertyTree.get_optional<Float64>( finalPath );

    // If property does not exist, add it
    if( !property ) {
        m_PropertyTree.put( finalPath, defaultValue );
        v = defaultValue;
    } else {
        v = *property;
    }

    return true;
}

Bool CMethodParameterStore::Get( const char* path, const char* name, Int64& v, Int64 defaultValue )
{
    std::string finalPath = GetQualifiedName( path, name );

    boost::optional<Int64> property = m_PropertyTree.get_optional<Int64>( finalPath );

    // If property does not exist, add it
    if( !property ) {
        m_PropertyTree.put( finalPath, defaultValue );
        v = defaultValue;
    } else {
        v = *property;
    }

    return true;
}

Bool CMethodParameterStore::Get( const char* path, const char* name, Bool& v, Bool defaultValue )
{
    std::string finalPath = GetQualifiedName( path, name );

    boost::optional<Bool> property = m_PropertyTree.get_optional<Bool>( finalPath );

    // If property does not exist, add it
    if( !property ) {
        m_PropertyTree.put( finalPath, defaultValue );
        v = defaultValue;
    } else {
        v = *property;
    }

    return true;
}

std::string CMethodParameterStore::GetQualifiedName( const char* path, const char* name )
{
    std::string finalPath;
    if( path ) {
        finalPath = path;
        finalPath += ".";
    }
    finalPath += name;

    return finalPath;
}

Bool CMethodParameterStore::Save( const char* path ) const
{
    bpt::json_parser::write_json( path, m_PropertyTree );
    return true;
}

Bool CMethodParameterStore::Load( const char* path )
{
    bpt::json_parser::read_json( path, m_PropertyTree );
    return true;
}
