#include <epic/redshift/processflow/parameterstore.h>

#include <boost/property_tree/json_parser.hpp>

using namespace NSEpic;

namespace bpt = boost::property_tree;


CParameterStore::CParameterStore()
{

}

CParameterStore::~CParameterStore()
{

}



Bool CParameterStore::Get( const char* path, const char* name, TBoolList& v, const TBoolList& defaultValue ) const
{
    std::string finalPath = GetQualifiedName( path, name );

    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( finalPath );

    // If property does not exist, add it
    if( !property ) {
        CParameterStore & self = const_cast<CParameterStore &>(*this);
        self.Set( path, name, defaultValue );
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

Bool CParameterStore::Get( const char* path, const char* name, TInt64List& v, const TInt64List& defaultValue ) const
{
    std::string finalPath = GetQualifiedName( path, name );

    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( finalPath );

    // If property does not exist, add it
    if( !property ) {
        CParameterStore & self = const_cast<CParameterStore &>(*this);
        self.Set( path, name, defaultValue );
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

Bool CParameterStore::Get( const char* path, const char* name, TFloat64List& v, const TFloat64List& defaultValue ) const
{
    std::string finalPath = GetQualifiedName( path, name );

    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( finalPath );

    // If property does not exist, add it
    if( !property ) {
        CParameterStore & self = const_cast<CParameterStore &>(*this);
        self.Set( path, name, defaultValue );
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

Bool CParameterStore::Get( const char* path, const char* name, std::string& v, std::string defaultValue ) const
{
    std::string finalPath = GetQualifiedName( path, name );

    boost::optional< std::string > property = m_PropertyTree.get_optional< std::string >( finalPath );

    // If property does not exist, add it
    if( !property ) {
        CParameterStore & self = const_cast<CParameterStore &>(*this);
        self.Set( path, name, defaultValue );
        v = defaultValue;
    } else {
        v = *property;
    }

    return true;
}

Bool CParameterStore::Get( const char* path, const char* name, Float64& v, Float64 defaultValue ) const
{
    std::string finalPath = GetQualifiedName( path, name );

    boost::optional<Float64> property = m_PropertyTree.get_optional<Float64>( finalPath );

    // If property does not exist, add it
    if( !property ) {
        CParameterStore & self = const_cast<CParameterStore &>(*this);
        self.Set( path, name, defaultValue );
        v = defaultValue;
    } else {
        v = *property;
    }

    return true;
}

Bool CParameterStore::Get( const char* path, const char* name, Int64& v, Int64 defaultValue ) const
{
    std::string finalPath = GetQualifiedName( path, name );

    boost::optional<Int64> property = m_PropertyTree.get_optional<Int64>( finalPath );

    // If property does not exist, add it
    if( !property ) {
        CParameterStore & self = const_cast<CParameterStore &>(*this);
        self.Set( path, name, defaultValue );
        v = defaultValue;
    } else {
        v = *property;
    }

    return true;
}

Bool CParameterStore::Get( const char* path, const char* name, Bool& v, Bool defaultValue ) const
{
    std::string finalPath = GetQualifiedName( path, name );

    boost::optional<Bool> property = m_PropertyTree.get_optional<Bool>( finalPath );

    // If property does not exist, add it
    if( !property ) {
        CParameterStore & self = const_cast<CParameterStore &>(*this);
        self.Set( path, name, defaultValue );
        v = defaultValue;
    } else {
        v = *property;
    }

    return true;
}

Bool CParameterStore::Set( const char* path, const char* name, const TFloat64List& v )
{
    std::string finalPath = GetQualifiedName( path, name );
    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( finalPath );

    bpt::ptree array;

    for( Int32 i=0; i<v.size(); i++ ){
        bpt::ptree item;
        item.put( "", v[i]);

        array.push_back( std::make_pair( "", item) );
    }

    m_PropertyTree.add_child( finalPath, array );

}

Bool CParameterStore::Set( const char* path, const char* name, const TInt64List& v )
{
    std::string finalPath = GetQualifiedName( path, name );
    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( finalPath );

    bpt::ptree array;

    for( Int32 i=0; i<v.size(); i++ ){
        bpt::ptree item;
        item.put( "", v[i]);

        array.push_back( std::make_pair( "", item) );
    }

    m_PropertyTree.add_child( finalPath, array );

}

Bool CParameterStore::Set( const char* path, const char* name, const TBoolList& v )
{
    std::string finalPath = GetQualifiedName( path, name );
    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( finalPath );

    bpt::ptree array;

    for( Int32 i=0; i<v.size(); i++ ){
        bpt::ptree item;
        item.put( "", v[i]);

        array.push_back( std::make_pair( "", item) );
    }

    m_PropertyTree.add_child( finalPath, array );

    return true;
}

Bool CParameterStore::Set( const char* path, const char* name, Float64 v )
{
    std::string finalPath = GetQualifiedName( path, name );
    boost::optional<Float64> property = m_PropertyTree.get_optional<Float64>( finalPath );

    m_PropertyTree.put( finalPath, v );
    return true;
}

Bool CParameterStore::Set( const char* path, const char* name, Int64 v )
{
    std::string finalPath = GetQualifiedName( path, name );
    boost::optional<Int64> property = m_PropertyTree.get_optional<Int64>( finalPath );

    m_PropertyTree.put( finalPath, v );
    return true;
}

Bool CParameterStore::Set( const char* path, const char* name, Bool v )
{
    std::string finalPath = GetQualifiedName( path, name );
    boost::optional<Bool> property = m_PropertyTree.get_optional<Bool>( finalPath );

    m_PropertyTree.put( finalPath, v );
    return true;
}

Bool CParameterStore::Set( const char* path, const char* name, const std::string& v )
{
    std::string finalPath = GetQualifiedName( path, name );
    boost::optional<std::string> property = m_PropertyTree.get_optional<std::string>( finalPath );

    m_PropertyTree.put( finalPath, v );
    return true;
}


std::string CParameterStore::GetQualifiedName( const char* path, const char* name ) const
{
    std::string finalPath;
    if( path ) {
        finalPath = path;
        finalPath += ".";
    }
    finalPath += name;

    return finalPath;
}


Bool CParameterStore::Save( const char* path ) const
{
    bpt::json_parser::write_json( path, m_PropertyTree );
    return true;
}

Bool CParameterStore::Load( const char* path )
{
    bpt::json_parser::read_json( path, m_PropertyTree );
    return true;
}
