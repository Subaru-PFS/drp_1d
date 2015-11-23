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


Bool CParameterStore::Get( const std::string& name, TBoolList& v, const TBoolList& defaultValue ) const
{
  //std::cout << "Get v1" << std::endl;
    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( name );

    // If property does not exist, add it
    if( !property ) {
        CParameterStore & self = const_cast<CParameterStore &>(*this);
        self.Set( name, defaultValue );
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

Bool CParameterStore::Get( const std::string& name, TInt64List& v, const TInt64List& defaultValue ) const
{
  //std::cout << "Get v2" << std::endl;
    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( name );

    // If property does not exist, add it
    if( !property ) {
        CParameterStore & self = const_cast<CParameterStore &>(*this);
        self.Set( name, defaultValue );
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

Bool CParameterStore::Get( const std::string& name, TFloat64List& v, const TFloat64List& defaultValue ) const
{
  //std::cout << "Get v3" << std::endl;
    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( name );

    // If property does not exist, add it
    if( !property ) {
        CParameterStore & self = const_cast<CParameterStore &>(*this);
        self.Set( name, defaultValue );
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

Bool CParameterStore::Get( const std::string& name, TStringList& v, const TStringList& defaultValue ) const
{
  //std::cout << "Get v4" << std::endl;
    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( name );

    // If property does not exist, add it
    if( !property ) {
        CParameterStore & self = const_cast<CParameterStore &>(*this);
        self.Set( name, defaultValue );
        v = defaultValue;
    } else {
        v.resize( (*property).size() );

        bpt::ptree::const_iterator it;
        Int32 i=0;
        for( it=property->begin(); it != property->end(); it++ ){
            v[i++] = it->second.get_value<std::string>();
        }
    }

    return true;
}

Bool CParameterStore::Get( const std::string& name, std::string& v, std::string defaultValue ) const
{
  //std::cout << "Get v5" << std::endl;
  boost::optional< std::string > property = m_PropertyTree.get_optional< std::string >( name );
  //std::cout << "got property tree" << std::endl;

  // If property does not exist, add it
  if( !property ) {
    //std::cout << "!property" << std::endl;
    CParameterStore & self = const_cast<CParameterStore &>(*this);
    self.Set( name, defaultValue );
    v = defaultValue;
  } else {
    //std::cout << "property" << std::endl;
    v = *property;
  }
  //std::cout << "Returning from CParameterStore::Get v5" << std::endl;
  return true;
}

Bool CParameterStore::Get( const std::string& name, Float64& v, Float64 defaultValue ) const
{
  //std::cout << "Get v6" << std::endl;
    boost::optional<Float64> property = m_PropertyTree.get_optional<Float64>( name );

    // If property does not exist, add it
    if( !property ) {
        CParameterStore & self = const_cast<CParameterStore &>(*this);
        self.Set( name, defaultValue );
        v = defaultValue;
    } else {
        v = *property;
    }

    return true;
}

Bool CParameterStore::Get( const std::string& name, Int64& v, Int64 defaultValue ) const
{
  //std::cout << "Get v7" << std::endl;
    boost::optional<Int64> property = m_PropertyTree.get_optional<Int64>( name );

    // If property does not exist, add it
    if( !property ) {
        CParameterStore & self = const_cast<CParameterStore &>(*this);
        self.Set( name, defaultValue );
        v = defaultValue;
    } else {
        v = *property;
    }

    return true;
}

Bool CParameterStore::Get( const std::string& name, Bool& v, Bool defaultValue ) const
{
  //std::cout << "Get v8" << std::endl;
    boost::optional<Bool> property = m_PropertyTree.get_optional<Bool>( name );

    // If property does not exist, add it
    if( !property ) {
        CParameterStore & self = const_cast<CParameterStore &>(*this);
        self.Set( name, defaultValue );
        v = defaultValue;
    } else {
        v = *property;
    }

    return true;
}


Bool CParameterStore::Get( const std::string& name, TFloat64Range& v, TFloat64Range defaultValue ) const
{
  //std::cout << "Get v9" << std::endl;
    TFloat64List listDefault( 2 );
    TFloat64List list( 2 );

    listDefault[0] = defaultValue.GetBegin();
    listDefault[1] = defaultValue.GetEnd();

    Bool r = Get( name, list, listDefault );

    v.Set( list[0], list[1] );

    return r;
}


Bool CParameterStore::Set( const std::string& name, const TFloat64List& v )
{
    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( name );

    bpt::ptree array;

    for( Int32 i=0; i<v.size(); i++ ){
        bpt::ptree item;
        item.put( "", v[i]);

        array.push_back( std::make_pair( "", item) );
    }

    m_PropertyTree.add_child( name, array );
    return true;
}

Bool CParameterStore::Set( const std::string& name, const TStringList& v )
{
    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( name );

    bpt::ptree array;

    for( Int32 i=0; i<v.size(); i++ ){
        bpt::ptree item;
        item.put( "", v[i]);

        array.push_back( std::make_pair( "", item) );
    }

    m_PropertyTree.add_child( name, array );
    return true;
}

Bool CParameterStore::Set( const std::string& name, const TFloat64Range& v )
{
    TFloat64List list( 2 );

    list[0] = v.GetBegin();
    list[1] = v.GetEnd();

    return Set( name, list );

}

Bool CParameterStore::Set( const std::string& name, const TInt64List& v )
{
    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( name );

    bpt::ptree array;

    for( Int32 i=0; i<v.size(); i++ ){
        bpt::ptree item;
        item.put( "", v[i]);

        array.push_back( std::make_pair( "", item) );
    }

    m_PropertyTree.add_child( name, array );
    return true;
}

Bool CParameterStore::Set( const std::string& name, const TBoolList& v )
{
    boost::optional< bpt::ptree & > property = m_PropertyTree.get_child_optional( name );

    bpt::ptree array;

    for( Int32 i=0; i<v.size(); i++ ){
        bpt::ptree item;
        item.put( "", v[i]);

        array.push_back( std::make_pair( "", item) );
    }

    m_PropertyTree.add_child( name, array );
    return true;
}

Bool CParameterStore::Set( const std::string& name, Float64 v )
{
    boost::optional<Float64> property = m_PropertyTree.get_optional<Float64>( name );

    m_PropertyTree.put( name, v );
    return true;
}

Bool CParameterStore::Set( const std::string& name, Int64 v )
{
    boost::optional<Int64> property = m_PropertyTree.get_optional<Int64>( name );

    m_PropertyTree.put( name, v );
    return true;
}

Bool CParameterStore::Set( const std::string& name, Bool v )
{
    boost::optional<Bool> property = m_PropertyTree.get_optional<Bool>( name );

    m_PropertyTree.put( name, v );
    return true;
}

Bool CParameterStore::Set( const std::string& name, const std::string& v )
{
    boost::optional<std::string> property = m_PropertyTree.get_optional<std::string>( name );

    m_PropertyTree.put( name, v );
    return true;
}


Bool CParameterStore::Save( const std::string& path ) const
{
    bpt::json_parser::write_json( path, m_PropertyTree );
    return true;
}

Bool CParameterStore::Load( const std::string& path )
{
    bpt::json_parser::read_json( path, m_PropertyTree );
    return true;
}
