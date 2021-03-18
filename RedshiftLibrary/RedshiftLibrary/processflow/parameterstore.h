#ifndef _REDSHIFT_PROCESSFLOW_PARAMETERSTORE_
#define _REDSHIFT_PROCESSFLOW_PARAMETERSTORE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/processflow/scopestore.h>
#include <boost/property_tree/ptree.hpp>
#include <RedshiftLibrary/common/exception.h>
#include <RedshiftLibrary/common/formatter.h>
#include <boost/property_tree/json_parser.hpp>

namespace bpt = boost::property_tree;
namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class CParameterStore : public CScopeStore
{

public:

    CParameterStore(const TScopeStack& stack);

    void Get( const std::string& name, TFloat64List& v, const TFloat64List& defaultValue = TFloat64List() ) const;
    void Get( const std::string& name, TInt64List& v, const TInt64List& defaultValue = TInt64List() ) const;
    void Get( const std::string& name, TBoolList& v, const TBoolList& defaultValue = TBoolList() ) const;
    void Get( const std::string& name, TStringList& v, const TStringList& defaultValue = TStringList() ) const;
    void Get( const std::string& name, Float64& v, Float64 defaultValue  = 0 ) const;
    void Get( const std::string& name, Int64& v, Int64 defaultValue = 0 ) const;
    void Get( const std::string& name, Bool& v, Bool defaultValue = true ) const;
    void Get( const std::string& name, TFloat64Range& v, const TFloat64Range& defaultValue = TFloat64Range( 0.0, 0.0 ) ) const;
    void Get( const std::string& name, std::string& v, const std::string& defaultValue = "" ) const;

    void GetScopedParam( const std::string& name, TFloat64Range& v, const TFloat64Range& defaultValue = TFloat64Range() ) const;
    void GetScopedParam( const std::string& name, TFloat64List& v, const TFloat64List& defaultValue = TFloat64List() ) const;
    void GetScopedParam( const std::string& name, TInt64List& v, const TInt64List& defaultValue = TInt64List() ) const;
    void GetScopedParam( const std::string& name, TBoolList& v, const TBoolList& defaultValue = TBoolList() ) const;
    void GetScopedParam( const std::string& name, Float64& v, Float64 defaultValue  = 0 ) const;
    void GetScopedParam( const std::string& name, Int64& v, Int64 defaultValue = 0 ) const;
    void GetScopedParam( const std::string& name, Bool& v, Bool defaultValue = true ) const;
    void GetScopedParam( const std::string& name, std::string& v, std::string defaultValue = "" ) const;

  template<typename T> bool Has(const std::string& name) const
  {
    boost::optional<T> property = m_PropertyTree.get_optional<T>( name );
    if(!property.is_initialized())
      {
        boost::optional<bpt::ptree &> propertree = m_PropertyTree.get_child_optional( name );
        if(!propertree.is_initialized()) return false;
      }
    return true;
  }
  
  template<typename T> T Get(const std::string& name) const
  {
    boost::optional<T> property = m_PropertyTree.get_optional<T>( name );
    if(!property.is_initialized()) throw GlobalException(UNKNOWN_PARAMETER,Formatter()<<"unknown parameter "<<name);
    return *property;
  }
  
  template<typename T> T GetScoped(const std::string& name) const
  {
    if (Has<T>(GetScopedName(name)))    return Get<T>(GetScopedName(name));    
    else return Get<T>(name);
  }

  
  template<typename T> std::vector<T> GetList(const std::string& name) const
  {
    boost::optional<bpt::ptree &> property = m_PropertyTree.get_child_optional( name );
    if(!property.is_initialized()) throw GlobalException(UNKNOWN_PARAMETER,Formatter()<<"unknown parameter "<<name);
    std::vector<T> ret = std::vector<T>( (*property).size());
  
    bpt::ptree::const_iterator it;
    Int32 i=0;
    for( it=property->begin(); it != property->end(); it++ ){
      ret[i++] = it->second.get_value<T>();
    }
    return ret;  
  }
  
  template<typename T> T GetListScoped(const std::string& name) const
  {
    if (Has<T>(GetScopedName(name)))    return GetList<T>(GetScopedName(name));
    else return GetList<T>(name);
  }


    
    void Set( const std::string& name, const TFloat64List& v );
    void Set( const std::string& name, const TInt64List& v );
    void Set( const std::string& name, const TBoolList& v );
    void Set( const std::string& name, const TStringList& v );
    void Set( const std::string& name, const TFloat64Range& v );
    void Set( const std::string& name, const std::string& v );
    void Set( const std::string& name, Float64 v );
    void Set( const std::string& name, Int64 v );
    void Set( const std::string& name, Bool v );

    void Save( const std::string& path ) const;
    void Load( const std::string& path );
    void FromString(const std::string& json);
private:

    mutable boost::property_tree::ptree m_PropertyTree;
  
};


}

#endif // PARAMETERSTORE_H
