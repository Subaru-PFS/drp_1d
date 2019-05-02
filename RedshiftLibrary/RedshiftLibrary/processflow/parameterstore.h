#ifndef _REDSHIFT_METHOD_PARAMETERSTORE_
#define _REDSHIFT_METHOD_PARAMETERSTORE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <boost/property_tree/ptree.hpp>

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class CParameterStore
{

public:

    CParameterStore();
    virtual ~CParameterStore();

    void Get( const std::string& name, TFloat64List& v, const TFloat64List& defaultValue = TFloat64List() ) const;
    void Get( const std::string& name, TInt64List& v, const TInt64List& defaultValue = TInt64List() ) const;
    void Get( const std::string& name, TBoolList& v, const TBoolList& defaultValue = TBoolList() ) const;
    void Get( const std::string& name, TStringList& v, const TStringList& defaultValue = TStringList() ) const;
    void Get( const std::string& name, Float64& v, Float64 defaultValue  = 0 ) const;
    void Get( const std::string& name, Int64& v, Int64 defaultValue = 0 ) const;
    void Get( const std::string& name, Bool& v, Bool defaultValue = true ) const;
    void Get( const std::string& name, TFloat64Range& v, const TFloat64Range& defaultValue = TFloat64Range( 0.0, 0.0 ) ) const;
    void Get( const std::string& name, std::string& v, const std::string& defaultValue = "" ) const;

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
