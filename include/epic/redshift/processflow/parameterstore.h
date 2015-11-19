#ifndef _REDSHIFT_METHOD_PARAMETERSTORE_
#define _REDSHIFT_METHOD_PARAMETERSTORE_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>
#include <boost/property_tree/ptree.hpp>

namespace NSEpic
{


class CParameterStore
{

public:

    CParameterStore();
    virtual ~CParameterStore();

    Bool Get( const std::string& name, TFloat64List& v, const TFloat64List& defaultValue = TFloat64List() ) const;
    Bool Get( const std::string& name, TInt64List& v, const TInt64List& defaultValue = TInt64List() ) const;
    Bool Get( const std::string& name, TBoolList& v, const TBoolList& defaultValue = TBoolList() ) const;
    Bool Get( const std::string& name, TStringList& v, const TStringList& defaultValue = TStringList() ) const;
    Bool Get( const std::string& name, Float64& v, Float64 defaultValue  = 0 ) const;
    Bool Get( const std::string& name, Int64& v, Int64 defaultValue = 0 ) const;
    Bool Get( const std::string& name, Bool& v, Bool defaultValue = true ) const;
    Bool Get( const std::string& name, TFloat64Range& v, TFloat64Range defaultValue = TFloat64Range( 0.0, 0.0 ) ) const;
    Bool Get( const std::string& name, std::string& v, std::string defaultValue = "" ) const;

    Bool Set( const std::string& name, const TFloat64List& v );
    Bool Set( const std::string& name, const TInt64List& v );
    Bool Set( const std::string& name, const TBoolList& v );
    Bool Set( const std::string& name, const TStringList& v );
    Bool Set( const std::string& name, const TFloat64Range& v );
    Bool Set( const std::string& name, const std::string& v );
    Bool Set( const std::string& name, Float64 v );
    Bool Set( const std::string& name, Int64 v );
    Bool Set( const std::string& name, Bool v );

    Bool Save( const std::string& path ) const;
    Bool Load( const std::string& path );

private:

    mutable boost::property_tree::ptree m_PropertyTree;

};


}

#endif // PARAMETERSTORE_H

