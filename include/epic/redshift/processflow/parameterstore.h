#ifndef _REDSHIFT_METHOD_PARAMETERSTORE_
#define _REDSHIFT_METHOD_PARAMETERSTORE_

#include <epic/core/common/datatypes.h>
#include <boost/property_tree/ptree.hpp>

namespace NSEpic
{


class CParameterStore
{

public:


    CParameterStore();
    virtual ~CParameterStore();

    Bool Get( const char* path, const char* name, TFloat64List& v, const TFloat64List& defaultValue = TFloat64List() ) const;
    Bool Get( const char* path, const char* name, TInt64List& v, const TInt64List& defaultValue = TInt64List() ) const;
    Bool Get( const char* path, const char* name, TBoolList& v, const TBoolList& defaultValue = TBoolList() ) const;
    Bool Get( const char* path, const char* name, Float64& v, Float64 defaultValue  = 0 ) const;
    Bool Get( const char* path, const char* name, Int64& v, Int64 defaultValue = 0 ) const;
    Bool Get( const char* path, const char* name, Bool& v, Bool defaultValue = true ) const;
    Bool Get( const char* path, const char* name, std::string& v, std::string defaultValue = "" ) const;

    Bool Set( const char* path, const char* name, const TFloat64List& v );
    Bool Set( const char* path, const char* name, const TInt64List& v );
    Bool Set( const char* path, const char* name, const TBoolList& v );
    Bool Set( const char* path, const char* name, Float64 v );
    Bool Set( const char* path, const char* name, Int64 v );
    Bool Set( const char* path, const char* name, Bool v );
    Bool Set( const char* path, const char* name, const std::string& v );

    Bool Save( const char* path ) const;
    Bool Load( const char* path );

private:

    std::string GetQualifiedName( const char* path, const char* name ) const;


    mutable boost::property_tree::ptree m_PropertyTree;

};


}

#endif // PARAMETERSTORE_H

