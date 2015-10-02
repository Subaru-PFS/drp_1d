#ifndef _REDSHIFT_METHOD_PARAMETERSTORE_
#define _REDSHIFT_METHOD_PARAMETERSTORE_

#include <epic/core/common/datatypes.h>
#include <boost/property_tree/ptree.hpp>

namespace NSEpic
{


class CMethodParameterStore
{

public:


    CMethodParameterStore();
    virtual ~CMethodParameterStore();

    Bool Get( const char* path, const char* name, TFloat64List& v, const TFloat64List& defaultValue = TFloat64List() );
    Bool Get( const char* path, const char* name, TInt64List& v, const TInt64List& defaultValue = TInt64List() );
    Bool Get( const char* path, const char* name, TBoolList& v, const TBoolList& defaultValue = TBoolList() );
    Bool Get( const char* path, const char* name, Float64& v, Float64 defaultValue  = 0 );
    Bool Get( const char* path, const char* name, Int64& v, Int64 defaultValue = 0 );
    Bool Get( const char* path, const char* name, Bool& v, Bool defaultValue = true );

    Bool Save( const char* path ) const;
    Bool Load( const char* path );

private:

    std::string GetQualifiedName( const char* path, const char* name );

    boost::property_tree::ptree m_PropertyTree;

};


}

#endif // PARAMETERSTORE_H

