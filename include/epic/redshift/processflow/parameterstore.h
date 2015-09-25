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

    Bool Get( const char* name, Float64& v, Float64 defaultValue  = 0 );
    Bool Get( const char* name, Int64& v, Float64 defaultValue = 0 );
    Bool Get( const char* name, Bool& v, Float64 defaultValue = true );

private:


    boost::property_tree::ptree m_PropertyTree;

};


}

#endif // PARAMETERSTORE_H

