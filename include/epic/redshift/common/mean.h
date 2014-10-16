#ifndef _REDSHIFT_COMMON_MEAN__
#define _REDSHIFT_COMMON_MEAN__

#include <epic/core/common/datatypes.h>

using namespace std;

namespace __NS__
{

template< typename T >
class CMean
{

public:
    
    CMean();
    ~CMean();
        
    T Find( const T * a, Int32 n );

private:
    
};
    
#include <epic/redshift/common/mean.hpp>

}

#endif
