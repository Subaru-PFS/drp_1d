#ifndef _REDSHIFT_COMMON_MEAN__
#define _REDSHIFT_COMMON_MEAN__

#include <epic/core/common/datatypes.h>

using namespace std;

namespace NSEpic
{

  /**
   * \ingroup Redshift
   * Statistical mean objects.
   **/
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
