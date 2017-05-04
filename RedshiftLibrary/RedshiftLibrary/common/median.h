#ifndef _REDSHIFT_COMMON_MEDIAN__
#define _REDSHIFT_COMMON_MEDIAN__

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/quicksort.h>

#include <algorithm>
#include <vector>

#define MEDIAN_FAST_OR_BEERS_THRESHOLD (1000)

using namespace std;

namespace NSEpic
{

  /**
   * \ingroup Redshift
   * Statistical median objects.
   */
template< typename T >
class CMedian
{

public:
    
    CMedian();
    ~CMedian();
        
    T Find( const T * a, Int32 n );

private:
    
    T FastFind( const T* a, Int32 n );
    T BeersFind( const T* a, Int32 n );
    T Opt3Find( const T  *a );
    T Opt5Find( const T  *a );
    T Opt7Find( const T  *a );
    T Opt9Find( const T  *a );


};
    
#include <RedshiftLibrary/common/median.hpp>

}

#endif
