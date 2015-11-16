#ifndef _CORE_COMMON_QUICKSORT_
#define _CORE_COMMON_QUICKSORT_

#include <epic/core/common/datatypes.h>

#include <algorithm>
#include <vector>

using namespace std;

namespace NSEpic
{

/**
 * \ingroup Core
 * Quicksort algorithm
 */
template< typename T >
class CQuickSort
{

public:
    
    CQuickSort();
    ~CQuickSort();
    
    Void Sort( T* arr, Int32 n ) const;
    Void SortIndexes( const T* arr, Int32* index, Int32 n ) const;
    
private:
    
    Void Sort( T* arr, Int32 beg, Int32 end ) const;
    Void SortIndexes( T* arr, Int32* index, Int32 beg, Int32 end ) const;

};
    
#include <epic/core/common/quicksort.hpp>

}

#endif
