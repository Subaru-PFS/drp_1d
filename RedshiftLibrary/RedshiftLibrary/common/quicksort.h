#ifndef _CORE_COMMON_QUICKSORT_
#define _CORE_COMMON_QUICKSORT_

#include <RedshiftLibrary/common/datatypes.h>

#include <algorithm>
#include <vector>

using namespace std;

namespace NSEpic {

/**
 * \ingroup Core
 * Quicksort algorithm
 */
template <typename T> class CQuickSort
{

  public:
    CQuickSort();
    ~CQuickSort();

    void Sort(T *arr, Int32 n) const;
    void SortIndexes(const T *arr, Int32 *index, Int32 n) const;

  private:
    void Sort(T *arr, Int32 beg, Int32 end) const;
    void SortIndexes(T *arr, Int32 *index, Int32 beg, Int32 end) const;
};

#include <RedshiftLibrary/common/quicksort.hpp>

} // namespace NSEpic

#endif
