#ifndef _REDSHIFT_COMMON_INDEX_
#define _REDSHIFT_COMMON_INDEX_
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/log/log.h>
#include <vector>

namespace NSEpic {

/**
 * \ingroup Redshift
 * Templated INDEX manipulation class
 */
template <typename T> class CIndexing
{

  public:

    static Int32 getIndex( std::vector<T>& list, T z)
    {
        typename std::vector<T>::iterator itr = std::find(list.begin(),list.end(), z); 
        if (itr == list.end())
        {
            size_t size = snprintf( nullptr, 0, "Could not find index for %f", z) + 1;
            std::unique_ptr<char[]> buf( new char[ size ] );                                                        
            snprintf( buf.get(), size, "Could not find index for %f", z);                          
            std::string _msg = std::string( buf.get(), buf.get() + size - 1 ); 
            throw std::runtime_error(_msg.c_str());
        }
        return (itr - list.begin()); 
    }

    //getIndex in orded_values corresponding to value:
    //value[index] can be equal or smaller than Z
    static bool getClosestLowerIndex(std::vector<T>& ordered_values, const T& value, Int32& i_min) 
    {
      if(value < ordered_values.front() || value > ordered_values.back())
        {
          return false;
        }
      
      typename std::vector<T>::iterator it_min = std::lower_bound(ordered_values.begin(), ordered_values.end(), value);
      if(*it_min > value) it_min = it_min -1; 

      i_min = it_min - ordered_values.begin();
      return true;
    }
};
typedef CIndexing<Int32> TInt32Index;
typedef CIndexing<Float64> TFloat64Index;



} // namespace NSEpic

#endif
