#ifndef _REDSHIFT_COMMON_INDEX_
#define _REDSHIFT_COMMON_INDEX_
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/log/log.h"
#include <vector>
#include<iostream>
namespace NSEpic {

/**
 * \ingroup Redshift
 * Templated INDEX manipulation class
 */
template <typename T> class CIndexing
{

  public:

    static Int32 getIndex(const std::vector<T>& list, const T z)
    {
        typename std::vector<T>::const_iterator itr = std::find(list.begin(),list.end(), z); 
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
    static bool getClosestLowerIndex(const std::vector<T>& ordered_values, const T& value, Int32& i_min) 
    {
      if(value < ordered_values.front())
        {
          return false;
        }
      typename std::vector<T>::const_iterator it_min = std::lower_bound(ordered_values.begin(), ordered_values.end(), value);
      if (it_min == ordered_values.end()) it_min--;
      if(*it_min > value) it_min = it_min -1; 

      i_min = it_min - ordered_values.begin();
      return true;
    }

    //the closet at left or right, at epsilon
    static Int32 getCloserIndex(const std::vector<T>& ordered_values, const T& value) 
    {
      typename std::vector<T>::const_iterator it = std::lower_bound(ordered_values.begin(),ordered_values.end(),value);

      //check if referring to the last element
      if(it == ordered_values.end())
          it = it-1;

      else if( it != ordered_values.begin()){
        //compare diff between value and it and it-1 --> select the it that gives the minimal difference
        if(std::abs(*it - value) > std::abs( *(it-1) - value)) 
          it = it -1;
      }

      Int32 i_min = it - ordered_values.begin();
      return i_min;
    }

};
typedef CIndexing<Int32> TInt32Index;
typedef CIndexing<Float64> TFloat64Index;



} // namespace NSEpic

#endif
