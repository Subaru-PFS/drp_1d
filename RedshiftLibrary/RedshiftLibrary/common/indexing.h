#ifndef _REDSHIFT_COMMON_INDEX_
#define _REDSHIFT_COMMON_INDEX_

//#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/log/log.h>

//#include <cmath>
#include <vector>

namespace NSEpic {

/**
 * \ingroup Redshift
 * Templated INDEX manipulation class
 */
template <typename T> class CIndexing
{

  public:
    CIndexing() {}

    Int32 getIndex(std::vector<T>& list, Float64 z)
    {
        typename std::vector<T>::iterator itr = std::lower_bound(list.begin(),list.end(), z);

        if (itr == list.end() || *itr != z)
        {
            size_t size = snprintf( nullptr, 0, "Could not find index for %f", z) + 1;
            std::unique_ptr<char[]> buf( new char[ size ] );                                                        
            snprintf( buf.get(), size, "Could not find index for %f", z);                          
            std::string _msg = std::string( buf.get(), buf.get() + size - 1 ); 
            throw std::runtime_error(_msg.c_str());
        }
        return (itr - list.begin()); 
    }

  /*private:
    T m_Begin;
    T m_End;*/
};
typedef CIndexing<Int32> TInt32Index;
typedef CIndexing<Float64> TFloat64Index;



} // namespace NSEpic

#endif
