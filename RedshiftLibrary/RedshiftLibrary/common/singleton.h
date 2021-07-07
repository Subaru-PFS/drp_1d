#ifndef _REDSHIFT_COMMON_SINGLETON_
#define _REDSHIFT_COMMON_SINGLETON_

#include "RedshiftLibrary/common/datatypes.h"

namespace NSEpic {

/**
 * \ingroup Redshift
 * Singleton base class
 */
template <typename T>
class CSingleton
{
public:

    // non copyable
    CSingleton(const CSingleton &) = delete;
    CSingleton & operator= (const CSingleton &) = delete;

    static T &GetInstance() { 
      static T m_Instance{};
      return m_Instance; 
      }

protected:
    CSingleton() = default;
    ~CSingleton() = default;

};

} // namespace NSEpic

#endif
