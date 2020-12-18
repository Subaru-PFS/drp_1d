#ifndef _REDSHIFT_COMMON_MUTEX_
#define _REDSHIFT_COMMON_MUTEX_

#include <RedshiftLibrary/common/datatypes.h>

#include <boost/thread.hpp>

namespace NSEpic {

/**
 * \ingroup Redshift
 * Mutex class
 */
class CMutex
{

  public:
    CMutex();
    virtual ~CMutex();

    void Lock();
    void Unlock();

  private:
    boost::mutex m_Mutex;
};

} // namespace NSEpic

#endif
