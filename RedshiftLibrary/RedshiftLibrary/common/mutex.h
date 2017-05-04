#ifndef _CORE_COMMON_MUTEX_
#define _CORE_COMMON_MUTEX_

#include <epic/core/common/datatypes.h>

#include <boost/thread.hpp>

namespace NSEpic
{

/**
 * \ingroup Core
 * Mutex class
 */
class CMutex
{

public:

    CMutex( );
    virtual ~CMutex();

    Void Lock();
    Void Unlock();

private:

    boost::mutex        m_Mutex;
};

}

#endif
