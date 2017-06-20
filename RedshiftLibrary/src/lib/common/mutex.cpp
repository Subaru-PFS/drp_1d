#include <RedshiftLibrary/common/mutex.h>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/debug/assert.h>

using namespace NSEpic;

/**
 *
 */
CMutex::CMutex()
{

}

/**
 *
 */
CMutex::~CMutex()
{

}

/**
 *
 */
Void CMutex::Lock()
{
    m_Mutex.lock();
}

/**
 *
 */
Void CMutex::Unlock()
{
    m_Mutex.unlock();
}
