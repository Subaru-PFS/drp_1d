#include "RedshiftLibrary/common/mutex.h"

#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/debug/assert.h"

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
void CMutex::Lock()
{
    m_Mutex.lock();
}

/**
 *
 */
void CMutex::Unlock()
{
    m_Mutex.unlock();
}
