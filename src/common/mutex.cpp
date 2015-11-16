#include <epic/core/common/mutex.h>

#include <epic/core/log/log.h>
#include <epic/core/debug/assert.h>

using namespace NSEpic;

CMutex::CMutex()
{

}

CMutex::~CMutex()
{

}

Void CMutex::Lock()
{
    m_Mutex.lock();
}

Void CMutex::Unlock()
{
    m_Mutex.unlock();
}
