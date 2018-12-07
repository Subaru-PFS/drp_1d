#ifndef _CORE_COMMON_THREADPOOL_
#define _CORE_COMMON_THREADPOOL_

#include <RedshiftLibrary/common/datatypes.h>

#include <boost/asio.hpp>
#include <boost/thread.hpp>

//#define THREAPPOOL_DEBUG

namespace NSEpic {

/**
 * \ingroup Core
 * ThreadPool class
 */
class CThreadPool
{

  public:
    CThreadPool(Int32 nbWorkers)
    {
        m_IOServiceWork = NULL;
        m_NbWorkers = nbWorkers;

        if (m_NbWorkers <= 0)
            return;

        m_IOServiceWork = new boost::asio::io_service::work(m_IOService);

        Int32 workers = boost::thread::hardware_concurrency();

        for (std::size_t i = 0; i < m_NbWorkers; ++i)
        {
            boost::thread *thread = new boost::thread(
                boost::bind(&boost::asio::io_service::run, &m_IOService));
            m_Threads.push_back(thread);
        }
    }

    virtual ~CThreadPool() { WaitForAllTaskToFinish(); }

    template <typename TTask> void AddTask(TTask task)
    {
        if (m_NbWorkers <= 0)
        {
            task();
            return;
        }
        m_IOService.dispatch(task);
    }

    Bool IsRunning()
    {
        for (std::list<boost::thread *>::iterator it = m_Threads.begin(),
                                                  end = m_Threads.end();
             it != end; ++it)
        {
            if ((*it)->joinable())
            {
                if ((*it)->timed_join(
                        boost::posix_time::time_duration(0, 0, 0)) == false)
                {
                    return true;
                }
            }
        }
        return false;
    }

    void Reset()
    {
        if (m_NbWorkers)
            return;

        if (m_IOServiceWork == NULL)
            m_IOServiceWork = new boost::asio::io_service::work(m_IOService);

        m_IOService.reset();
    }
    void Stop()
    {
        if (m_IOServiceWork)
        {
            delete m_IOServiceWork;
            m_IOServiceWork = NULL;
        }
    }
    void WaitForAllTaskToFinish()
    {
        if (m_NbWorkers <= 0)
            return;

        Stop();

        for (std::list<boost::thread *>::iterator it = m_Threads.begin(),
                                                  end = m_Threads.end();
             it != end; ++it)
        {
            if ((*it)->joinable())
                (*it)->join();
        }
    }

  private:
    Int32 m_NbWorkers;
    boost::asio::io_service m_IOService;

    std::list<boost::thread *> m_Threads;

    boost::asio::io_service::work *m_IOServiceWork;
};

} // namespace NSEpic

#endif
