// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#ifndef _REDSHIFT_COMMON_THREADPOOL_
#define _REDSHIFT_COMMON_THREADPOOL_

#include "RedshiftLibrary/common/datatypes.h"

#include <boost/asio.hpp>
#include <boost/thread.hpp>

//#define THREAPPOOL_DEBUG

namespace NSEpic {

/**
 * \ingroup Redshift
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

    bool IsRunning()
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
