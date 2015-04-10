#ifndef _REDSHIFT_PROCESSFLOW_PROCESSFLOW_
#define _REDSHIFT_PROCESSFLOW_PROCESSFLOW_

#include <epic/core/common/ref.h>
#include <epic/core/common/managedobject.h>
#include <epic/core/common/threadpool.h>
#include <epic/redshift/processflow/processflow.h>

#include <boost/thread.hpp>

namespace NSEpic
{

class CProcessFlowContext;

class CProcessFlow : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( CProcessFlow )

public:

    CProcessFlow( Int32 nbThread = 0 );
    ~CProcessFlow();

    bool Process( CProcessFlowContext& ctx );

private:

    bool ProcessWithoutEL( CProcessFlowContext& ctx );
    bool ProcessWithEL( CProcessFlowContext& ctx );
    bool ELSearch( CProcessFlowContext& ctx );
    CThreadPool         m_ThreadPool;
    boost::mutex        m_SyncMutex;

};


}

#endif
