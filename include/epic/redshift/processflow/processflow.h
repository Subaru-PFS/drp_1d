#ifndef _REDSHIFT_PROCESSFLOW_PROCESSFLOW_
#define _REDSHIFT_PROCESSFLOW_PROCESSFLOW_

#include <epic/core/common/ref.h>
#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>
#include <epic/core/common/managedobject.h>
#include <epic/redshift/processflow/processflow.h>

#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/template/catalog.h>
namespace NSEpic
{

class CProcessFlowContext;
class CTemplate;
class CSpectrum;

class CProcessFlow : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( CProcessFlow )

public:

    CProcessFlow();
    ~CProcessFlow();

    Bool Process( CProcessFlowContext& ctx );

private:

    Bool ProcessWithoutEL(CProcessFlowContext& ctx , NSEpic::CTemplate::ECategory CategoryFilter = NSEpic::CTemplate::nCategory_None);
    Bool Blindsolve(CProcessFlowContext& ctx , NSEpic::CTemplate::ECategory CategoryFilter = NSEpic::CTemplate::nCategory_None);
    Bool Fullsolve(CProcessFlowContext& ctx , NSEpic::CTemplate::ECategory CategoryFilter = NSEpic::CTemplate::nCategory_None);
    //Bool Chisolve(CProcessFlowContext& ctx , TFloat64List& redshifts, NSEpic::CTemplate::ECategory CategoryFilter = NSEpic::CTemplate::nCategory_None);
    Bool Correlation(CProcessFlowContext& ctx , NSEpic::CTemplate::ECategory CategoryFilter = NSEpic::CTemplate::nCategory_None);
    Bool LineMatching( CProcessFlowContext& ctx );
    Bool LineMatching2( CProcessFlowContext& ctx );
    Bool DecisionalTree7( CProcessFlowContext& ctx );
    Bool DecisionalTreeA( CProcessFlowContext& ctx );

};


}

#endif
