#ifndef _REDSHIFT_PROCESSFLOW_PROCESSFLOW_
#define _REDSHIFT_PROCESSFLOW_PROCESSFLOW_

#include <epic/core/common/ref.h>
#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>
#include <epic/core/common/managedobject.h>
#include <epic/redshift/processflow/processflow.h>

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

    Bool ProcessWithoutEL( CProcessFlowContext& ctx );
    Bool ProcessLineMatching( CProcessFlowContext& ctx );
    Bool ProcessDecisionalTree7( CProcessFlowContext& ctx );

    bool ComputeMerits( CProcessFlowContext& ctx, const TFloat64List& redshifts);

    Bool BlindSolve( CProcessFlowContext& ctx, const CTemplate& tpl, const CTemplate& tplWithoutCont  );
    Bool FullSolve( CProcessFlowContext& ctx, const CTemplate& tpl, const CTemplate& tplWithoutCont  );
    Bool FullSolveBrute( CProcessFlowContext& ctx, const CTemplate& tpl, const CTemplate& tplWithoutCont  );
};


}

#endif
