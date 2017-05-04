#ifndef _REDSHIFT_LINEMODEL_RULES
#define _REDSHIFT_LINEMODEL_RULES

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>
#include <epic/redshift/processflow/processflow.h>

#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/template/catalog.h>

namespace NSEpic
{

class CProcessFlowContext;
class CTemplate;
class CSpectrum;

/**
 * \ingroup Redshift
 */
class CProcessFlow
{

public:

    CProcessFlow();
    ~CProcessFlow();

    Bool Process( CProcessFlowContext& ctx );

private:

    Bool LineModelSolve( CProcessFlowContext& ctx );
};


}

#endif
