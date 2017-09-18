#ifndef _REDSHIFT_PROCESSFLOW_PROCESSFLOW_
#define _REDSHIFT_PROCESSFLOW_PROCESSFLOW_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/processflow/processflow.h>
#include <RedshiftLibrary/operator/operator.h>


#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>

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

    Bool Blindsolve(CProcessFlowContext& ctx , const std::string& CategoryFilter = "all");
    Bool Fullsolve(CProcessFlowContext& ctx ,  const std::string& CategoryFilter = "all");

    Bool Correlation(CProcessFlowContext& ctx , const std::string& CategoryFilter = "all");
    Bool LineMatching( CProcessFlowContext& ctx );
    Bool LineMatching2( CProcessFlowContext& ctx );
    Bool LineModelTplshapeSolve(CProcessFlowContext& ctx , const std::string &CategoryFilter = "all");
    Bool DecisionalTree7( CProcessFlowContext& ctx );
    Bool DecisionalTreeA( CProcessFlowContext& ctx );
    Bool DecisionalTreeB( CProcessFlowContext& ctx );
    Bool DecisionalTreeC( CProcessFlowContext& ctx );

    Bool isPdfValid(CProcessFlowContext &ctx) const;
};


}

#endif
