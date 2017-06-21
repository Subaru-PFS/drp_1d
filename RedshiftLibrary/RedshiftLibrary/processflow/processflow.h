#ifndef _REDSHIFT_PROCESSFLOW_PROCESSFLOW_
#define _REDSHIFT_PROCESSFLOW_PROCESSFLOW_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/processflow/processflow.h>

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

    Bool ProcessWithoutEL(CProcessFlowContext& ctx , const std::string& CategoryFilter = "all");
    Bool Blindsolve(CProcessFlowContext& ctx , const std::string& CategoryFilter = "all");
    Bool Fullsolve(CProcessFlowContext& ctx ,  const std::string& CategoryFilter = "all");
    Bool Chisquare(CProcessFlowContext& ctx , const std::string& CategoryFilter = "all");
    Bool ChisquareLog(CProcessFlowContext& ctx , const std::string& CategoryFilter = "all");

    Bool Correlation(CProcessFlowContext& ctx , const std::string& CategoryFilter = "all");
    Bool LineMatching( CProcessFlowContext& ctx );
    Bool LineMatching2( CProcessFlowContext& ctx );
    Bool LineModelSolve( CProcessFlowContext& ctx );
    Bool LineModelTplshapeSolve(CProcessFlowContext& ctx , const std::string &CategoryFilter = "all");
    Bool DecisionalTree7( CProcessFlowContext& ctx );
    Bool DecisionalTreeA( CProcessFlowContext& ctx );
    Bool DecisionalTreeB( CProcessFlowContext& ctx );
    Bool DecisionalTreeC( CProcessFlowContext& ctx );

};


}

#endif
