#ifndef _REDSHIFT_OPERATOR_CHISQUARE2SOLVERESULT_
#define _REDSHIFT_OPERATOR_CHISQUARE2SOLVERESULT_

#include <RedshiftLibrary/method/chisquaresolveresult.h>
#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/ray/catalog.h>

#include <vector>

namespace NSEpic
{

class CProcessFlowContext;

/**
 * \ingroup Redshift
 */
class CChisquare2SolveResult : public CChisquareSolveResult
{

public:

    CChisquare2SolveResult();
    virtual ~CChisquare2SolveResult();

};


}

#endif
