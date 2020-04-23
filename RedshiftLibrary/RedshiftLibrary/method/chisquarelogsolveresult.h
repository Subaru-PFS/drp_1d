#ifndef _REDSHIFT_OPERATOR_CHISQUARELOGSOLVERESULT_
#define _REDSHIFT_OPERATOR_CHISQUARELOGSOLVERESULT_

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
class CChisquareLogSolveResult : public CChisquareSolveResult
{

public:

    CChisquareLogSolveResult();
    virtual ~CChisquareLogSolveResult();

};


}

#endif
