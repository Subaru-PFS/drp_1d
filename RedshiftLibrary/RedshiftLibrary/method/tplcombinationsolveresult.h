#ifndef _REDSHIFT_OPERATOR_TPLCOMBINATIONSOLVERESULT_
#define _REDSHIFT_OPERATOR_TPLCOMBINATIONSOLVERESULT_

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
class CTplcombinationSolveResult : public CChisquareSolveResult
{

public:

    CTplcombinationSolveResult();
    virtual ~CTplcombinationSolveResult();

};


}

#endif
