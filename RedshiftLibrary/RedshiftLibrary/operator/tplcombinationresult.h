#ifndef _REDSHIFT_OPERATOR_TPLCOMBINATIONRESULT_
#define _REDSHIFT_OPERATOR_TPLCOMBINATIONRESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/chisquareresult.h>

namespace NSEpic
{

class CTplcombinationResult : public CChisquareResult
{

public:

    CTplcombinationResult();
    virtual ~CTplcombinationResult();

};


}

#endif
