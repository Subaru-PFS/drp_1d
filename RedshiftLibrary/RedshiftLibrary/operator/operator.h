#ifndef _REDSHIFT_OPERATOR_OPERATOR_
#define _REDSHIFT_OPERATOR_OPERATOR_

//#include <RedshiftLibrary/common/datatypes.h>
//#include <RedshiftLibrary/common/range.h>
//#include <RedshiftLibrary/processflow/result.h>
//#include <RedshiftLibrary/common/mask.h>
//#include <RedshiftLibrary/statistics/priorhelper.h>

#include <vector>

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class COperator
{

public:

    enum EStatus
    {
        nStatus_OK = 0,
        nStatus_DataError,
        nStatus_LoopError,
        nStatus_InvalidProductsError,
        nStatus_NoOverlap
    };

    typedef std::vector<EStatus> TStatusList;

    COperator();
    virtual ~COperator()=0;


protected:

};

}

#endif
