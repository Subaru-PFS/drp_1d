#ifndef _REDSHIFT_OPERATOR_OPERATOR_
#define _REDSHIFT_OPERATOR_OPERATOR_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>
#include <epic/redshift/processflow/result.h>

#include <vector>

namespace NSEpic
{

class CSpectrum;
class CTemplate;
class COperatorResult;

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
        nStatus_NoOverlap
    };

    typedef std::vector<EStatus> TStatusList;

    COperator();
    virtual ~COperator();

    virtual const COperatorResult* Compute( const CSpectrum& spectrum, const CTemplate& tpl,
                                            const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold ) = 0;

protected:

};


}

#endif
