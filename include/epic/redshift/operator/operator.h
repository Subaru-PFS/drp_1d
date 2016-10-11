#ifndef _REDSHIFT_OPERATOR_OPERATOR_
#define _REDSHIFT_OPERATOR_OPERATOR_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>
#include <epic/redshift/processflow/result.h>
#include <epic/redshift/common/mask.h>

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

    virtual  std::shared_ptr<COperatorResult> Compute( const CSpectrum& spectrum,
                                                       const CTemplate& tpl,
                                            const TFloat64Range& lambdaRange,
                                                       const TFloat64List& redshifts,
                                                       Float64 overlapThreshold,
                                                       std::vector<CMask> additional_spcMasks,
                                                       std::string opt_interp,
                                                       Int32 opt_extinction) = 0;

protected:

};


}

#endif
