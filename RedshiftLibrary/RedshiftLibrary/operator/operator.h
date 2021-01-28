#ifndef _REDSHIFT_OPERATOR_OPERATOR_
#define _REDSHIFT_OPERATOR_OPERATOR_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/statistics/priorhelper.h>

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
        nStatus_LoopError,
        nStatus_InvalidProductsError,
        nStatus_NoOverlap
    };

    typedef std::vector<EStatus> TStatusList;

    COperator();
    virtual ~COperator()=0;

    virtual  std::shared_ptr<COperatorResult> Compute( const CSpectrum& spectrum,
                                                       const CTemplate& tpl,
                                                       const TFloat64Range& lambdaRange,
                                                       const TFloat64List& redshifts,
                                                       Float64 overlapThreshold,
                                                       std::vector<CMask> additional_spcMasks,
                                                       std::string opt_interp,
                                                       Int32 opt_extinction,
                                                       Int32 opt_dustFitting,
                                                       CPriorHelper::TPriorZEList logprior=CPriorHelper::TPriorZEList(),
                                                       Bool keepigmism = false,
                                                       Float64 FitDustCoeff=-1,
                                                       Float64 FitMeiksinIdx=-1) = 0;

protected:

};

}

#endif
