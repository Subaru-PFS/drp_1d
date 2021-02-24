#ifndef _REDSHIFT_OPERATOR_TEMPLATE_FITTING_BASE_
#define _REDSHIFT_OPERATOR_TEMPLATE_FITTING_BASE_

#include <RedshiftLibrary/operator/operator.h>

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
class COperatorTemplateFittingBase : public COperator
{

public:

    COperatorTemplateFittingBase();
    virtual ~COperatorTemplateFittingBase()=0;

    virtual  std::shared_ptr<COperatorResult> Compute( const CSpectrum& spectrum,
                                                       const CTemplate& tpl,
                                                      const TFloat64Range& lambdaRange,
                                                       const TFloat64List& redshifts,
                                                       Float64 overlapThreshold,
                                                       std::vector<CMask> additional_spcMasks,
                                                       std::string opt_interp,
                                                       Int32 opt_extinction,
                                                       Int32 opt_dustFitting,
                                                       CPriorHelper::TPriorZEList logprior,
                                                       Bool keepigmism = false,
                                                       Float64 FitDustCoeff=-1,
                                                       Float64 FitMeiksinIdx=-1) = 0;


};


}

#endif
