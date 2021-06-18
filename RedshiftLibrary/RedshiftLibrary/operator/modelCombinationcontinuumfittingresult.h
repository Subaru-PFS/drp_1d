#ifndef _REDSHIFT_LINEMODEL_TPLCOMB_MODELCONTINUUMFITTINGRESULT_
#define _REDSHIFT_LINEMODEL_TPLCOMB_MODELCONTINUUMFITTINGRESULT_


#include "RedshiftLibrary/processflow/result.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/operator/modelcontinuumfittingresult.h"

namespace NSEpic
{
class CModelContinuumFittingResult;
  /**
   * \ingroup Redshift
   */
class CModelContinuumTplCombinationFittingResult : public CModelContinuumFittingResult
{

public:

    CModelContinuumTplCombinationFittingResult(Float64 _redshift,
                                 Float64 _merit,
                                 TFloat64List _amp,
                                 TFloat64List _amp_err,
                                 Float64 _ismCoeff,
                                 Int32 _igmIndex,
                                 Float64 _fitting_snr=NAN);

    virtual ~CModelContinuumTplCombinationFittingResult(){};

private:

    TFloat64List Amp;
    TFloat64List AmpErr;
};


}

#endif
