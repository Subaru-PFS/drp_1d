#ifndef _REDSHIFT_OPERATOR_TEMPLATECOMBINATIONRESULT_
#define _REDSHIFT_OPERATOR_TEMPLATECOMBINATIONRESULT_

#include "RedshiftLibrary/operator/templatefittingresult.h"

#include "RedshiftLibrary/processflow/result.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/operator/operator.h"

namespace NSEpic
{

class CTplCombinationResult : public CTemplateFittingResult
{

public:

    void Init( UInt32 n, Int32 nISM, Int32 nIGM, Int32 componentSize);

    //best fit results
    std::vector<TFloat64List>   FitAmplitude;
    std::vector<TFloat64List>   FitAmplitudeError;
    std::vector<TFloat64List>   FitAmplitudeSigma;
    //TFloat64List            FitDtM;//not computed for tplCombination
    std::vector<std::vector<TFloat64List>>            FitMtM; //corresponding to the covariance matrix
    //TFloat64List            LogPrior;//not yet used 
};

}

#endif
