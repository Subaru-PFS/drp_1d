#ifndef _REDSHIFT_OPERATOR_TEMPLATE_FITTING_
#define _REDSHIFT_OPERATOR_TEMPLATE_FITTING_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/operator/templatefittingBase.h>
#include <RedshiftLibrary/operator/templatefittingresult.h>
#include <RedshiftLibrary/common/mask.h>

#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h>
#include <RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h>
#include <RedshiftLibrary/statistics/priorhelper.h>
#include <RedshiftLibrary/operator/modelspectrumresult.h>
namespace NSEpic
{

class COperatorTemplateFitting : public COperatorTemplateFittingBase
{

public:
    explicit COperatorTemplateFitting() = default;
    ~COperatorTemplateFitting() = default;

     std::shared_ptr<COperatorResult> Compute(const CSpectrum& spectrum,
                                              const CTemplate& tpl,
                                              const TFloat64Range& lambdaRange,
                                              const TFloat64List& redshifts,
                                              Float64 overlapThreshold,
                                              std::vector<CMask> additional_spcMasks,
                                              std::string opt_interp,
                                              Int32 opt_extinction=0,
                                              Int32 opt_dustFitting=0,
                                              CPriorHelper::TPriorZEList logpriorze=CPriorHelper::TPriorZEList(),
                                              Bool keepigmism = false,
                                              Float64 FitDustCoeff=-1,
                                              Float64 FitMeiksinIdx=-1);

    const COperatorResult* ExportChi2versusAZ( const CSpectrum& spectrum, const CTemplate& tpl,
                                    const TFloat64Range& lambdaRange, const TFloat64List& redshifts,
                                    Float64 overlapThreshold );


private:

    void BasicFit(const CSpectrum& spectrum,
                  const CTemplate& tpl,
                  const TFloat64Range& lambdaRange,
                  Float64 redshift,
                  Float64 overlapThreshold,
                  Float64& overlapRate,
                  Float64& chiSquare,
                  Float64 &fittingAmplitude,
                  Float64& fittingAmplitudeError,
                  Bool& fittingAmplitudeNegative,
                  Float64& fittingDtM,
                  Float64& fittingMtM,
                  Float64 &fittingLogprior,
                  Float64 &fittingDustCoeff,
                  Float64 &fittingMeiksinIdx,
                  EStatus& status,
                  std::vector<TFloat64List>& ChiSquareInterm,
                  std::vector<TFloat64List>& IsmCalzettiCoeffInterm,
                  std::vector<TInt32List>& IgmMeiksinIdxInterm,
                  std::string opt_interp,
                  Float64 forcedAmplitude=-1,
                  Int32 opt_extinction=0,
                  Int32 opt_dustFitting=0,
                  CMask spcMaskAdditional=CMask(),
                  CPriorHelper::TPriorEList logpriore=CPriorHelper::TPriorEList(),
                  bool keepigmism=false);


    Int32    GetSpcSampleLimits(const TAxisSampleList & Xspc,  Float64 lbda_min, Float64 lbda_max, Int32& kStart, Int32& kEnd);

    //Likelihood
    Float64 EstimateLikelihoodCstLog(const CSpectrum& spectrum, const TFloat64Range& lambdaRange);

    // buffers for the interpolated axis (template & spectrum)




};


}

#endif
