#ifndef _REDSHIFT_OPERATOR_TPLCOMBINATION_
#define _REDSHIFT_OPERATOR_TPLCOMBINATION_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/operator/operator.h>
#include <RedshiftLibrary/operator/templatefittingresult.h>
#include <RedshiftLibrary/operator/modelspectrumresult.h>
#include <RedshiftLibrary/operator/modelcontinuumfittingresult.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/processflow/resultstore.h>

#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h>
#include <RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h>

namespace NSEpic
{

class COperatorTplcombination
{
public:

    std::shared_ptr<COperatorResult> Compute(const CSpectrum& spectrum,
                                             const std::vector<CTemplate> &tplList,
                                             const TFloat64Range& lambdaRange,
                                             const TFloat64List& redshifts,
                                             Float64 overlapThreshold,
                                             std::vector<CMask> additional_spcMasks, 
                                             std::string opt_interp, 
                                             Int32 opt_extinction=0, 
                                             Int32 opt_dustFitting=0);

  void SaveSpectrumResults(std::shared_ptr<COperatorResultStore> resultStore);

private:

    struct STplcombination_basicfitresult
    {
        COperator::EStatus status;
        Float64     overlapRate;
        Float64     chisquare;
        CSpectrum/*CModelSpectrumResult*/   modelSpectrum;
        TFloat64List    fittingAmplitudes;
        TFloat64List    fittingErrors;
        std::vector<TFloat64List> ChiSquareInterm;
        std::vector<std::string> tplNames; //cause combination of templates
        Int32 igmIdx;
        Float64 ebmvCoeff;
        Float64 snr;
    };

    std::vector<std::shared_ptr<CModelSpectrumResult>  > m_savedModelSpectrumResults;
    std::vector<std::shared_ptr<CModelContinuumFittingResult> > m_savedModelContinuumFittingResults;
    void BasicFit_preallocateBuffers(const CSpectrum& spectrum,
                                     const std::vector<CTemplate>& tplList);

    void BasicFit(const CSpectrum& spectrum,
                  const std::vector<CTemplate>& tplList,
                  const TFloat64Range& lambdaRange,
                  Float64 redshift,
                  Float64 overlapThreshold,
                  STplcombination_basicfitresult& fittingResults,
                  std::string opt_interp, Float64 forcedAmplitude=-1, Int32 opt_extinction=0, Int32 opt_dustFitting=0, CMask spcMaskAdditional=CMask() );
    Int32 RebinTemplate(const CSpectrum& spectrum,
                        const std::vector<CTemplate>& tplList,
                        Float64 redshift,
                        const TFloat64Range& lambdaRange,
                        std::string opt_interp,
                        TFloat64Range& currentRange,
                        Float64& overlapRate,
                        Float64 overlapThreshold);
    // buffers for the interpolated axis (templates & spectrum)
    std::vector<CTemplate>   m_templatesRebined_bf; //vector of buffer
    std::vector<CMask>       m_masksRebined_bf; //vector of buffer
    CSpectrumSpectralAxis    m_spcSpectralAxis_restframe; //buffer

    //Likelihood
    Float64 EstimateLikelihoodCstLog(const CSpectrum& spectrum, const TFloat64Range& lambdaRange);

};


}

#endif
