#ifndef _REDSHIFT_OPERATOR_TPLCOMBINATION_
#define _REDSHIFT_OPERATOR_TPLCOMBINATION_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/operator/operator.h"

#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/operator/modelcontinuumfittingresult.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/processflow/resultstore.h"

#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"
#include <gsl/gsl_matrix_double.h>
namespace NSEpic
{
class CSpectrum;
class COperatorResult;
class CModelSpectrumResult;

class COperatorTplcombination
{
public:

    std::shared_ptr<COperatorResult> Compute(const CSpectrum& spectrum,
                                             const TTemplateConstRefList& tplList,
                                             const TFloat64Range& lambdaRange,
                                             const TFloat64List& redshifts,
                                             Float64 overlapThreshold,
                                             std::vector<CMask> additional_spcMasks, 
                                             std::string opt_interp, 
                                             Int32 opt_extinction=0, 
                                             Int32 opt_dustFitting=0,
                                             CPriorHelper::TPriorZEList logpriorze=CPriorHelper::TPriorZEList(),
                                             Bool keepigmism = false,
                                             Float64 FitEbmvCoeff=-1,
                                             Float64 FitMeiksinIdx=-1);

    Float64 ComputeDtD(const CSpectrumFluxAxis& spcFluxAxis, const TInt32Range& range); //could be also made static
    Int32   ComputeSpectrumModel(   const CSpectrum& spectrum,
                                    const TTemplateConstRefList& tplList,
                                    Float64 redshift,
                                    Float64 EbmvCoeff,
                                    Int32 meiksinIdx,
                                    const TFloat64List& amplitudes,
                                    std::string opt_interp,
                                    const TFloat64Range& lambdaRange,
                                    Float64 overlapThreshold,
                                    std::shared_ptr<CModelSpectrumResult> & spcPtr);
private:

    struct STplcombination_basicfitresult
    {
        COperator::EStatus status;
        Float64     overlapRate;
        Float64     chisquare;
        CSpectrum   modelSpectrum;
        TFloat64List    fittingAmplitudes;
        TFloat64List    fittingAmplitudeErrors;
        TFloat64List    fittingAmplitudeSigmas;
        std::vector<TFloat64List>    ChiSquareInterm;
        std::vector<TFloat64List>    IsmCalzettiCoeffInterm;
        std::vector<TInt32List>      IgmMeiksinIdxInterm;
        std::vector<std::vector<TFloat64List>>    fittingAmplitudesInterm; //intermediate amplitudes
        std::vector<std::string> tplNames; //cause combination of templates
        Int32 IGMIdx;
        Float64 EbmvCoeff;
        Float64 SNR;
        std::vector<TFloat64List> COV; 
    };

    void BasicFit_preallocateBuffers(const CSpectrum& spectrum,
                                     const TTemplateConstRefList& tplList);

    void BasicFit(const CSpectrum& spectrum,
                  const TTemplateConstRefList& tplList,
                  const TFloat64Range& lambdaRange,
                  Float64 redshift,
                  Float64 overlapThreshold,
                  STplcombination_basicfitresult& fittingResults,
                  std::string opt_interp, Float64 forcedAmplitude=-1, 
                  Int32 opt_extinction=0, 
                  Int32 opt_dustFitting=0, 
                  CMask spcMaskAdditional=CMask(),
                  CPriorHelper::TPriorEList logpriore=CPriorHelper::TPriorEList(),
                  bool keepigmism=false,
                  const TInt32List& MeiksinList=TInt32List(-1));
    Int32 RebinTemplate(const CSpectrum& spectrum,
                        const TTemplateConstRefList& tplList,
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

    Float64 ComputeXi2_bruteForce(const CSpectrumFluxAxis& correctedFlux, 
                                  const CSpectrumFluxAxis& spcFluxAxis,
                                  const Int32 imin_lbda);
    Float64 GetNormFactor(const CSpectrumFluxAxis spcFluxAxis, UInt32 kStart, UInt32 n);
};


}

#endif
