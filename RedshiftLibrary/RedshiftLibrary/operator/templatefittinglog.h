#ifndef _REDSHIFT_OPERATOR_CHISQUARELOGLAMBDA_
#define _REDSHIFT_OPERATOR_CHISQUARELOGLAMBDA_

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
#include <RedshiftLibrary/operator/modelcontinuumfittingresult.h>

#include <fftw3.h>

namespace NSEpic
{

class COperatorTemplateFittingLog : public COperatorTemplateFittingBase
{

public:
    COperatorTemplateFittingLog() = delete;
    COperatorTemplateFittingLog();
    ~COperatorTemplateFittingLog();

    COperatorTemplateFittingLog(COperatorTemplateFittingLog const& other) = delete;
    COperatorTemplateFittingLog& operator=(COperatorTemplateFittingLog const& other) = delete;  

    std::shared_ptr<COperatorResult> Compute( const CSpectrum& spectrum,
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
                                              Float64 FitEbmvCoeff=-1,
                                              Float64 FitMeiksinIdx=-1);

    void enableSpcLogRebin(Bool enable);
    Int32   ComputeSpectrumModel(const CSpectrum& spectrum,
                                const CTemplate& tpl,
                                Float64 redshift,
                                Float64 DustCoeff,
                                Int32 meiksinIdx,
                                Float64 amplitude,
                                std::string opt_interp,
                                std::string opt_extinction,
                                const TFloat64Range& lambdaRange,
                                Float64 overlapThreshold,
                                CModelSpectrumResult& spc);
        
    Float64 Computelogstep( const CSpectrum &spectrum,
                            const CTemplate &tpl,
                            const TFloat64Range &lambdaRange, 
                            const TFloat64List &sortedRedshifts, 
                            TFloat64Range& newZrange);
    Int32   LoglambdaRebinSpectrum(const CSpectrum &spectrum, const TFloat64Range &lambdaRange, Float64 logGridStep, std::string errorRebinMethod="");
    Int32   CheckLoglambdaRebinSpectrum(const CSpectrum &spectrum, const TFloat64Range &lambdaRange, Float64 logGridStep);
    Int32   LogRebinTemplate(const CTemplate &tpl, 
                            Float64 logGridStep,
                            TFloat64Range& zrange);
    void    BasicFit_preallocateBuffers(const CTemplate & tpl, const Int32 logGridCount);
    //made public for unit-testing
    TInt32Range FindTplSpectralIndex(const TFloat64Range redshiftrange,
                                    const Float64 redshiftStep);
    Int32 computeTargetLogSpectralAxis(const CSpectrumSpectralAxis &ref_axis, 
                                        Float64 logGridStep,
                                        Float64 tgt_loglbdamin,
                                        CSpectrumSpectralAxis& targetSpectralAxis);
private:

    //hardcoded config: REBIN
    Bool verboseLogRebin = 0;
    Bool verboseExportLogRebin = 0;
    const std::string rebinMethod = "lin";

    //hardcoded config: FIT_RANGEZ
    bool verboseLogFitFitRangez = false;
    bool verboseExportFitRangez = false;
    bool verboseExportFitRangez_model = false;
    UInt32 exportIGMIdx = 5;
    UInt32 exportISMIdx = -1;

    //hardcoded config: XTY_FFT
    bool verboseLogXtYFFT = false;
    bool verboseExportXtYFFT = false;

    Int32 FitAllz(const TFloat64Range& lambdaRange,
                  std::shared_ptr<CTemplateFittingResult> result,
                  std::vector<Int32> igmMeiksinCoeffs=std::vector<Int32>(1, 0),
                  std::vector<Int32> ismEbmvCoeffs=std::vector<Int32>(1, 0),
                  CMask spcMaskAdditional=CMask(),
                  CPriorHelper::TPriorZEList logpriorze=CPriorHelper::TPriorZEList());

    Int32 FitRangez(const TFloat64List & inv_err2,
                    TInt32Range& range,
                    std::shared_ptr<CTemplateFittingResult> result,
                    std::vector<Int32> igmMeiksinCoeffs,
                    std::vector<Int32> ismEbmvCoeffs,
                    const Float64& dtd);

    Int32 EstimateXtY(const std::vector<Float64>& X, const std::vector<Float64>& Y,
                      UInt32 nshifts, std::vector<Float64>& XtY, Int32 precomputedFFT=-1);
    Int32 InitFFT(Int32 n);
    Int32 EstimateXtYSlow(const std::vector<Float64>& X, const std::vector<Float64>& Y, UInt32 nShifts,
                          std::vector<Float64>& XtY);
    Int32 EstimateMtMFast(const std::vector<Float64>& X, const std::vector<Float64>& Y, UInt32 nShifts, std::vector<Float64>& XtY);


    Int32 InterpolateResult(const std::vector<Float64>& in, std::vector<Float64>& inGrid,
                            const std::vector<Float64>& tgtGrid, std::vector<Float64>& out, Float64 defaultValue);

    void freeFFTPlans();
    void freeFFTPrecomputedBuffers();

    bool m_opt_spcrebin;

    //log grid data
    CTemplate       m_templateRebinedLog;
    CMask           m_mskRebinedLog;
    CSpectrum       m_spectrumRebinedLog;

    Bool m_enableISM = 1;
    Bool m_enableIGM = 1; 

    //buffers for fft computation
    Int32 m_nPaddedSamples;
    Float64 *inSpc;
    fftw_complex *outSpc;
    fftw_plan pSpc;
    Float64 *inTpl;
    Float64 *inTpl_padded;
    fftw_complex *outTpl;
    fftw_plan pTpl;
    fftw_complex* outCombined;
    Float64* inCombined;
    fftw_plan pBackward ;
    fftw_complex* precomputedFFT_spcFluxOverErr2;
    fftw_complex* precomputedFFT_spcOneOverErr2;

    //Likelihood
    Float64 EstimateLikelihoodCstLog(const CSpectrum& spectrum, const TFloat64Range& lambdaRange);
};


}

#endif
