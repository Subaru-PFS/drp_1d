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
    COperatorTemplateFittingLog();
    ~COperatorTemplateFittingLog();

    COperatorTemplateFittingLog(const COperatorTemplateFittingLog & other) = delete; 
    COperatorTemplateFittingLog(COperatorTemplateFittingLog && other) = delete; 
    COperatorTemplateFittingLog& operator=(const COperatorTemplateFittingLog& other) = delete;  
    COperatorTemplateFittingLog& operator=(COperatorTemplateFittingLog&& other) = delete; 

    std::shared_ptr<COperatorResult> Compute( const CSpectrum& rebinnedSpectrum,
                                              const CTemplate& rebinnedTpl,
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

    inline  bool IsFFTProcessing() override{return true;}; 

    //made public for unit-testing
    TInt32Range FindTplSpectralIndex(const TFloat64Range & redshiftrange) const;
    TInt32Range FindTplSpectralIndex( const CSpectrumSpectralAxis & spcSpectrailAxis, 
                                      const CSpectrumSpectralAxis& tplSpectralAxis,
                                      const TFloat64Range & redshiftrange) const;
    //log grid data
    CTemplate       m_templateRebinedLog;
    CSpectrum       m_spectrumRebinedLog;

private:
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

};


}

#endif
