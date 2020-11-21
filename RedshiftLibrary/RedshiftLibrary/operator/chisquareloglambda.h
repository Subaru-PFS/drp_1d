#ifndef _REDSHIFT_OPERATOR_CHISQUARELOGLAMBDA_
#define _REDSHIFT_OPERATOR_CHISQUARELOGLAMBDA_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/operator/operator.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <RedshiftLibrary/common/mask.h>

#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h>
#include <RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h>
#include <RedshiftLibrary/statistics/priorhelper.h>

#include <fftw3.h>

namespace NSEpic
{

class COperatorChiSquareLogLambda : public COperator
{

public:

    COperatorChiSquareLogLambda(std::string calibrationPath);
    ~COperatorChiSquareLogLambda();

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
                                              Float64 FitDustCoeff=-1,
                                              Float64 FitMeiksinIdx=-1);

    void enableSpcLogRebin(Bool enable);

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
                  std::shared_ptr<CChisquareResult> result,
                  std::vector<Int32> igmMeiksinCoeffs=std::vector<Int32>(1, 0),
                  std::vector<Int32> ismEbmvCoeffs=std::vector<Int32>(1, 0),
                  CMask spcMaskAdditional=CMask(),
                  CPriorHelper::TPriorZEList logpriorze=CPriorHelper::TPriorZEList());

    Int32 FitRangez(const TAxisSampleList &  spectrumRebinedLambda,
                    const TAxisSampleList &  spectrumRebinedFluxRaw,
                    const TAxisSampleList &  error,
                    const TAxisSampleList &  tplRebinedLambda,
                    const TAxisSampleList &  tplRebinedFluxRaw,
                    std::shared_ptr<CChisquareResult> result,
                    std::vector<Int32> igmMeiksinCoeffs,
                    std::vector<Int32> ismEbmvCoeffs);

    Int32 EstimateXtY(const std::vector<Float64>& X, const std::vector<Float64>& Y,
                      UInt32 nshifts, std::vector<Float64>& XtY, Int32 precomputedFFT=-1);
    Int32 InitFFT(Int32 n);
    Int32 EstimateXtYSlow(const std::vector<Float64>& X, const std::vector<Float64>& Y, UInt32 nShifts,
                          std::vector<Float64>& XtY);
    Int32 EstimateMtMFast(const std::vector<Float64>& X, const std::vector<Float64>& Y, UInt32 nShifts, std::vector<Float64>& XtY);

    TInt32Range FindTplSpectralIndex(const TAxisSampleList & spcLambda,
                                     const TAxisSampleList & tplLambda,
                                     TFloat64Range redshiftrange, Float64 redshiftStep);

    Int32 InterpolateResult(const std::vector<Float64>& in, std::vector<Float64>& inGrid,
                            const std::vector<Float64>& tgtGrid, std::vector<Float64>& out, Float64 defaultValue);

    void freeFFTPlans();
    void freeFFTPrecomputedBuffers();

    bool m_opt_spcrebin;

    //log grid data
    CTemplate       m_templateRebinedLog;
    CMask           m_mskRebinedLog;
    CSpectrum       m_spectrumRebinedLog;

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


    //ISM Calzetti
    std::unique_ptr<CSpectrumFluxCorrectionCalzetti> m_ismCorrectionCalzetti;

    //IGM meiksin
    std::unique_ptr<CSpectrumFluxCorrectionMeiksin> m_igmCorrectionMeiksin;

    //Likelihood
    Float64 EstimateLikelihoodCstLog(const CSpectrum& spectrum, const TFloat64Range& lambdaRange);
};


}

#endif
