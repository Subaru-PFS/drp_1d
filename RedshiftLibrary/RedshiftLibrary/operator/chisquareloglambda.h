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

#include <fftw3.h>

namespace NSEpic
{

class COperatorChiSquareLogLambda : public COperator
{

public:

    COperatorChiSquareLogLambda( std::string calibrationPath );
    ~COperatorChiSquareLogLambda();

     std::shared_ptr<COperatorResult> Compute(const CSpectrum& spectrum, const CTemplate& tpl,
                                    const TFloat64Range& lambdaRange, const TFloat64List& redshifts,
                                    Float64 overlapThreshold, std::vector<CMask> additional_spcMasks, std::string opt_interp, Int32 opt_extinction=0, Int32 opt_dustFitting=0);

    const Float64*  getDustCoeff(Float64 dustCoeff, Float64 maxLambda);
    const Float64*  getMeiksinCoeff(Int32 meiksinIdx, Float64 redshift, Float64 maxLambda);


private:

    Int32 FitAllz(const TFloat64Range& lambdaRange,
                  std::shared_ptr<CChisquareResult> result,
                  std::vector<Int32> igmMeiksinCoeffs=std::vector<Int32>(1, 0),
                  std::vector<Int32> ismEbmvCoeffs=std::vector<Int32>(1, 0),
                  CMask spcMaskAdditional=CMask());
    Int32 EstimateXtY(const Float64 *X, const Float64 *Y, UInt32 nx, UInt32 ny, UInt32 nshifts, std::vector<Float64>& XtY);
    Int32 InitFFT(Int32 n);
    Int32 EstimateXtYSlow(const Float64* X, const Float64* Y, UInt32 nX, UInt32 nShifts, std::vector<Float64>& XtY);
    Int32 EstimateMtMFast(const Float64* X, const Float64* Y, UInt32 nX, UInt32 nShifts, std::vector<Float64>& XtY);

    Int32 InterpolateResult(const Float64* in, const Float64* inGrid, const Float64* tgtGrid, Int32 n, Int32 tgtn, std::vector<Float64>& out, Float64 defaultValue);


    //log grid template
    CTemplate       m_templateRebinedLog;
    CMask           m_mskRebinedLog;
    CSpectrum       m_spectrumRebinedLog;
    CSpectrumFluxAxis m_errorRebinedLog;

    //buffers for fft computation
    Int32 m_nPaddedSamples;
    fftw_complex *inSpc;
    fftw_complex *outSpc;
    fftw_plan pSpc;
    fftw_complex *inTpl;
    fftw_complex *inTpl_padded;
    fftw_complex *outTpl;
    fftw_plan pTpl;
    fftw_complex* outCombined;
    fftw_complex* inCombined;
    fftw_plan pBackward ;


    //ISM Calzetti
    CSpectrumFluxCorrectionCalzetti* m_ismCorrectionCalzetti;

    //IGM meiksin
    CSpectrumFluxCorrectionMeiksin* m_igmCorrectionMeiksin;
};


}

#endif
