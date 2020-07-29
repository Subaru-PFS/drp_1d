#ifndef _REDSHIFT_OPERATOR_CHISQUARE2_
#define _REDSHIFT_OPERATOR_CHISQUARE2_

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
#include <RedshiftLibrary/linemodel/modelspectrumresult.h>
namespace NSEpic
{

class COperatorChiSquare2 : public COperator
{

public:

    explicit COperatorChiSquare2( std::string calibrationPath );
    ~COperatorChiSquare2();

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
    const Float64*  getDustCoeff(Float64 dustCoeff, Float64 maxLambda);
    const Float64*  getMeiksinCoeff(Int32 meiksinIdx, Float64 redshift, Float64 maxLambda);
    Int32           GetSpectrumModel(const CSpectrum& spectrum,
                                      const CTemplate& tpl,
                                      Float64 redshift,
                                      Float64 IdxDustCoeff,
                                      Int32 meiksinIdx,
                                      std::string opt_interp,
                                      std::string opt_extinction,
                                      const TFloat64Range& lambdaRange,
                                      std::shared_ptr<CModelSpectrumResult> resultspcmodel);
    void            SaveSpectrumResults(CDataStore &dataStore);
private:

    void BasicFit_preallocateBuffers(const CSpectrum& spectrum);

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

    Int32  RebinTemplate( const CSpectrum& spectrum,
                          const CTemplate& tpl, 
                          Float64 redshift,
                          const TFloat64Range& lambdaRange,
                          std::string opt_interp,
                          TFloat64Range& currentRange,
                          Float64& overlaprate,
                          Float64 overlapThreshold);// const;

    Int32    GetSpcSampleLimits(const TAxisSampleList & Xspc, TFloat64Range& currentRange, Int32& kStart, Int32& kEnd);
    // buffers for the interpolated axis (template & spectrum)
    CTemplate       m_templateRebined_bf; //buffer
    CMask           m_mskRebined_bf; //buffer
    CSpectrumSpectralAxis m_spcSpectralAxis_restframe; //buffer

    //ISM Calzetti
    Float64* m_YtplRawBuffer;
    Int32 m_YtplRawBufferMaxBufferSize;
    std::vector<std::shared_ptr<CModelSpectrumResult>>  m_savedModelSpectrumResults;

    std::shared_ptr<CSpectrumFluxCorrectionCalzetti> m_ismCorrectionCalzetti;
    //IGM meiksin
    std::shared_ptr<CSpectrumFluxCorrectionMeiksin> m_igmCorrectionMeiksin;

    //Likelihood
    Float64 EstimateLikelihoodCstLog(const CSpectrum& spectrum, const TFloat64Range& lambdaRange);
    void SaveSpectrumResults(CDataStore &dataStore, std::shared_ptr<COperatorResult> result);
};


}

#endif
