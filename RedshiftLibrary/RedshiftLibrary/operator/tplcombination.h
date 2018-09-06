#ifndef _REDSHIFT_OPERATOR_TPLCOMBINATION_
#define _REDSHIFT_OPERATOR_TPLCOMBINATION_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/operator/operator.h>
#include <RedshiftLibrary/operator/tplcombinationresult.h>
#include <RedshiftLibrary/linemodel/modelspectrumresult.h>
#include <RedshiftLibrary/common/mask.h>

#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h>
#include <RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h>

namespace NSEpic
{

class COperatorTplcombination
{
public:

    COperatorTplcombination( std::string calibrationPath );
    ~COperatorTplcombination();

    std::shared_ptr<COperatorResult> Compute(const CSpectrum& spectrum,
                                             std::vector<CTemplate> tplList,
                                             const TFloat64Range& lambdaRange,
                                             const TFloat64List& redshifts,
                                             Float64 overlapThreshold,
                                             std::vector<CMask> additional_spcMasks, std::string opt_interp, Int32 opt_extinction=0, Int32 opt_dustFitting=0);

    const Float64*  getDustCoeff(Float64 dustCoeff, Float64 maxLambda);
    const Float64*  getMeiksinCoeff(Int32 meiksinIdx, Float64 redshift, Float64 maxLambda);

    void SaveSpectrumResults(CDataStore &dataStore);

private:

    struct STplcombination_basicfitresult
    {
        COperator::EStatus status;
        Float64     overlapRate;
        Float64     chisquare;
        CSpectrum   modelSpectrum;
        std::vector<Float64>    fittingAmplitudes;
        std::vector<Float64>    fittingErrors;
        std::vector<TFloat64List> ChiSquareInterm;
    };

    std::vector<std::shared_ptr<CModelSpectrumResult>  > m_savedModelSpectrumResults;

    void BasicFit(const CSpectrum& spectrum,
                  std::vector<CTemplate>& tplList,
                  Float64 *pfgTplBuffer,
                  const TFloat64Range& lambdaRange,
                  Float64 redshift,
                  Float64 overlapThreshold,
                  STplcombination_basicfitresult& fittingResults,
                  std::string opt_interp, Float64 forcedAmplitude=-1, Int32 opt_extinction=0, Int32 opt_dustFitting=0, CMask spcMaskAdditional=CMask() );

    // buffers for the precomputed fine grid templates
    std::vector<std::shared_ptr<CTemplate>>  m_templatesRebined_bf; //buffer
    std::vector<std::shared_ptr<CMask>>      m_masksRebined_bf; //buffer
    std::vector<std::shared_ptr<CSpectrumSpectralAxis>>   m_shiftedTemplatesSpectralAxis_bf; //buffer

    //ISM Calzetti
    Float64* m_YtplRawBuffer;
    Int32 m_YtplRawBufferMaxBufferSize;

    CSpectrumFluxCorrectionCalzetti* m_ismCorrectionCalzetti;

    //IGM meiksin
    CSpectrumFluxCorrectionMeiksin* m_igmCorrectionMeiksin;

    //Likelihood
    Float64 EstimateLikelihoodCstLog(const CSpectrum& spectrum, const TFloat64Range& lambdaRange);
};


}

#endif
