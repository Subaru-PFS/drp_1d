#ifndef _REDSHIFT_SPECTRUM_LOGREBINNING_
#define _REDSHIFT_SPECTRUM_LOGREBINNING_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/common/indexing.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>



namespace NSEpic
{

class CSpectrumLogRebinning 
{

public:
    //applying Rule of zero
    TFloat64Range&  Computelogstep( const CSpectrum &spectrum,
                                    const TFloat64Range &lambdaRange, 
                                    const Float64 zInputStep,
                                    const TFloat64Range zInputRange);
    std::shared_ptr<const CSpectrum>   LoglambdaRebinSpectrum(const CSpectrum &spectrum, const TFloat64Range &lambdaRange, std::string errorRebinMethod="rebinVariance");
    std::shared_ptr<const CSpectrum>   CheckLoglambdaSampling(const CSpectrum &spectrum);
    std::shared_ptr<CTemplate>         LoglambdaRebinTemplate(const CTemplate &tpl);

    CSpectrumSpectralAxis  computeTargetLogSpectralAxis(const CSpectrumSpectralAxis &ref_axis,
                                                        Float64 tgt_loglbdamin,
                                                        Float64 gridCount);
    Float64 GetLogGridStep();
private:
    void InferTemplateRebinningSetup(TFloat64Range rebinnedSpectrumLambdaRange);
    const std::string m_rebinMethod = "lin";
    CMask           m_mskRebinedLog;

    Float64 m_logGridStep, m_loglambdaReference; 
    UInt32 m_spectrumLogGridCount; //spectrum grid count
    UInt32 m_templateLogGridCount; //template grid count
    TFloat64Range m_zrange;
    TFloat64Range m_tplLambdaRange; //template range
};
}
#endif