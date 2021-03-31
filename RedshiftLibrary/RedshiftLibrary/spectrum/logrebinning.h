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
    void  Computelogstep( CSpectrum &spectrum,
                        const TFloat64Range &lambdaRange, 
                        const Float64 zInputStep,
                        const TFloat64Range zInputRange);
    //todo: merge these two functions once all is stable:
    //pass shared_ptr
    std::shared_ptr<const CSpectrum>   LoglambdaRebinSpectrum(std::shared_ptr<CSpectrum> spectrum, 
                                                            const TFloat64Range &lambdaRange, 
                                                            std::string errorRebinMethod="rebinVariance");
    std::shared_ptr<CTemplate>         LoglambdaRebinTemplate(std::shared_ptr<CTemplate>  tpl);

    CSpectrumSpectralAxis  computeTargetLogSpectralAxis(TFloat64Range lambdarange,
                                                        UInt32 gridCount);
    Bool verifyLogRebinningResults(const CSpectrumSpectralAxis& spcAxis, const CSpectrumSpectralAxis& tplAxis);
    TFloat64Range m_zrange;
    Float64 m_logGridStep; 
    TFloat64Range m_lambdaRange_spc, m_lambdaRange_tpl;
private:
    void InferTemplateRebinningSetup(TFloat64Range lambdaRange);
    const std::string m_rebinMethod = "lin";

    UInt32 m_loglambda_count_spc, m_loglambda_count_tpl;
};

}
#endif