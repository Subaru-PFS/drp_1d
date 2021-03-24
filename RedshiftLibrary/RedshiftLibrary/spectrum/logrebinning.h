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
    std::shared_ptr<const CSpectrum>   LoglambdaRebinSpectrum(std::shared_ptr<CSpectrum> spectrum, 
                                                            const TFloat64Range &lambdaRange, 
                                                            std::string errorRebinMethod="rebinVariance");
    std::shared_ptr<CTemplate>         LoglambdaRebinTemplate(const CTemplate &tpl);

    CSpectrumSpectralAxis  computeTargetLogSpectralAxis(TFloat64Range lambdarange,
                                                        UInt32 gridCount);
    Float64 GetLogGridStep();
    TFloat64Range GetRedshiftRange();
private:
    void InferTemplateRebinningSetup(TFloat64Range lambdaRange);
    const std::string m_rebinMethod = "lin";

    Float64 m_logGridStep; 
    Float64 m_log_lambda_start; 

    UInt32 m_loglambda_count, m_loglambda_count_ref;

    TFloat64Range m_loglambdaRange, m_loglambdaRange_ref;
    TFloat64Range m_zrange;
};
inline
Float64 CSpectrumLogRebinning::GetLogGridStep()
{
    return m_logGridStep;
}
inline
TFloat64Range CSpectrumLogRebinning::GetRedshiftRange()
{
    return m_zrange;
}
}
#endif