#ifndef _REDSHIFT_SPECTRUM_LOGREBINNING_
#define _REDSHIFT_SPECTRUM_LOGREBINNING_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/common/indexing.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/processflow/inputcontext.h>

namespace NSEpic
{

class CSpectrumLogRebinning 
{

public:
    //applying Rule of zero
    void  RebinInputs(CInputContext& inputContext); 
    TFloat64Range m_zrange;
    Float64 m_logGridStep; 
    TFloat64Range m_lambdaRange_spc, m_lambdaRange_tpl;
private:
    void  SetupRebinning( CSpectrum &spectrum,
                        const TFloat64Range &lambdaRange, 
                        Float64 zInputStep,
                        const TFloat64Range & zInputRange,
                        UInt32 SSratio);
    std::shared_ptr<CSpectrum>   LoglambdaRebinSpectrum(std::shared_ptr<const CSpectrum> spectrum,
                                                            std::string errorRebinMethod="rebinVariance");
    std::shared_ptr<CTemplate>         LoglambdaRebinTemplate(std::shared_ptr<const CTemplate> tpl);

    CSpectrumSpectralAxis  computeTargetLogSpectralAxis(TFloat64Range lambdarange, UInt32 gridCount);
    void InferTemplateRebinningSetup(const TFloat64Range & lambdaRange_Ref);
    const std::string m_rebinMethod = "lin";

    UInt32 m_loglambda_count_spc, m_loglambda_count_tpl;
};

}
#endif