#ifndef _REDSHIFT_OPERATOR_CORRELATION_
#define _REDSHIFT_OPERATOR_CORRELATION_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/spectrum/fluxaxis.h>
#include <RedshiftLibrary/spectrum/spectralaxis.h>
#include <RedshiftLibrary/operator/operator.h>
#include <RedshiftLibrary/operator/correlationresult.h>
#include <RedshiftLibrary/statistics/priorhelper.h>
#include <RedshiftLibrary/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplate;

/**
 * \ingroup Redshift
 */
class COperatorCorrelation : public COperator
{

public:

    COperatorCorrelation();
    ~COperatorCorrelation();

    std::shared_ptr<COperatorResult> Compute( const CSpectrum& s1,
                                              const CTemplate& s2,
                                              const TFloat64Range& r,
                                              const TFloat64List& redhisfts,
                                              Float64 overlap,
                                              std::vector<CMask> additional_spcMasks_unused,
                                              std::string opt_interp_unused="lin",
                                              Int32 opt_extinction_unused=0 ,
                                              Int32 opt_dustFitting_unused=0,
                                              CPriorHelper::TPriorZEList logpriorze_unused=CPriorHelper::TPriorZEList(),
                                              Bool keepigmism = false,
                                              Float64 FitDustCoeff=-1,
                                              Float64 FitMeiksinIdx=-1);

    Float64 GetComputationDuration() const;

private:

    Float64                 m_TotalDuration;

    // buffers for the interpolated axis (template & spectrum)
    CTemplate       m_templateRebined_bf; //buffer
    CMask           m_mskRebined_bf; //buffer
    CSpectrumSpectralAxis m_spcSpectralAxis_restframe; //buffer

};


}

#endif
