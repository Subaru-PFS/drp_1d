#ifndef _REDSHIFT_OPERATOR_CORRELATION_
#define _REDSHIFT_OPERATOR_CORRELATION_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>
#include <epic/redshift/spectrum/fluxaxis.h>
#include <epic/redshift/spectrum/spectralaxis.h>
#include <epic/redshift/operator/operator.h>
#include <epic/redshift/operator/correlationresult.h>


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
                                              Int32 opt_extinction_unused=0 );

    Float64 GetComputationDuration() const;

private:

    Float64                 m_TotalDuration;
};


}

#endif
