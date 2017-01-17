#ifndef _REDSHIFT_OPERATOR_CHISQUARE_
#define _REDSHIFT_OPERATOR_CHISQUARE_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>
#include <epic/redshift/operator/operator.h>
#include <epic/redshift/operator/chisquareresult.h>

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class COperatorChiSquare : public COperator
{

public:

    COperatorChiSquare();
    ~COperatorChiSquare();

    std::shared_ptr<COperatorResult> Compute(const CSpectrum& spectrum, const CTemplate& tpl,
                                    const TFloat64Range& lambdaRange, const TFloat64List& redshifts,
                                    Float64 overlapThreshold, std::vector<CMask> additional_spcMasks_unused, std::string opt_interp_unused="lin", Int32 opt_extinction_unused=0,
                                              Int32 opt_dustFitting_unused=0 );


private:

    Void BasicFit(const CSpectrum& spectrum, const CTemplate& tpl,
                   const TFloat64Range& lambdaRange, Float64 redshift, Float64 overlapThreshold,
                   Float64& overlapRate, Float64& chiSquare, Float64 &fitamplitude, EStatus& status  );


};


}

#endif
