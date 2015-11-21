#ifndef _REDSHIFT_OPERATOR_CHISQUARE_
#define _REDSHIFT_OPERATOR_CHISQUARE_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>
#include <epic/redshift/operator/operator.h>
#include <epic/redshift/operator/chisquareresult.h>

namespace NSEpic
{

class COperatorChiSquare : public COperator
{

public:

    COperatorChiSquare();
    ~COperatorChiSquare();

    const COperatorResult* Compute( const CSpectrum& spectrum, const CTemplate& tpl,
                                    const TFloat64Range& lambdaRange, const TFloat64List& redshifts,
                                    Float64 overlapThreshold );


private:

    Void BasicFit(const CSpectrum& spectrum, const CTemplate& tpl,
                   const TFloat64Range& lambdaRange, Float64 redshift, Float64 overlapThreshold,
                   Float64& overlapRate, Float64& chiSquare, Float64 &fitamplitude, EStatus& status  );


};


}

#endif
