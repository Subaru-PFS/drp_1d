#ifndef _REDSHIFT_OPERATOR_OPERATOR_
#define _REDSHIFT_OPERATOR_OPERATOR_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>

#include <vector>

namespace __NS__
{

class CSpectrum;
class CTemplate;
class CRedshifts;

class COperator
{


public:

    COperator();
    virtual ~COperator();

    virtual Bool Compute( const CSpectrum& spectrum, const CTemplate& tpl,
                          const TFloat64Range& lambdaRange, const CRedshifts& redshifts, Float64 overlapThreshold ) = 0;

    virtual const TFloat64List& GetResults() const = 0;

private:



};


}

#endif
