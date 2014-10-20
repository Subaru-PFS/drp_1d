#ifndef _REDSHIFT_OPERATOR_CHISQUARE_
#define _REDSHIFT_OPERATOR_CHISQUARE_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>
#include <epic/redshift/operator/operator.h>

#include <vector>

namespace NSEpic
{

class CSpectrum;
class CTemplate;
class CRedshifts;

class COperatorChiSquare : public COperator
{

public:

    COperatorChiSquare();
    ~COperatorChiSquare();

    Bool Compute( const CSpectrum& spectrum, const CTemplate& tpl,
                  const TFloat64Range& lambdaRange, const CRedshifts& redshifts, Float64 overlapThreshold );

    const TFloat64List& GetResults() const;

private:

    Float64 BasicFit( const CSpectrum& spectrum, const CTemplate& tpl,
                   const TFloat64Range& lambdaRange, Float64 redshift, Float64 overlapThreshold,
                   Float64& ampl, Float64& overlapRate, Int32& numDevs );

    TFloat64List    m_Chisquare;
    TFloat64List    m_Overlap;
};


}

#endif
