#ifndef _REDSHIFT_OPERATOR_CHISQUARE2_
#define _REDSHIFT_OPERATOR_CHISQUARE2_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>
#include <epic/redshift/operator/operator.h>
#include <epic/redshift/operator/chisquareresult.h>
#include <epic/redshift/common/mask.h>

#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/template/template.h>

namespace NSEpic
{

class COperatorChiSquare2 : public COperator
{

public:

    COperatorChiSquare2();
    ~COperatorChiSquare2();

    const COperatorResult* Compute( const CSpectrum& spectrum, const CTemplate& tpl,
                                    const TFloat64Range& lambdaRange, const TFloat64List& redshifts,
                                    Float64 overlapThreshold );


private:

    Void BasicFit(const CSpectrum& spectrum, const CTemplate& tpl, CTemplate &tplFineBuffer,
                   const TFloat64Range& lambdaRange, Float64 redshift, Float64 overlapThreshold,
                   Float64& overlapRate, Float64& chiSquare, EStatus& status  );

    CTemplate       m_templateRebined_bf; //buffer
    CMask           m_mskRebined_bf; //buffer
    CSpectrumSpectralAxis m_shiftedTplSpectralAxis_bf; //buffer

};


}

#endif
