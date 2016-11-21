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

     std::shared_ptr<COperatorResult> Compute(const CSpectrum& spectrum, const CTemplate& tpl,
                                    const TFloat64Range& lambdaRange, const TFloat64List& redshifts,
                                    Float64 overlapThreshold, std::vector<CMask> additional_spcMasks, std::string opt_interp, Int32 opt_extinction=0);

    const COperatorResult* ExportChi2versusAZ( const CSpectrum& spectrum, const CTemplate& tpl,
                                    const TFloat64Range& lambdaRange, const TFloat64List& redshifts,
                                    Float64 overlapThreshold );



private:

    Void BasicFit(const CSpectrum& spectrum, const CTemplate& tpl, Float64 *pfgTplBuffer,
                   const TFloat64Range& lambdaRange, Float64 redshift, Float64 overlapThreshold,
                   Float64& overlapRate, Float64& chiSquare, Float64 &fittingAmplitude, Float64 &fittingDustCoeff, EStatus& status, std::string opt_interp, Float64 forcedAmplitude=-1, Int32 opt_extinction=0, CMask spcMaskAdditional=CMask() );

    // buffers for the precomputed fine grid template
    CTemplate       m_templateRebined_bf; //buffer
    CMask           m_mskRebined_bf; //buffer
    CSpectrumSpectralAxis m_shiftedTplSpectralAxis_bf; //buffer

    Float64 *m_dataCalzetti;
    
};


}

#endif
