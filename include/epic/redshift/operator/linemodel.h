#ifndef _REDSHIFT_OPERATOR_LINEMODEL_
#define _REDSHIFT_OPERATOR_LINEMODEL_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>
#include <epic/redshift/operator/operator.h>
#include <epic/redshift/operator/linemodelresult.h>
#include <epic/redshift/linemodel/elementlist.h>
#include <epic/redshift/common/mask.h>

#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/ray/catalog.h>

namespace NSEpic
{

class COperatorLineModel : public CManagedObject
{
    DEFINE_MANAGED_OBJECT( COperatorLineModel )
public:

    COperatorLineModel();
    virtual ~COperatorLineModel();

    const COperatorResult* Compute(const CSpectrum& spectrum, const CSpectrum &spectrumContinuum, const CRayCatalog& restraycatalog,
                                    const TFloat64Range& lambdaRange, const TFloat64List& redshifts , Int32 lineWidthType);

private:

    Void ModelFit(const CSpectrum& spectrum, NSEpic::CLineModelElementList &model, const CRayCatalog::TRayVector &restraycatalog,
                   const TFloat64Range& lambdaRange, Float64 redshift,
                  Float64& chiSquare, CLineModelResult::SLineModelSolution &modelSolution);

    void ComputeArea1(CLineModelResult* results);
    void ComputeArea2(CLineModelResult* results);
    Float64 FitBayesWidth( CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis, Float64 z, Int32 start, Int32 end);

    Float64 PrecomputeLogErr(const CSpectrum& spectrum);

    Float64 mSumLogErr;
};


}

#endif
