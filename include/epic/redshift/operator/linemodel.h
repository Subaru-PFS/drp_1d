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

/**
 * \ingroup Redshift
 */
class COperatorLineModel : public CManagedObject
{
    DEFINE_MANAGED_OBJECT( COperatorLineModel )
public:

    COperatorLineModel();
    virtual ~COperatorLineModel();

    const COperatorResult* Compute(CDataStore &dataStore, const CSpectrum& spectrum, const CSpectrum &spectrumContinuum, const CRayCatalog& restraycatalog,
                                    const TFloat64Range& lambdaRange, const TFloat64List& redshifts , const Int32 opt_extremacount, const std::string &opt_continuumcomponent, const std::string& opt_lineWidthType, const std::string &opt_continuumreest="no");

private:

    Void ModelFit(NSEpic::CLineModelElementList &model, const TFloat64Range& lambdaRange, Float64 redshift,
                  Float64& chiSquare, CLineModelResult::SLineModelSolution &modelSolution, Int32 contreest_iterations);

    void ComputeArea1(CLineModelResult* results);
    void ComputeArea2(CLineModelResult* results);
    Float64 FitBayesWidth( CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis, Float64 z, Int32 start, Int32 end);

    Float64 PrecomputeLogErr(const CSpectrum& spectrum);

    Float64 mSumLogErr;
};


}

#endif
