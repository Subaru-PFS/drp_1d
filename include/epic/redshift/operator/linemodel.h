#ifndef _REDSHIFT_OPERATOR_LINEMODEL_
#define _REDSHIFT_OPERATOR_LINEMODEL_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>
#include <epic/redshift/operator/operator.h>
#include <epic/redshift/operator/linemodelresult.h>
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

    const COperatorResult* Compute( const CSpectrum& spectrum, const CRayCatalog& restraycatalog,
                                    const TFloat64Range& lambdaRange, const TFloat64List& redshifts );

private:

    Void ModelFit(const CSpectrum& spectrum, CSpectrum &model, const CRayCatalog::TRayVector &restraycatalog,
                   const TFloat64Range& lambdaRange, Float64 redshift,
                  Float64& chiSquare, CLineModelResult::SLineModelSolution &modelSolution);

    Float64 FitAmplitude( const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64 lambda, Float64 width, Int32 start, Int32 end);
    Float64 FitAmplitudeIterative( const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64 lambda, Float64 width, Int32 start, Int32 end);
    Void AddGaussLine( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 lambda, Float64 width, Float64 amplitude, Int32 start, Int32 end);
    Void Apply2LinesAmplitudeRule(const CRayCatalog::TRayVector& restRayList,
                                                      std::vector<Float64>& Amplitudes,
                                                      std::vector<Bool> outsidelambdarange,
                                                      std::string lineA, std::string lineB, Float64 coeff );

    void ComputeGaussAreaForExtrema(CLineModelResult* results);
    Float64 FitBayesWidth( CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis, Float64 z, Int32 start, Int32 end);
};


}

#endif
