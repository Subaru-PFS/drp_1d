#ifndef LINEMODEL_ELEMENT_MULTILINE_H
#define LINEMODEL_ELEMENT_MULTILINE_H

#include <epic/redshift/linemodel/element.h>
#include <epic/redshift/ray/catalog.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/core/common/range.h>
#include <epic/redshift/common/datatypes.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace NSEpic
{

  /**
   * \ingroup Redshift
   */
class CMultiLine:public CLineModelElement
{

public:

    CMultiLine(std::vector<CRay> rs, const std::string& widthType, const Float64 resolution, const Float64 velocityEmission, const Float64 velocityAbsorption, std::vector<Float64> nominalAmplitudes, Float64 nominalWidth, std::vector<Int32> catalogIndexes);
    ~CMultiLine();

    std::string GetRayName(Int32 subeIdx);
    Float64 GetWidth(Int32 subeIdx, Float64 redshift);
    Float64 GetSignFactor(Int32 subeIdx);
    std::vector<CRay> GetRays();

    void prepareSupport(const CSpectrumSpectralAxis& spectralAxis, Float64 redshift, const TFloat64Range& lambdaRange);
    TInt32RangeList getSupport();
    TInt32Range getSupportSubElt(Int32 subeIdx);

    void fitAmplitude(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64  redshift);
    Float64 getModelAtLambda(Float64 lambda, Float64 redshift , Int32 kRaySupport);
    Float64 GetModelDerivAmplitudeAtLambda( Float64 lambda, Float64 redshift );
    Float64 GetModelDerivSigmaAtLambda( Float64 lambda, Float64 redshift );

    void addToSpectrumModel( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift );
    void addToSpectrumModelDerivSigma( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift );
    void initSpectrumModel(CSpectrumFluxAxis &modelfluxAxis , CSpectrumFluxAxis &continuumfluxAxis);
    Float64 GetFittedAmplitude(Int32 subeIdx);
    Float64 GetFittedAmplitudeErrorSigma(Int32 subeIdx);
    Float64 GetNominalAmplitude(Int32 subeIdx);
    Float64 GetElementAmplitude();
    void SetFittedAmplitude(Float64 A, Float64 SNR);
    void LimitFittedAmplitude(Int32 subeIdx, Float64 limit);
    bool IsOutsideLambdaRange(Int32 subeIdx);

private:
    Int32 FindElementIndex(std::string LineTagStr);

    std::vector<CRay>       m_Rays;
    std::vector<std::vector<Int32>>     m_RayIsActiveOnSupport;
    std::vector<Float64>    m_SignFactors;
    std::vector<Float64>        m_FittedAmplitudes;
    std::vector<Float64>        m_FittedAmplitudeErrorSigmas;
    std::vector<Float64>        m_NominalAmplitudes;

    std::vector<Int32>          m_Start;
    std::vector<Int32>          m_End;

    std::vector<bool>           m_OutsideLambdaRangeList;
};

}


#endif // ELEMENT_H

