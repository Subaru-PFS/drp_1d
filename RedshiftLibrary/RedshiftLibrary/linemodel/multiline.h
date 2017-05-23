#ifndef LINEMODEL_ELEMENT_MULTILINE_H
#define LINEMODEL_ELEMENT_MULTILINE_H

#include <RedshiftLibrary/linemodel/element.h>
#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/common/datatypes.h>

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
    TInt32Range EstimateTheoreticalSupport(Int32 subeIdx, const CSpectrumSpectralAxis& spectralAxis, Float64 redshift,  const TFloat64Range &lambdaRange);
    TInt32Range getSupportSubElt(Int32 subeIdx);
    TInt32Range getTheoreticalSupportSubElt(Int32 subeIdx);

    void fitAmplitude(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, const CSpectrumFluxAxis &continuumfluxAxis, Float64  redshift, Int32 lineIdx=-1 );
    Float64 getModelAtLambda(Float64 lambda, Float64 redshift, Float64 continuumFlux, Int32 kRaySupport);
    Float64 GetModelDerivAmplitudeAtLambda( Float64 lambda, Float64 redshift );
    Float64 GetModelDerivSigmaAtLambda( Float64 lambda, Float64 redshift );

    void addToSpectrumModel( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, CSpectrumFluxAxis &continuumfluxAxis, Float64 redshift, Int32 lineIdx=-1 );
    void addToSpectrumModelDerivSigma( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift );
    void initSpectrumModel(CSpectrumFluxAxis &modelfluxAxis , CSpectrumFluxAxis &continuumfluxAxis, Int32 lineIdx=-1 );
    Float64 GetFittedAmplitude(Int32 subeIdx);
    Float64 GetFittedAmplitudeErrorSigma(Int32 subeIdx);
    Float64 GetNominalAmplitude(Int32 subeIdx);
    bool SetNominalAmplitude(Int32 subeIdx, Float64 nominalamp);
    Float64 GetElementAmplitude();
    void SetFittedAmplitude(Float64 A, Float64 SNR);
    void LimitFittedAmplitude(Int32 subeIdx, Float64 limit);
    bool IsOutsideLambdaRange(Int32 subeIdx);

private:
    Int32 FindElementIndex(std::string LineTagStr);

    std::vector<std::vector<Int32>>     m_RayIsActiveOnSupport;
    std::vector<Float64>    m_SignFactors;
    std::vector<Float64>        m_FittedAmplitudes;
    std::vector<Float64>        m_FittedAmplitudeErrorSigmas;
    std::vector<Float64>        m_NominalAmplitudes;

    std::vector<Float64>        mBuffer_mu;
    std::vector<Float64>        mBuffer_c;
    std::vector<std::string>    m_profile;


    std::vector<Int32>          m_StartNoOverlap;
    std::vector<Int32>          m_EndNoOverlap;
    std::vector<Int32>          m_StartTheoretical;
    std::vector<Int32>          m_EndTheoretical;

    std::vector<bool>           m_OutsideLambdaRangeList;

    //constant
    Float64 m_c_kms;
};

}


#endif // ELEMENT_H

