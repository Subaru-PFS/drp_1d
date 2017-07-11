#ifndef ELEMENT_H
#define ELEMENT_H


#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/common/datatypes.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/ray/catalog.h>

namespace NSEpic
{

  /**
   * /ingroup Redshift
   */
class CLineModelElement
{

public:

    CLineModelElement(const std::string& widthType, const Float64 resolution, const Float64 velocityEmission, const Float64 velocityAbsorption);
    ~CLineModelElement();

    std::string GetElementTypeTag();

    virtual void prepareSupport(const CSpectrumSpectralAxis& spectralAxis, Float64 redshift, const TFloat64Range& lambdaRange)=0;
    virtual TInt32RangeList getSupport()=0;
    virtual TInt32Range getSupportSubElt(Int32 subeIdx)=0;
    virtual TInt32Range getTheoreticalSupportSubElt(Int32 subeIdx)=0;
    virtual Float64 GetContinuumAtCenterProfile(Int32 subeIdx, const CSpectrumSpectralAxis& spectralAxis, Float64 redshift, CSpectrumFluxAxis &continuumfluxAxis)=0;


    virtual void fitAmplitude(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, const CSpectrumFluxAxis &continuumfluxAxis, Float64  redshift, Int32 lineIdx=-1 ) =0;
    virtual Float64 getModelAtLambda( Float64 lambda, Float64 redshift, Float64 continuumFlux, Int32 kRaySupport=-1 )=0;
    virtual Float64 GetModelDerivAmplitudeAtLambda( Float64 lambda, Float64 redshift )=0;
    virtual Float64 GetModelDerivSigmaAtLambda( Float64 lambda, Float64 redshift )=0;
    virtual void addToSpectrumModel( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, CSpectrumFluxAxis &continuumfluxAxis, Float64 redshift, Int32 lineIdx=-1 )=0;
    virtual void addToSpectrumModelDerivSigma( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift )=0;

    virtual void initSpectrumModel( CSpectrumFluxAxis& modelfluxAxis, CSpectrumFluxAxis& continuumfluxAxis, Int32 lineIdx=-1 )=0;

    virtual Float64 GetNominalAmplitude(Int32 subeIdx)=0;
    virtual bool SetNominalAmplitude(Int32 subeIdx, Float64 nominalamp)=0;
    virtual Float64 GetFittedAmplitude(Int32 subeIdx)=0;
    virtual Float64 GetFittedAmplitudeErrorSigma(Int32 subeIdx)=0;
    virtual Float64 GetElementAmplitude()=0;
    virtual void SetFittedAmplitude(Float64 A, Float64 SNR)=0;
    virtual void LimitFittedAmplitude(Int32 subeIdx, Float64 limit)=0;

    virtual bool SetAbsLinesLimit(Float64 limit)=0;

    void SetVelocityEmission(Float64 vel);
    Float64 GetVelocityEmission();
    void SetVelocityAbsorption(Float64 vel);
    Float64 GetVelocityAbsorption();

    void SetAsymfitWidthCoeff(Float64 coeff);
    Float64 GetAsymfitWidthCoeff();
    void SetAsymfitAlphaCoeff(Float64 coeff);
    Float64 GetAsymfitAlphaCoeff();
    void SetAsymfitDelta(Float64 coeff);
    Float64 GetAsymfitDelta();

    virtual Float64 GetSignFactor(Int32 subeIdx)=0;
    virtual Float64 GetWidth(Int32 subeIdx, Float64 redshift)=0;
    Int32 GetSize();
    virtual std::vector<CRay> GetRays()=0;
    virtual std::string GetRayName(Int32 subeIdx)=0;
    bool IsOutsideLambdaRange();
    virtual bool IsOutsideLambdaRange(Int32 subeIdx)=0;
    Int32 FindElementIndex(Int32 LineCatalogIndex);
    virtual Int32 FindElementIndex(std::string LineTagStr)=0;

    std::vector<Int32> m_LineCatalogIndexes;
    Float64 GetLineWidth(Float64 lambda, Float64 z, Bool isEmission, std::string profile);
    Float64 GetLineProfile(std::string profile, Float64 x, Float64 x0, Float64 c);
    Float64 GetLineProfileDerivSigma(std::string profile, Float64 x, Float64 x0, Float64 sigma);
    Float64 GetNSigmaSupport(std::string profile);

    Float64 GetSumCross();
    void SetSumCross(Float64 val);
    Float64 GetSumGauss();
    void SetSumGauss(Float64 val);
    Float64 GetFitAmplitude();

    std::vector<CRay>       m_Rays; //only used in multiline for now... tbd: should be moved elsewhere ?
protected:
    Bool LoadDataExtinction();


    std::string m_LineWidthType;
    Float64 m_NominalWidth;
    Float64 m_Resolution;
    Float64 m_VelocityEmission;
    Float64 m_VelocityAbsorption;
    Float64 m_instrumentResolutionEmpiricalFactor;

    Float64 m_OutsideLambdaRangeOverlapThreshold;
    bool m_OutsideLambdaRange;
    std::string m_ElementType;

    Float64 m_asym_sigma_coeff;
    Float64 m_asym_alpha;
    Float64 m_symxl_sigma_coeff;

    Float64 m_asym2_sigma_coeff;
    Float64 m_asym2_alpha;
    Float64 m_asymfit_sigma_coeff;
    Float64 m_asymfit_alpha;
    Float64 m_asymfit_delta;


    Float64* m_dataExtinctionFlux;
    Float64 m_dataStepLambda = 0.1;
    Float64 m_dataN = 3000;

    Float64 m_sumCross = 0.0;
    Float64 m_sumGauss = 0.0;
    Float64 m_fitAmplitude = 0.0;

private:


};

}







#endif // ELEMENT_H

