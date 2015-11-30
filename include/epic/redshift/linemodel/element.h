#ifndef ELEMENT_H
#define ELEMENT_H


#include <epic/core/common/range.h>
#include <epic/redshift/common/datatypes.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/ray/catalog.h>

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

    virtual void fitAmplitude(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64  redshift) =0;
    virtual Float64 getModelAtLambda( Float64 lambda, Float64 redshift )=0;
    virtual void addToSpectrumModel( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift )=0;
    virtual void initSpectrumModel( CSpectrumFluxAxis& modelfluxAxis, CSpectrumFluxAxis& continuumfluxAxis )=0;

    virtual Float64 GetNominalAmplitude(Int32 subeIdx)=0;
    virtual Float64 GetFittedAmplitude(Int32 subeIdx)=0;
    virtual Float64 GetFittedAmplitudeErrorSigma(Int32 subeIdx)=0;
    virtual Float64 GetElementAmplitude()=0;
    virtual void SetFittedAmplitude(Float64 A, Float64 SNR)=0;
    virtual void LimitFittedAmplitude(Int32 subeIdx, Float64 limit)=0;

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
    Float64 GetLineWidth(Float64 lambda, Float64 z, Bool isEmission);

    Float64 GetNSigmaSupport();

protected:


    std::string m_LineWidthType;
    Float64 m_NominalWidth;
    Float64 m_Resolution;
    Float64 m_VelocityEmission;
    Float64 m_VelocityAbsorption;
    Float64 m_FWHM_factor;

    Float64 m_NSigmaSupport;
    Float64 m_OutsideLambdaRangeOverlapThreshold;
    bool m_OutsideLambdaRange;
    std::string m_ElementType;

private:



};

}







#endif // ELEMENT_H

