#ifndef LINEMODEL_ELEMENT_SINGLELINE_H
#define LINEMODEL_ELEMENT_SINGLELINE_H

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
   * Model predicting a single spectral line to be present in the data.
   */
class CSingleLine:public CLineModelElement
{

public:

    CSingleLine(const CRay &r, const std::string& widthType, const Float64 resolution, const Float64 velocityEmission, const Float64 velocityAbsorption, Float64 nominalWidth, std::vector<Int32> catalogIndexes);
    ~CSingleLine();

    std::string GetRayName(Int32 subeIdx);
    Float64 GetSignFactor(Int32 subeIdx);
    Float64 GetWidth(Int32 subeIdx, Float64 redshift);
    std::vector<CRay> GetRays();

    void prepareSupport(const CSpectrumSpectralAxis& spectralAxis, Float64 redshift, const TFloat64Range& lambdaRange);
    TInt32RangeList getSupport();
    TInt32Range getSupportSubElt(Int32 subeIdx);
    TInt32Range getTheoreticalSupportSubElt(Int32 subeIdx);


    void fitAmplitude(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64  redshift, Int32 lineIdx );
    //Float64 FitAmplitudeIterative( const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64 lambda, Float64 width, Int32 start, Int32 end); //deprecated
    Float64 getModelAtLambda( Float64 lambda, Float64 redshift, Int32 kRaySupport=-1 );
    Float64 GetModelDerivAmplitudeAtLambda( Float64 lambda, Float64 redshift );
    Float64 GetModelDerivSigmaAtLambda( Float64 lambda, Float64 redshift );

    void addToSpectrumModel(const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift, Int32 lineIdx=-1 );
    void addToSpectrumModelDerivSigma(const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift );
    void initSpectrumModel( CSpectrumFluxAxis& modelfluxAxis, CSpectrumFluxAxis& continuumfluxAxis, Int32 lineIdx=-1 );
    Float64 GetNominalAmplitude(Int32 subeIdx);
    Float64 GetFittedAmplitude(Int32 subeIdx);
    Float64 GetFittedAmplitudeErrorSigma(Int32 subeIdx);
    Float64 GetElementAmplitude();
    void SetFittedAmplitude(Float64 A, Float64 SNR);
    void LimitFittedAmplitude(Int32 subeIdx, Float64 limit);
    bool IsOutsideLambdaRange(Int32 subeIdx);



private:
    Int32 FindElementIndex(std::string LineTagStr);

    CRay    m_Ray;
    Float64 m_SignFactor;
    Float64 m_FittedAmplitude;
    Float64 m_FittedAmplitudeErrorSigma;

    Int32 m_Start;
    Int32 m_End;

};

}



#endif // ELEMENT_H

