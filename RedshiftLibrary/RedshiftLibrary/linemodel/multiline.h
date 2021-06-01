#ifndef _REDSHIFT_LINEMODEL_MULTILINE_
#define _REDSHIFT_LINEMODEL_MULTILINE_

#include <RedshiftLibrary/linemodel/element.h>
#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/ray/lineprofile.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace NSEpic
{
class CLineProfileASYM;
class CLineProfileASYMFIXED;
  /**
   * \ingroup Redshift
   */
class CMultiLine:public CLineModelElement
{

public:

    CMultiLine(std::vector<CRay> rs, 
               const std::string& widthType,
               const Float64 velocityEmission, 
               const Float64 velocityAbsorption, 
               TFloat64List nominalAmplitudes, 
               Float64 nominalWidth, 
               TUInt32List catalogIndexes);
    ~CMultiLine();

    std::string GetRayName(Int32 subeIdx);
    Float64 GetObservedPosition(Int32 subeIdx, Float64 redshift, Bool doAsymfitdelta=true);
    Float64 GetLineProfileAtRedshift(Int32 subeIdx, Float64 redshift, Float64 x);
    Float64 GetWidth(Int32 subeIdx, Float64 redshift);
    Float64 GetSignFactor(Int32 subeIdx);
    std::vector<CRay> GetRays();

    void prepareSupport(const CSpectrumSpectralAxis& spectralAxis, Float64 redshift, const TFloat64Range& lambdaRange);
    TInt32RangeList getSupport();
    TInt32RangeList getTheoreticalSupport();
    TInt32Range EstimateTheoreticalSupport(Int32 subeIdx, const CSpectrumSpectralAxis& spectralAxis, Float64 redshift,  const TFloat64Range &lambdaRange);
    TInt32Range getSupportSubElt(Int32 subeIdx);
    TInt32Range getTheoreticalSupportSubElt(Int32 subeIdx);

    TInt32Range EstimateIndexRange(const CSpectrumSpectralAxis& spectralAxis, Float64 mu,  const TFloat64Range &lambdaRange, Float64 winsizeAngstrom);

    Float64 GetContinuumAtCenterProfile(Int32 subeIdx, const CSpectrumSpectralAxis& spectralAxis, Float64 redshift, CSpectrumFluxAxis &continuumfluxAxis);

    void fitAmplitude(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, const CSpectrumFluxAxis &continuumfluxAxis, Float64  redshift, Int32 lineIdx=-1 );
    void fitAmplitudeAndLambdaOffset(const CSpectrumSpectralAxis& spectralAxis,
                                     const CSpectrumFluxAxis& fluxAxis,
                                     const CSpectrumFluxAxis &continuumfluxAxis,
                                     Float64  redshift,
                                     Int32 lineIdx=-1 ,
                                     bool enableOffsetFitting=true,
                                     Float64 step=25.,
                                     Float64 min=-400.,
                                     Float64 max=400.);
    Float64 getModelAtLambda(Float64 lambda, Float64 redshift, Float64 continuumFlux, Int32 kRaySupport);
    Float64 GetModelDerivAmplitudeAtLambda( Float64 lambda, Float64 redshift, Float64 continuumFlux  );
    Float64 GetModelDerivContinuumAmpAtLambda(Float64 lambda, Float64 redshift, Float64 continuumFluxUnscale );
    Float64 GetModelDerivZAtLambdaNoContinuum(Float64 lambda, Float64 redshift, Float64 continuumFlux);
    Float64 GetModelDerivZAtLambda(Float64 lambda, Float64 redshift, Float64 continuumFlux,  Float64 continuumFluxDerivZ);


    void addToSpectrumModel( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, CSpectrumFluxAxis &continuumfluxAxis, Float64 redshift, Int32 lineIdx=-1 );
    void addToSpectrumModelDerivVel( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, CSpectrumFluxAxis& continuumFluxAxis, Float64 redshift, bool emissionRay );
    void initSpectrumModel(CSpectrumFluxAxis &modelfluxAxis , const CSpectrumFluxAxis &continuumfluxAxis, Int32 lineIdx=-1 );
    Float64 GetFittedAmplitude(Int32 subeIdx);
    Float64 GetFittedAmplitudeErrorSigma(Int32 subeIdx);
    Float64 GetNominalAmplitude(Int32 subeIdx);
    bool SetNominalAmplitude(Int32 subeIdx, Float64 nominalamp);
    Float64 GetElementAmplitude();
    Float64 GetElementError();
    void SetFittedAmplitude(Int32 subeIdx, Float64 A, Float64 SNR);
    void SetFittedAmplitude(Float64 A, Float64 SNR);
    void LimitFittedAmplitude(Int32 subeIdx, Float64 limit);
    bool IsOutsideLambdaRange(Int32 subeIdx);

    bool SetAbsLinesLimit(Float64 limit);

private:
    Int32 FindElementIndex(std::string LineTagStr);

    std::vector<TInt32List>     m_RayIsActiveOnSupport;
    TFloat64List        m_SignFactors;
    TFloat64List        m_FittedAmplitudes;
    TFloat64List        m_FittedAmplitudeErrorSigmas;
    TFloat64List        m_NominalAmplitudes;

    Float64 m_absLinesLimit;

    TProfileList  m_profile;


    TInt32List          m_StartNoOverlap;
    TInt32List          m_EndNoOverlap;
    TInt32List          m_StartTheoretical;
    TInt32List          m_EndTheoretical;

    TBoolList           m_OutsideLambdaRangeList;

    const bool m_verbose=false;
};

}

#endif // _REDSHIFT_LINEMODEL_MULTILINE_
