// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#ifndef _REDSHIFT_LINEMODEL_MULTILINE_
#define _REDSHIFT_LINEMODEL_MULTILINE_

#include "RedshiftLibrary/linemodel/element.h"
#include "RedshiftLibrary/ray/catalog.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/ray/lineprofile.h"
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
