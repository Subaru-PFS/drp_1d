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
#ifndef _REDSHIFT_LINEMODEL_ELEMENT_
#define _REDSHIFT_LINEMODEL_ELEMENT_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/ray/lineprofile.h"
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "RedshiftLibrary/ray/catalog.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

namespace NSEpic {

/**
 * \ingroup Redshift
 */
class CLineModelElement
{
  enum TLineWidthType {
    INSTRUMENTDRIVEN,
    COMBINED,
    VELOCITYDRIVEN
  };



  public:
    CLineModelElement(const std::string &widthType,
                      const Float64 velocityEmission,
                      const Float64 velocityAbsorption);
    ~CLineModelElement();

    std::string GetElementTypeTag();

    virtual void prepareSupport(const CSpectrumSpectralAxis &spectralAxis,
                                Float64 redshift,
                                const TFloat64Range &lambdaRange) = 0;
    virtual TInt32RangeList getSupport() = 0;
    virtual TInt32RangeList getTheoreticalSupport() = 0;
    virtual TInt32Range getSupportSubElt(Int32 subeIdx) = 0;
    virtual TInt32Range getTheoreticalSupportSubElt(Int32 subeIdx) = 0;

    virtual TInt32Range
    EstimateIndexRange(const CSpectrumSpectralAxis &spectralAxis,
                       Float64 mu, const TFloat64Range &lambdaRange,
                       Float64 winsizeAngstrom) = 0;

    virtual Float64 GetContinuumAtCenterProfile(
        Int32 subeIdx, const CSpectrumSpectralAxis &spectralAxis,
        Float64 redshift, CSpectrumFluxAxis &continuumfluxAxis) = 0;

    virtual void fitAmplitude(const CSpectrumSpectralAxis &spectralAxis,
                              const CSpectrumFluxAxis &fluxAxis,
                              const CSpectrumFluxAxis &continuumfluxAxis,
                              Float64 redshift, Int32 lineIdx = -1) = 0;
    virtual void fitAmplitudeAndLambdaOffset(
        const CSpectrumSpectralAxis &spectralAxis,
        const CSpectrumFluxAxis &fluxAxis,
        const CSpectrumFluxAxis &continuumfluxAxis, Float64 redshift,
        Int32 lineIdx = -1, bool enableOffsetFitting = true, Float64 step = 25.,
        Float64 min = -400., Float64 max = 400.) = 0;
    virtual Float64 getModelAtLambda(Float64 lambda, Float64 redshift,
                                     Float64 continuumFlux,
                                     Int32 kRaySupport = -1) = 0;
    virtual Float64 GetModelDerivAmplitudeAtLambda(Float64 lambda,
                                                   Float64 redshift,
                                                   Float64 continuumFlux) = 0;
    virtual Float64
    GetModelDerivContinuumAmpAtLambda(Float64 lambda, Float64 redshift,
                                      Float64 continuumFluxUnscale) = 0;
    virtual Float64
    GetModelDerivZAtLambdaNoContinuum(Float64 lambda, Float64 redshift,
                                      Float64 continuumFlux) = 0;
    virtual Float64 GetModelDerivZAtLambda(Float64 lambda, Float64 redshift,
                                           Float64 continuumFlux,
                                           Float64 continuumFluxDerivZ) = 0;

    virtual void
    addToSpectrumModel(const CSpectrumSpectralAxis &modelspectralAxis,
                       CSpectrumFluxAxis &modelfluxAxis,
                       CSpectrumFluxAxis &continuumfluxAxis, Float64 redshift,
                       Int32 lineIdx = -1) = 0;
    virtual void
    addToSpectrumModelDerivVel(const CSpectrumSpectralAxis &modelspectralAxis,
                               CSpectrumFluxAxis &modelfluxAxis,
                               CSpectrumFluxAxis &continuumFluxAxis,
                               Float64 redshift, bool emissionRay) = 0;

    virtual void initSpectrumModel(CSpectrumFluxAxis &modelfluxAxis,
                                   const CSpectrumFluxAxis &continuumfluxAxis,
                                   Int32 lineIdx = -1) = 0;

    virtual Float64 GetNominalAmplitude(Int32 subeIdx) = 0;
    virtual bool SetNominalAmplitude(Int32 subeIdx, Float64 nominalamp) = 0;
    virtual Float64 GetFittedAmplitude(Int32 subeIdx) = 0;
    virtual Float64 GetFittedAmplitudeErrorSigma(Int32 subeIdx) = 0;
    virtual Float64 GetElementAmplitude() = 0;
    virtual Float64 GetElementError() = 0;

    virtual void SetFittedAmplitude(Int32 subeIdx, Float64 A, Float64 SNR) = 0;
    virtual void SetFittedAmplitude(Float64 A, Float64 SNR) = 0;
    virtual void LimitFittedAmplitude(Int32 subeIdx, Float64 limit) = 0;

    virtual bool SetAbsLinesLimit(Float64 limit) = 0;

    void SetVelocityEmission(Float64 vel);
    Float64 GetVelocityEmission();
    void SetVelocityAbsorption(Float64 vel);
    Float64 GetVelocityAbsorption();
    Float64 GetVelocity();
  void setVelocity(Float64 vel);

    void SetSourcesizeDispersion(Float64 sigma);

    void SetLSF(const std::shared_ptr<const CLSF> & lsf);

    void SetAsymfitParams(TAsymParams params, Int32 indx=-99);//-99 means setting for all
    const TAsymParams GetAsymfitParams(UInt32 asymIdx=0);
    void resetAsymfitParams();

    virtual Float64 GetSignFactor(Int32 subeIdx) = 0;
    virtual Float64 GetObservedPosition(Int32 subeIdx, Float64 redshift, Bool doAsymfitdelta=true) = 0;
    virtual Float64 GetWidth(Int32 subeIdx, Float64 redshift) = 0;
    Int32 GetSize();
    virtual std::vector<CRay> GetRays() = 0;
    virtual std::string GetRayName(Int32 subeIdx) = 0;
    bool IsOutsideLambdaRange();
    virtual bool IsOutsideLambdaRange(Int32 subeIdx) = 0;
    Int32 FindElementIndex(Int32 LineCatalogIndex);
    virtual Int32 FindElementIndex(std::string LineTagStr) = 0;

    std::vector<Int32> m_LineCatalogIndexes;
    Float64 GetLineWidth(Float64 lambda, Float64 z = 0., Bool isEmission=0);
    Float64 GetLineProfileDerivVel(std::shared_ptr<CLineProfile>& profile, Float64 x, Float64 x0,
                                   Float64 sigma, Bool isEmission);

    Float64 GetSumCross();
    void SetSumCross(Float64 val);
    Float64 GetDtmFree();
    void SetDtmFree(Float64 val);
    Float64 GetSumGauss();
    void SetSumGauss(Float64 val);
    Float64 GetFitAmplitude();

    std::vector<CRay> m_Rays; // only used in multiline for now... tbd: should
                              // be moved elsewhere ?
    std::string m_fittingGroupInfo;

  protected:
    Bool LoadDataExtinction();

    TLineWidthType m_LineWidthType;
    Float64 m_NominalWidth;//relevant only for LSF GaussianConstantWidth

    Float64 m_VelocityEmission;
    Float64 m_VelocityAbsorption;

    Float64 m_OutsideLambdaRangeOverlapThreshold;
    bool m_OutsideLambdaRange;
    std::string m_ElementType;

    TUInt32List         m_asymLineIndices;//corresponds to indices of asymmetric lines, mainly LyA. Currently max 1 asymfit is found per linecatalog

    Float64 *m_dataExtinctionFlux = NULL;
    Float64 m_dataStepLambda = 0.1;
    Float64 m_dataN = 3000.0;

    Float64 m_sumCross = 0.0;
    Float64 m_sumGauss = 0.0;
    Float64 m_dtmFree =  0.0; // dtmFree is the non-positive-constrained version of sumCross
    Float64 m_fitAmplitude = 0.0;

    // Constant
    const Float64 m_speedOfLightInVacuum = GSL_CONST_MKSA_SPEED_OF_LIGHT / 1000.0; // km.s^-1

    std::shared_ptr<const CLSF> m_LSF;

  private:
};

} // namespace NSEpic

#endif // _REDSHIFT_LINEMODEL_ELEMENT_
