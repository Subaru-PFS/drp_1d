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
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/line/lineprofile.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

namespace NSEpic {

/**
 * \ingroup Redshift
 */
class CLineModelElement {
  enum TLineWidthType { INSTRUMENTDRIVEN, COMBINED, VELOCITYDRIVEN };

public:
  CLineModelElement(std::vector<CLine> rs, const std::string &widthType,
                    const Float64 velocityEmission,
                    const Float64 velocityAbsorption,
                    TFloat64List nominalAmplitudes, Float64 nominalWidth,
                    TInt32List catalogIndexes);

  Float64 GetObservedPosition(Int32 subeIdx, Float64 redshift,
                              bool doAsymfitdelta = true) const;
  Float64 GetLineProfileAtRedshift(Int32 subeIdx, Float64 redshift,
                                   Float64 x) const;
  void getObservedPositionAndLineWidth(Int32 subeIdx, Float64 redshift,
                                       Float64 &mu, Float64 &sigma,
                                       bool doAsymfitdelta = true) const;

  const std::string &GetElementTypeTag() const;

  void prepareSupport(const CSpectrumSpectralAxis &spectralAxis,
                      Float64 redshift, const TFloat64Range &lambdaRange);
  TInt32RangeList getSupport() const;
  TInt32RangeList getTheoreticalSupport() const;
  void EstimateTheoreticalSupport(Int32 subeIdx,
                                  const CSpectrumSpectralAxis &spectralAxis,
                                  Float64 redshift,
                                  const TFloat64Range &lambdaRange);
  TInt32Range getSupportSubElt(Int32 subeIdx) const;
  TInt32Range getTheoreticalSupportSubElt(Int32 subeIdx) const;

  static TInt32Range
  EstimateIndexRange(const CSpectrumSpectralAxis &spectralAxis, Float64 mu,
                     const TFloat64Range &lambdaRange, Float64 winsizeAngstrom);

  Float64 GetContinuumAtCenterProfile(
      Int32 subeIdx, const CSpectrumSpectralAxis &spectralAxis,
      Float64 redshift, const CSpectrumFluxAxis &continuumfluxAxis,
      const TPolynomCoeffs &line_polynomCoeffs) const;

  void fitAmplitude(const CSpectrumSpectralAxis &spectralAxis,
                    const CSpectrumFluxAxis &fluxAxis,
                    const CSpectrumFluxAxis &continuumfluxAxis,
                    Float64 redshift, Int32 lineIdx = undefIdx);
  void fitAmplitudeAndLambdaOffset(const CSpectrumSpectralAxis &spectralAxis,
                                   const CSpectrumFluxAxis &fluxAxis,
                                   const CSpectrumFluxAxis &continuumfluxAxis,
                                   Float64 redshift, Int32 lineIdx = undefIdx,
                                   bool enableOffsetFitting = true,
                                   Float64 step = 25., Float64 min = -400.,
                                   Float64 max = 400.);
  Float64 getModelAtLambda(Float64 lambda, Float64 redshift,
                           Float64 continuumFlux,
                           Int32 kLineSupport = -1) const;
  Float64 GetModelDerivAmplitudeAtLambda(Float64 lambda, Float64 redshift,
                                         Float64 continuumFlux) const;
  Float64 GetModelDerivContinuumAmpAtLambda(Float64 lambda, Float64 redshift,
                                            Float64 continuumFluxUnscale) const;
  Float64 GetModelDerivZAtLambdaNoContinuum(Float64 lambda, Float64 redshift,
                                            Float64 continuumFlux) const;
  Float64 GetModelDerivZAtLambda(Float64 lambda, Float64 redshift,
                                 Float64 continuumFlux,
                                 Float64 continuumFluxDerivZ) const;

  void addToSpectrumModel(const CSpectrumSpectralAxis &modelspectralAxis,
                          CSpectrumFluxAxis &modelfluxAxis,
                          const CSpectrumFluxAxis &continuumfluxAxis,
                          Float64 redshift, Int32 lineIdx = undefIdx) const;
  void
  addToSpectrumModelDerivVel(const CSpectrumSpectralAxis &modelspectralAxis,
                             CSpectrumFluxAxis &modelfluxAxis,
                             const CSpectrumFluxAxis &continuumFluxAxis,
                             Float64 redshift, bool emissionLine) const;

  void initSpectrumModel(CSpectrumFluxAxis &modelfluxAxis,
                         const CSpectrumFluxAxis &continuumfluxAxis,
                         Int32 lineIdx = undefIdx) const;

  Float64 GetNominalAmplitude(Int32 subeIdx) const;
  bool SetNominalAmplitude(Int32 subeIdx, Float64 nominalamp);
  Float64 GetFittedAmplitude(Int32 subeIdx) const;
  Float64 GetFittedAmplitudeErrorSigma(Int32 subeIdx) const;
  Float64 GetElementAmplitude() const;
  Float64 GetElementError() const;

  void SetFittedAmplitude(Int32 subeIdx, Float64 A, Float64 SNR);
  void SetFittedAmplitude(Float64 A, Float64 SNR);
  void LimitFittedAmplitude(Int32 subeIdx, Float64 limit);

  bool SetAbsLinesLimit(Float64 limit);

  void SetVelocityEmission(Float64 vel);
  Float64 GetVelocityEmission() const;
  void SetVelocityAbsorption(Float64 vel);
  Float64 GetVelocityAbsorption() const;
  Float64 GetVelocity() const;
  void setVelocity(Float64 vel);

  void SetLSF(const std::shared_ptr<const CLSF> &lsf);

  void SetAsymfitParams(TAsymParams params,
                        Int32 indx = -99); //-99 means setting for all
  const TAsymParams GetAsymfitParams(Int32 asymIdx = 0) const;
  void resetAsymfitParams();
  Int32 findElementIndex(Int32 LineCatalogIndex) const;
  Int32 findElementIndex(const std::string &LineTagStr) const;
  const CLineProfile &getLineProfile(Int32 lineIdx) const;

  Float64 GetSignFactor(Int32 subeIdx) const;

  Int32 GetSize() const;
  std::vector<CLine> GetLines() const;
  const std::string &GetLineName(Int32 subeIdx) const;
  bool IsOutsideLambdaRange() const;
  Float64 GetLineWidth(Float64 lambda, Float64 z = 0.,
                       bool isEmission = 0) const;
  bool IsOutsideLambdaRange(Int32 subeIdx) const;

  Float64 GetLineProfileDerivVel(const CLineProfile &profile, Float64 x,
                                 Float64 x0, Float64 sigma,
                                 bool isEmission) const;

  Float64 GetSumCross() const;
  void SetSumCross(Float64 val);
  Float64 GetDtmFree() const;
  void SetDtmFree(Float64 val);
  Float64 GetSumGauss() const;
  void SetSumGauss(Float64 val);
  Float64 GetFitAmplitude() const;

  void debug(std::ostream &os) const;

  std::vector<CLine> m_Lines; // only used in multiline for now... tbd: should
                              // be moved elsewhere ?
  TInt32List m_LineCatalogIndexes;
  std::string m_fittingGroupInfo;

protected:
  TLineWidthType m_LineWidthType;
  Float64 m_NominalWidth; // relevant only for LSF GaussianConstantWidth

  Float64 m_VelocityEmission;
  Float64 m_VelocityAbsorption;

  Float64 m_OutsideLambdaRangeOverlapThreshold;
  bool m_OutsideLambdaRange;
  std::string m_ElementType;

  TInt32List m_asymLineIndices; // corresponds to indices of asymmetric lines,
                                // mainly LyA. Currently max 1 asymfit is found
                                // per linecatalog

  Float64 *m_dataExtinctionFlux = NULL;
  Float64 m_dataStepLambda = 0.1;
  Float64 m_dataN = 3000.0;

  Float64 m_sumCross = 0.0;
  Float64 m_sumGauss = 0.0;
  Float64 m_dtmFree =
      0.0; // dtmFree is the non-positive-constrained version of sumCross
  Float64 m_fitAmplitude = 0.0;
  const Float64 m_speedOfLightInVacuum = SPEED_OF_LIGHT_IN_VACCUM;
  std::shared_ptr<const CLSF> m_LSF;

  std::vector<TInt32List> m_LineIsActiveOnSupport;
  TFloat64List m_SignFactors;
  TFloat64List m_FittedAmplitudes;
  TFloat64List m_FittedAmplitudeErrorSigmas;
  TFloat64List m_NominalAmplitudes;

  Float64 m_absLinesLimit;

  TInt32List m_StartNoOverlap;
  TInt32List m_EndNoOverlap;
  TInt32List m_StartTheoretical;
  TInt32List m_EndTheoretical;

  TBoolList m_OutsideLambdaRangeList;

  const bool m_verbose = false;
};

} // namespace NSEpic

#endif // _REDSHIFT_LINEMODEL_ELEMENT_
