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
#include "RedshiftLibrary/common/polynom.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/line/lineprofile.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace NSEpic {

struct TLineModelElementParam {

  TLineModelElementParam(TLineVector lines, Float64 velocityEmission,
                         Float64 velocityAbsorption,
                         TInt32List lineCatalogIndexes);

  Float64 m_VelocityEmission = NAN;
  Float64 m_VelocityAbsorption = NAN;
  TFloat64List m_FittedAmplitudes;
  TFloat64List m_FittedAmplitudeErrorSigmas;
  TFloat64List m_NominalAmplitudes;
  TFloat64List m_Offsets;
  TLineVector m_Lines;
  TInt32List m_LineCatalogIndexes;
  std::string m_fittingGroupInfo;
  TPolynomCoeffs m_ampOffsetsCoeffs;
};

using TLineModelElementParam_ptr = std::shared_ptr<TLineModelElementParam>;

/**
 * \ingroup Redshift
 */
class CLineModelElement {
  enum TLineWidthType { INSTRUMENTDRIVEN, COMBINED, VELOCITYDRIVEN };

public:
  CLineModelElement(const TLineModelElementParam_ptr elementParam,
                    const std::string &widthType);

  void reset();

  Float64 GetObservedPosition(Int32 subeIdx, Float64 redshift,
                              bool doAsymfitdelta = true) const;
  Float64 GetLineProfileAtRedshift(Int32 subeIdx, Float64 redshift,
                                   Float64 x) const;
  void getObservedPositionAndLineWidth(Int32 subeIdx, Float64 redshift,
                                       Float64 &mu, Float64 &sigma,
                                       bool doAsymfitdelta = true) const;
  Int32 GetElementType() const { return m_type; };
  bool GetIsEmission() const { return m_isEmission; };
  void prepareSupport(const CSpectrumSpectralAxis &spectralAxis,
                      Float64 redshift, const TFloat64Range &lambdaRange);
  TInt32RangeList getSupport() const;
  TInt32RangeList getTheoreticalSupport() const;
  void EstimateTheoreticalSupport(Int32 subeIdx,
                                  const CSpectrumSpectralAxis &spectralAxis,
                                  Float64 redshift,
                                  const TFloat64Range &lambdaRange);
  void SetOutsideLambdaRange();

  TInt32Range getSupportSubElt(Int32 subeIdx) const;
  TInt32Range getTheoreticalSupportSubElt(Int32 subeIdx) const;

  static TInt32Range
  EstimateIndexRange(const CSpectrumSpectralAxis &spectralAxis, Float64 mu,
                     const TFloat64Range &lambdaRange, Float64 winsizeAngstrom);

  Float64 GetContinuumAtCenterProfile(
      Int32 subeIdx, const CSpectrumSpectralAxis &spectralAxis,
      Float64 redshift, const CSpectrumFluxAxis &continuumfluxAxis,
      bool enableAmplitudeOffsets = false) const;

  Float64 getModelAtLambda(Float64 lambda, Float64 redshift,
                           Float64 continuumFlux,
                           Int32 kLineSupport = -1) const;
  Float64 GetModelDerivAmplitudeAtLambda(Float64 lambda, Float64 redshift,
                                         Float64 continuumFlux) const;
  Float64 GetModelDerivVelAtLambda(Float64 lambda, Float64 redshift,
                                   Float64 continuumFlux) const;
  Float64 GetModelDerivContinuumAmpAtLambda(Float64 lambda, Float64 redshift,
                                            Float64 continuumFluxUnscale) const;
  Float64 GetModelDerivZAtLambda(Float64 lambda, Float64 redshift,
                                 Float64 continuumFlux,
                                 Float64 continuumFluxDerivZ = 0.0) const;

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
  void initSpectrumModelPolynomial(CSpectrumFluxAxis &modelfluxAxis,
                                   const CSpectrumSpectralAxis &spcAxis,
                                   Int32 lineIdx) const;

  Float64 GetNominalAmplitude(Int32 subeIdx) const;
  Float64 GetMaxNominalAmplitude() const;
  bool SetNominalAmplitude(Int32 subeIdx, Float64 nominalamp);
  Float64 GetFittedAmplitude(Int32 subeIdx) const;
  Float64 GetFittedAmplitudeErrorSigma(Int32 subeIdx) const;
  Float64 GetElementAmplitude() const;
  Float64 GetElementError() const;

  void SetFittedAmplitude(Int32 subeIdx, Float64 A, Float64 SNR);
  void SetElementAmplitude(Float64 A, Float64 SNR);
  void LimitFittedAmplitude(Int32 subeIdx, Float64 limit);

  bool SetAbsLinesLimit(Float64 limit);
  Float64 GetAbsLinesLimit() const;

  void SetVelocityEmission(Float64 vel);
  Float64 getVelocityEmission() const;
  void SetVelocityAbsorption(Float64 vel);
  Float64 getVelocityAbsorption() const;
  Float64 getVelocity() const;
  void setVelocity(Float64 vel);

  void SetLSF(const std::shared_ptr<const CLSF> &lsf);

  void
  SetAsymfitParams(const TAsymParams &params,
                   Int32 indx = undefIdx); // undefIdx means setting for all
  void SetSymIgmParams(const TSymIgmParams &params,
                       Int32 indx = undefIdx); // undefIdx means setting for all
  void SetSymIgmFit(bool val = true,
                    Int32 indx = undefIdx); // undefIdx means setting for all
  TAsymParams GetAsymfitParams(Int32 asymIdx = 0) const;
  TSymIgmParams GetSymIgmParams(Int32 asymIdx = 0) const;

  void resetAsymfitParams();
  Int32 findElementIndex(Int32 LineCatalogIndex) const;
  Int32 findElementIndex(const std::string &LineTagStr) const;
  Int32 getLineIndexInCatalog(Int32 idxLine, const TLineVector &catalog) const;

  const TInt32List &getIgmLinesIndices() const { return m_asymLineIndices; };
  const CLineProfile_ptr &getLineProfile(Int32 lineIdx) const;

  Float64 GetSignFactor(Int32 subeIdx) const;

  Int32 GetSize() const;
  const TLineVector &GetLines() const { return m_ElementParam->m_Lines; };

  const std::string &GetLineName(Int32 subeIdx) const;
  bool IsOutsideLambdaRange() const;
  Float64 GetLineWidth(Float64 lambda, bool isEmission = 0) const;
  bool IsOutsideLambdaRange(Int32 subeIdx) const;
  void SetOutsideLambdaRangeList(Int32 subeIdx);

  Float64 GetLineProfileDerivVel(const CLineProfile &profile, Float64 x,
                                 Float64 x0, Float64 sigma,
                                 bool isEmission) const;

  Float64 GetSumCross() const;
  void SetSumCross(Float64 val);
  Float64 GetDtmFree() const;
  void SetDtmFree(Float64 val);
  Float64 GetSumGauss() const;
  void SetSumGauss(Float64 val);
  const std::string &GetFittingGroupInfo() const;
  void SetFittingGroupInfo(const std::string &val);

  const TPolynomCoeffs &GetPolynomCoeffs() const;
  void SetPolynomCoeffs(TPolynomCoeffs pCoeffs);

  void SetAllOffsetsEnabled(Float64 val);
  void SetOffset(Int32 lineIdx, Float64 val);
  Float64 GetOffset(Int32 lineIdx) const {
    return m_ElementParam->m_Offsets[lineIdx];
  }

  void SetLineProfile(Int32 lineIdx, CLineProfile_ptr &&profile);

  bool isLineActiveOnSupport(Int32 line, Int32 lineIdx) const;
  Int32 getStartNoOverlap(Int32 lineIdx) const;
  Int32 getEndNoOverlap(Int32 lineIdx) const;
  Int32 getSignFactor(Int32 lineIdx) const;

  void debug(std::ostream &os) const;
  void dumpElement(std::ostream &os) const;

  void computeCrossProducts(Float64 redshift,
                            const CSpectrumSpectralAxis &spectralAxis,
                            const CSpectrumFluxAxis &noContinuumfluxAxis,
                            const CSpectrumFluxAxis &continuumfluxAxis,
                            Int32 lineIdx);

  void fitAmplitude(Float64 redshift, const CSpectrumSpectralAxis &spectralAxis,
                    const CSpectrumFluxAxis &noContinuumfluxAxis,
                    const CSpectrumFluxAxis &continuumfluxAxis,
                    Int32 lineIdx = undefIdx);

protected:
  const TLineModelElementParam_ptr m_ElementParam;

  TLineWidthType m_LineWidthType;

  Float64 m_OutsideLambdaRangeOverlapThreshold;
  bool m_OutsideLambdaRange;

  Float64 m_sumCross = 0.0;
  Float64 m_sumGauss = 0.0;
  Float64 m_dtmFree =
      0.0; // dtmFree is the non-positive-constrained version of sumCross

  const Float64 m_speedOfLightInVacuum = SPEED_OF_LIGHT_IN_VACCUM;
  std::shared_ptr<const CLSF> m_LSF;

  std::vector<TInt32List> m_LineIsActiveOnSupport;
  TFloat64List m_SignFactors;

  TInt32List m_asymLineIndices;

  Float64 m_absLinesLimit;

  TInt32List m_StartNoOverlap;
  TInt32List m_EndNoOverlap;
  TInt32List m_StartTheoretical;
  TInt32List m_EndTheoretical;

  TBoolList m_OutsideLambdaRangeList;
  Int32 m_size;
  Int32 m_type;
  bool m_isEmission;
};

inline bool CLineModelElement::IsOutsideLambdaRange() const {
  return m_OutsideLambdaRange;
}
/**
 * \brief Returns whether the line with index subeIdx is outside the lambda
 *range.
 **/
inline bool CLineModelElement::IsOutsideLambdaRange(
    Int32 subeIdx) const // IsOutsideLambdaRange
{
  return m_OutsideLambdaRangeList[subeIdx];
}

inline bool CLineModelElement::isLineActiveOnSupport(Int32 lineIdxA,
                                                     Int32 lineIdxB) const {
  return m_LineIsActiveOnSupport[lineIdxA][lineIdxB] == 1;
}

/**
 * \brief Returns the nominal amplitude with index subeIdx.
 **/
inline Float64 CLineModelElement::GetNominalAmplitude(Int32 subeIdx) const {
  return m_ElementParam->m_NominalAmplitudes[subeIdx];
}

/**
 * \brief Set the nominal amplitude with index subeIdx.
 **/
inline bool CLineModelElement::SetNominalAmplitude(Int32 subeIdx,
                                                   Float64 nominalamp) {
  m_ElementParam->m_NominalAmplitudes[subeIdx] = nominalamp;
  return true;
}

inline Int32 CLineModelElement::getStartNoOverlap(Int32 lineIdx) const {
  return m_StartNoOverlap[lineIdx];
}
inline Int32 CLineModelElement::getEndNoOverlap(Int32 lineIdx) const {
  return m_EndNoOverlap[lineIdx];
}
inline Int32 CLineModelElement::getSignFactor(Int32 lineIdx) const {
  return m_SignFactors[lineIdx];
}
inline Int32 CLineModelElement::GetSize() const { return m_size; }

/**
 * \brief Returns the fitted amplitude for the argument index.
 **/
inline Float64 CLineModelElement::GetFittedAmplitude(Int32 subeIdx) const {
  return m_ElementParam->m_FittedAmplitudes[subeIdx];
}

/**
 * \brief Returns the fitted amplitude error for the argument index.
 **/
inline Float64
CLineModelElement::GetFittedAmplitudeErrorSigma(Int32 subeIdx) const {
  return m_ElementParam->m_FittedAmplitudeErrorSigmas[subeIdx];
}

inline void CLineModelElement::SetVelocityEmission(Float64 vel) {
  m_ElementParam->m_VelocityEmission = vel;
}

inline void CLineModelElement::SetVelocityAbsorption(Float64 vel) {
  m_ElementParam->m_VelocityAbsorption = vel;
}

inline Float64 CLineModelElement::getVelocityEmission() const {
  return m_ElementParam->m_VelocityEmission;
}

inline Float64 CLineModelElement::getVelocityAbsorption() const {
  return m_ElementParam->m_VelocityAbsorption;
}

inline void CLineModelElement::SetLSF(const std::shared_ptr<const CLSF> &lsf) {
  m_LSF = lsf;
}

inline const CLineProfile_ptr &
CLineModelElement::getLineProfile(Int32 lineIdx) const {
#ifdef DEBUG
  if (lineIdx > GetSize() - 1)
    THROWG(INTERNAL_ERROR, "out-of-bound index");
#endif
  return m_ElementParam->m_Lines[lineIdx].GetProfile();
}

inline void CLineModelElement::setVelocity(Float64 vel) {
#ifdef DEBUG
  if (!GetSize())
    THROWG(INTERNAL_ERROR, "Empty line model element, could not set velocity");
#endif

  if (GetIsEmission()) {
    m_ElementParam->m_VelocityEmission = vel;
  } else {
    m_ElementParam->m_VelocityAbsorption = vel;
  }
}

// wrapper function
inline void CLineModelElement::SetAsymfitParams(const TAsymParams &params,
                                                Int32 idx) {
  if (!m_asymLineIndices.size())
    return;
  if (idx >= 0)
    m_ElementParam->m_Lines[idx].SetAsymParams(params);
  else
    for (auto i : m_asymLineIndices)
      m_ElementParam->m_Lines[i].SetAsymParams(params);
}

// wrapper function
inline void CLineModelElement::SetSymIgmParams(const TSymIgmParams &params,
                                               Int32 idx) {
  if (!m_asymLineIndices.size())
    return;
  if (idx >= 0)
    m_ElementParam->m_Lines[idx].SetSymIgmParams(params);
  else
    for (auto i : m_asymLineIndices)
      m_ElementParam->m_Lines[i].SetSymIgmParams(params);
}

// wrapper function
inline void CLineModelElement::SetSymIgmFit(bool val, Int32 idx) {
  if (!m_asymLineIndices.size())
    return;
  if (idx >= 0)
    m_ElementParam->m_Lines[idx].SetSymIgmFit(val);
  else
    for (auto i : m_asymLineIndices)
      m_ElementParam->m_Lines[i].SetSymIgmFit(val);
}

// wrapper function
inline void CLineModelElement::resetAsymfitParams() {
  for (auto i : m_asymLineIndices)
    m_ElementParam->m_Lines[i].resetAsymFitParams();
}

// wrapper function
inline TAsymParams CLineModelElement::GetAsymfitParams(Int32 idx) const {
  if (!m_asymLineIndices.size())
    return TAsymParams(); // case where no asymprofile in linecatalog
  return m_ElementParam->m_Lines[m_asymLineIndices[idx]].GetAsymParams();
}

// wrapper function
inline TSymIgmParams CLineModelElement::GetSymIgmParams(Int32 idx) const {
  if (!m_asymLineIndices.size())
    return TSymIgmParams(); // case where no asymprofile in linecatalog
  return m_ElementParam->m_Lines[m_asymLineIndices[idx]].GetSymIgmParams();
}

inline Float64 CLineModelElement::GetSumCross() const { return m_sumCross; }

inline void CLineModelElement::SetSumCross(Float64 val) { m_sumCross = val; }

inline Float64 CLineModelElement::GetDtmFree() const { return m_dtmFree; }

inline void CLineModelElement::SetDtmFree(Float64 val) { m_dtmFree = val; }

inline Float64 CLineModelElement::GetSumGauss() const { return m_sumGauss; }

inline void CLineModelElement::SetSumGauss(Float64 val) { m_sumGauss = val; }

inline const std::string &CLineModelElement::GetFittingGroupInfo() const {
  return m_ElementParam->m_fittingGroupInfo;
}

inline void CLineModelElement::SetFittingGroupInfo(const std::string &val) {
  m_ElementParam->m_fittingGroupInfo = val;
}

inline const TPolynomCoeffs &CLineModelElement::GetPolynomCoeffs() const {
  return m_ElementParam->m_ampOffsetsCoeffs;
}

inline void CLineModelElement::SetPolynomCoeffs(TPolynomCoeffs pCoeffs) {
  m_ElementParam->m_ampOffsetsCoeffs = std::move(pCoeffs);
}

inline void CLineModelElement::SetAllOffsetsEnabled(Float64 val) {
  for (size_t i = 0; i < GetSize(); ++i) {
    if (m_ElementParam->m_Lines[i].GetOffsetFitEnabled())
      m_ElementParam->m_Offsets[i] = val;
  }
}

inline void CLineModelElement::SetOffset(Int32 lineIdx, Float64 val) {
  m_ElementParam->m_Offsets[lineIdx] = val;
}

inline void CLineModelElement::SetLineProfile(Int32 lineIdx,
                                              CLineProfile_ptr &&profile) {
  m_ElementParam->m_Lines[lineIdx].SetProfile(std::move(profile));
}

/**
 * \brief  returns a call to the m_Lines GetName.
 **/
inline const std::string &CLineModelElement::GetLineName(Int32 subeIdx) const {
#ifdef DEBUG
  if (subeIdx >= GetSize())
    THROWG(INTERNAL_ERROR, "invalid index");
#endif

  return m_ElementParam->m_Lines[subeIdx].GetName();
}

/**
 * \brief Returns the content of the m_SignFactors with index equal to the
 *argument.
 **/
inline Float64 CLineModelElement::GetSignFactor(Int32 subeIdx) const {
  return m_SignFactors[subeIdx];
}

inline bool CLineModelElement::SetAbsLinesLimit(Float64 limit) {
  m_absLinesLimit = limit;
  return true;
}

inline Float64 CLineModelElement::GetAbsLinesLimit() const {
  return m_absLinesLimit;
}

} // namespace NSEpic

#endif // _REDSHIFT_LINEMODEL_ELEMENT_
