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

#ifndef _REDSHIFT_LINEMODEL_ELEMENT_PARAM_
#define _REDSHIFT_LINEMODEL_ELEMENT_PARAM_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/line/line.h"
#include "RedshiftLibrary/line/lineprofile.h"

namespace NSEpic {

// TODO should be defined elsewhere
enum TLineWidthType { INSTRUMENTDRIVEN, COMBINED, VELOCITYDRIVEN };

enum class ElementComposition {
  Default, // follow input linecatalog using AmplitudeGroupName: one amplitude
           // by element
  EmissionAbsorption, // for tpl_ratio: all emission lines in one element, all
                      // abs lines in a 2nd element
  OneLine // for linemeas, one unique line by element, ie. don't follow input
          // linecatalog AmplitudeGroupName
};

struct TLineModelElementParam {

  TLineModelElementParam(CLineVector lines, Float64 velocity,
                         const std::string &lineWidthType);

  CLineVector m_Lines;
  Float64 m_Velocity = NAN;
  Float64 m_VelocityStd = NAN;
  TFloat64List m_FittedAmplitudes;
  TFloat64List m_FittedAmplitudesStd;
  TFloat64List m_NominalAmplitudes;
  TFloat64List m_Offsets;
  TFloat64List m_OffsetsStd;
  TInt32Map m_LinesIds;
  std::string m_fittingGroupInfo;
  TPolynomCoeffs m_ampOffsetsCoeffs;
  Float64 m_sumCross = 0.0;
  Float64 m_sumGauss = 0.0;
  Float64 m_dtmFree =
      0.0; // dtmFree is the non-positive-constrained version of sumCross

  TLineWidthType m_LineWidthType;
  TFloat64List m_SignFactors;

  TInt32List m_asymLineIndices;
  Float64 m_absLinesLimit =
      1.0; //-1: disable the ABS lines amplitude cut, any other value is
           // used as a limit for the abs line coeff (typically: 1.0)
  Float64 m_defaultVelocity;
  CLine::EType m_type;
  bool m_isEmission;

  // global (considering all spectra) validity of the element:
  // element is outside range (all element lines, ie no line is visible)
  bool m_globalOutsideLambdaRange = false;
  // visibility of each line of the element
  std::vector<bool> m_globalOutsideLambdaRangeList;
  // null continuum under all element absorption lines
  bool m_absLinesNullContinuum = false;
  // all visible element lines have null nominal amplitudes
  bool m_nullNominalAmplitudes = false;
  // the profile is null on the element (all lines) support
  bool m_nullLineProfiles = false;

  void init(const std::string &widthType);
  const Float64 &getSumGauss() const { return m_sumGauss; }
  const Float64 &getSumCross() const { return m_sumCross; }

  void SetSumCross(Float64 val) { m_sumCross = val; }
  Float64 GetDtmFree() const { return m_dtmFree; }
  void SetDtmFree(Float64 val) { m_dtmFree = val; }
  void SetSumGauss(Float64 val) { m_sumGauss = val; }

  const std::string &getFittingGroupInfo() const { return m_fittingGroupInfo; }
  // const Float64 &getSumGauss() const {return m_sumGauss;}

  // TODO this is ugly, and maybe m_SignFactor should not exist, knowing
  // m_type/m_isEmission should be enough
  Int32 getSignFactor(Int32 line_index) const;
  Float64 GetSignFactor(Int32 line_index) const;

  Float64 getVelocity() const { return m_Velocity; }
  Float64 getVelocityStd() const { return m_VelocityStd; };

  TAsymParams GetAsymfitParams(Int32 asym_line_index = 0) const {
    if (!m_asymLineIndices.size())
      return TAsymParams(); // case where no asymprofile in linecatalog
    return m_Lines[m_asymLineIndices[asym_line_index]].GetAsymParams();
  }
  TSymIgmParams GetSymIgmParams(Int32 asym_line_index = 0) const {
    if (!m_asymLineIndices.size())
      return TSymIgmParams(); // case where no asymprofile in linecatalog
    return m_Lines[m_asymLineIndices[asym_line_index]].GetSymIgmParams();
  }

  void SetAsymfitParams(const TAsymParams &params,
                        Int32 line_index = undefIdx) {
    if (!m_asymLineIndices.size())
      return;
    if (line_index >= 0)
      m_Lines[line_index].SetAsymParams(params);
    else
      for (auto index : m_asymLineIndices)
        m_Lines[index].SetAsymParams(params);
  }

  void SetSymIgmParams(const TSymIgmParams &params,
                       Int32 line_index = undefIdx) {
    if (!m_asymLineIndices.size())
      return;
    if (line_index >= 0)
      m_Lines[line_index].SetSymIgmParams(params);
    else
      for (auto index : m_asymLineIndices)
        m_Lines[index].SetSymIgmParams(params);
  }

  void SetSymIgmFit(bool val = true, Int32 line_index = undefIdx) {
    if (!m_asymLineIndices.size())
      return;
    if (line_index >= 0)
      m_Lines[line_index].SetSymIgmFit(val);
    else
      for (auto index : m_asymLineIndices)
        m_Lines[index].SetSymIgmFit(val);
  }

  void setFittedAmplitude(Int32 index, Float64 fittedAmp, Float64 fittedAmpStd,
                          Float64 nominalAmplitude = 1.0) {
    auto &fa = m_FittedAmplitudes;
    auto &fastd = m_FittedAmplitudesStd;
    fa[index] = fittedAmp * nominalAmplitude;
    // limit the absorption to 0.0-1.0, so that it's never <0
    if (m_SignFactors[index] == -1 && m_absLinesLimit > 0.0 &&
        fa[index] > m_absLinesLimit) {
      fa[index] = m_absLinesLimit;
    }
    fastd[index] = fittedAmpStd * nominalAmplitude;
  }

  void setAmplitudes(Float64 A, Float64 AStd) {
    auto &fa = m_FittedAmplitudes;
    auto &fastd = m_FittedAmplitudesStd;
    auto &na = m_NominalAmplitudes;

    if (std::isnan(A) || isNotFittable()) {
      fa.assign(size(), NAN);
      fastd.assign(size(), NAN);
      return;
    }

    for (Int32 index = 0; index != size(); ++index) {
      if (m_globalOutsideLambdaRangeList[index]) {
        fa[index] = NAN;
        fastd[index] = NAN;
        continue;
      }
      setFittedAmplitude(index, A, AStd, na[index]);
    }
  }

  Int32 size() const { return m_Lines.size(); }

  Int32 getLineIndex(Int32 line_id) const;
  Int32 getLineIndex(const std::string &LineTagStr) const;
  bool hasLine(Int32 line_id) const {
    return getLineIndex(line_id) != undefIdx;
  };
  bool hasLine(const std::string &LineTagStr) const {
    return getLineIndex(LineTagStr) != undefIdx;
  };

  const CLineProfile_ptr &getLineProfile(Int32 line_index) const {
    return m_Lines[line_index].GetProfile();
  }

  void SetLineProfile(Int32 line_index, CLineProfile_ptr &&profile) {
    m_Lines[line_index].SetProfile(std::move(profile));
  }

  Int32 GetElementValidLine() const {
    if (isNotFittable())
      return undefIdx;
    for (Int32 line_idx = 0; line_idx != size(); ++line_idx) {
      auto const &nominal_amplitude = m_NominalAmplitudes[line_idx];
      if (!std::isnan(m_FittedAmplitudes[line_idx]) && nominal_amplitude != 0.0)
        return line_idx;
    }
    return undefIdx;
  }

  Float64 GetElementAmplitude() const {
    Int32 line_idx = GetElementValidLine();
    if (line_idx == undefIdx)
      return NAN;
    return m_FittedAmplitudes[line_idx] / m_NominalAmplitudes[line_idx];
  }

  Float64 GetElementAmplitudeError() const {
    Int32 line_idx = GetElementValidLine();
    if (line_idx == undefIdx)
      return NAN;
    return m_FittedAmplitudesStd[line_idx] / m_NominalAmplitudes[line_idx];
  }

  void resetLambdaOffsets() {
    for (Int32 line_idx = 0; line_idx != size(); ++line_idx)
      setLambdaOffset(line_idx, m_Lines[line_idx].GetOffset());
  }

  void setLambdaOffset(Int32 line_index, Float64 val) {
    m_Offsets[line_index] = val;
  }
  void setLambdaOffsetStd(Int32 line_index, Float64 val) {
    m_OffsetsStd[line_index] = val;
  }
  Float64 getLambdaOffset(Int32 line_index) { return m_Offsets[line_index]; }
  Float64 getLambdaOffsetStd(Int32 line_index) {
    return m_OffsetsStd[line_index];
  }

  void resetFittingParams();
  void resetAsymfitParams() {
    for (auto index : m_asymLineIndices)
      m_Lines[index].resetAsymFitParams();
  }
  void SetFittingGroupInfo(const std::string &val) { m_fittingGroupInfo = val; }
  void resetAmplitudeOffset() { m_ampOffsetsCoeffs = TPolynomCoeffs(); }
  void SetPolynomCoeffs(const TPolynomCoeffs &coeffs) {
    m_ampOffsetsCoeffs = coeffs;
  }
  /**
   * @brief Look for polynom coeffs corresponding to one specific Line
   *
   * @param eIdx
   * @return TPolynomCoeffs
   */
  const TPolynomCoeffs &GetPolynomCoeffs() const { return m_ampOffsetsCoeffs; }

  CLine::EType GetElementType() const { return m_type; };
  bool IsEmission() const { return m_isEmission; };
  bool IsAbsorption() const { return !m_isEmission; };

  bool isFittable() const { return !isNotFittable(); }

  bool isNotFittable() const {
    return m_globalOutsideLambdaRange || m_nullNominalAmplitudes ||
           m_absLinesNullContinuum || m_nullLineProfiles;
  }

  bool isOutsideLambdaRangeLine(Int32 line_index) const {
    return m_globalOutsideLambdaRangeList[line_index];
  }

  const CLineVector &GetLines() const { return m_Lines; };
  const TInt32Map &GetLinesIndices() const { return m_LinesIds; };

  const std::string &GetLineName(Int32 line_index) const;

  Float64 GetNominalAmplitude(Int32 line_index) const;
  bool SetNominalAmplitude(Int32 line_index, Float64 nominalamp);
  Float64 GetMaxNominalAmplitude() const;
  Float64 GetFittedAmplitude(Int32 line_index) const;
  Float64 GetFittedAmplitudeStd(Int32 line_index) const;

  bool isAllAmplitudesNull() const;
  void setNullNominalAmplitudesNotFittable();

  bool LimitFittedAmplitude(Int32 line_index, Float64 limit);
  void SetAllOffsetsEnabled(Float64 val);
  bool SetAbsLinesLimit(Float64 limit);
  Float64 GetAbsLinesLimit() const;

  void setVelocity(Float64 vel);
  void setVelocityStd(Float64 velStd);

  Float64 GetLineProfileDerivVel(const CLineProfile &profile, Float64 x,
                                 Float64 x0, Float64 sigma,
                                 bool isEmission) const;

  Float64 GetLineProfileDerivVel(Int32 index, Float64 x, Float64 x0,
                                 Float64 sigma) const;
};

} // namespace NSEpic
#endif
