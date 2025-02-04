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

#include <algorithm>
#include <cfloat>
#include <climits>

#include "RedshiftLibrary/linemodel/element.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace std;
using namespace NSEpic;

TLineModelElementParam::TLineModelElementParam(CLineVector lines,
                                               Float64 velocity,
                                               const std::string &lineWidthType)
    : m_Lines(std::move(lines)), m_Velocity(velocity), m_VelocityStd(NAN),
      m_FittedAmplitudes(m_Lines.size(), NAN),
      m_FittedAmplitudesStd(m_Lines.size(), NAN),
      m_OffsetsStd(m_Lines.size(), NAN), m_fittingGroupInfo(undefStr),
      m_defaultVelocity(velocity),
      m_globalOutsideLambdaRangeList(m_Lines.size(), false) {
  m_NominalAmplitudes.reserve(m_Lines.size());
  m_Offsets.reserve(m_Lines.size());
  for (Int32 index = 0; index != m_Lines.size(); ++index) {
    auto const &line = m_Lines[index];
    m_NominalAmplitudes.push_back(line.GetNominalAmplitude());
    m_Offsets.push_back(line.GetOffset());
    m_LinesIds[line.GetID()] = index;
  }
  if (m_Lines.empty())
    THROWG(ErrorCode::INTERNAL_ERROR, "Empty line vector");
  auto const &first_line = m_Lines.front();

  m_type = first_line.GetType();
  m_isEmission = first_line.IsEmission();
  init(lineWidthType);
}

// for unit test

void TLineModelElementParam::init(const std::string &widthType) {
  if (widthType == "instrumentDriven") {
    m_LineWidthType = INSTRUMENTDRIVEN;
  } else if (widthType == "combined") {
    m_LineWidthType = COMBINED;
  } else if (widthType == "velocityDriven") {
    m_LineWidthType = VELOCITYDRIVEN;
  } else {
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Unknown LineWidthType" << widthType);
  }

  Float64 sign = IsEmission() ? 1.0 : -1.0;
  m_SignFactors.assign(size(), sign);
  for (Int32 index = 0; index != size(); ++index) {
    if (getLineProfile(index)->isAsym() || getLineProfile(index)->isSymIgm())
      m_asymLineIndices.push_back(index);
  }
}
/**
 * \brief Returns the nominal amplitude of line with id line_id.
 **/
Float64 TLineModelElementParam::GetNominalAmplitude(Int32 line_index) const {
  return m_NominalAmplitudes[line_index];
}

/**
 * \brief Set the nominal amplitude of the line with id line_id.
 **/
bool TLineModelElementParam::SetNominalAmplitude(Int32 line_index,
                                                 Float64 nominalamp) {
  m_NominalAmplitudes[line_index] = nominalamp;
  return true;
}

Int32 TLineModelElementParam::getSignFactor(Int32 line_index) const {
  return m_SignFactors[line_index];
}

/**
 * \brief Returns the fitted amplitude for the argument index.
 **/
Float64 TLineModelElementParam::GetFittedAmplitude(Int32 line_index) const {
  return m_FittedAmplitudes[line_index];
}

bool TLineModelElementParam::isAllAmplitudesNull() const {
  auto const &amps = m_FittedAmplitudes;
  return std::find_if(amps.cbegin(), amps.cend(),
                      [](Float64 a) { return a > 0.0; }) == amps.cend();
}

void TLineModelElementParam::setNullNominalAmplitudesNotFittable() {
  m_nullNominalAmplitudes = true;
  for (size_t idx = 0; idx < size(); ++idx) {
    if (!isOutsideLambdaRangeLine(idx) && GetNominalAmplitude(idx) > 0.0) {
      m_nullNominalAmplitudes = false;
      break;
    }
  }
}

/**
 * \brief Returns the fitted amplitude error for the argument index.
 **/
Float64 TLineModelElementParam::GetFittedAmplitudeStd(Int32 line_index) const {
  return m_FittedAmplitudesStd[line_index];
}

/*
const CLineProfile_ptr &
TLineModelElementParam::getLineProfile(Int32 line_index) const {
  return getLineProfile(line_index);
}

*/

void TLineModelElementParam::setVelocity(Float64 vel) {
#ifdef DEBUG
  if (!GetSize())
    THROWG(ErrorCode::INTERNAL_ERROR,
           "Empty line model element, could not set velocity");
#endif
  m_Velocity = isfinite(vel) ? vel : m_defaultVelocity;
}

void TLineModelElementParam::setVelocityStd(Float64 velStd) {
#ifdef DEBUG
  if (!GetSize())
    THROWG(ErrorCode::INTERNAL_ERROR,
           "Empty line model element, could not set velocity");
#endif
  m_VelocityStd = velStd;
}

/*
Float64 TLineModelElementParam::GetSumCross() const {
  return m_sumCross;
}

void TLineModelElementParam::SetSumCross(Float64 val) {
  m_sumCross = val;
}

Float64 TLineModelElementParam::GetDtmFree() const {
  return m_dtmFree;
}

void TLineModelElementParam::SetDtmFree(Float64 val) {
  m_dtmFree = val;
}

Float64 TLineModelElementParam::GetSumGauss() const {
  return m_sumGauss;
}

void TLineModelElementParam::SetSumGauss(Float64 val) {
  m_sumGauss = val;
}

const std::string &TLineModelElementParam::GetFittingGroupInfo() const {
  return m_fittingGroupInfo;
}

const TPolynomCoeffs &TLineModelElementParam::GetPolynomCoeffs() const {
  return m_ampOffsetsCoeffs;
}

void TLineModelElementParam::SetLineProfile(Int32 line_index,
                                              CLineProfile_ptr &&profile) {
  SetLineProfile(line_index, std::move(profile));
}
*/
/**
 * \brief  returns a call to the m_Lines GetName.
 **/
const std::string &TLineModelElementParam::GetLineName(Int32 line_index) const {
  return m_Lines[line_index].GetName();
}

/**
 * \brief Returns the content of the m_SignFactors with index equal to the
 *argument.
 **/
Float64 TLineModelElementParam::GetSignFactor(Int32 line_index) const {
  return m_SignFactors[line_index];
}

bool TLineModelElementParam::SetAbsLinesLimit(Float64 limit) {
  m_absLinesLimit = limit;
  return true;
}

Float64 TLineModelElementParam::GetAbsLinesLimit() const {
  return m_absLinesLimit;
}

/**
 * \brief If the fitted amplitude of line with id line_id is above the limit,
 *sets it to either that limit or zero, whichever is greater.
 **/
bool TLineModelElementParam::LimitFittedAmplitude(Int32 line_index,
                                                  Float64 limit) {

  bool limited = false;
  if (GetFittedAmplitude(line_index) > limit) {
    limited = true;
    m_FittedAmplitudes[line_index] = std::max(0.0, limit);

    // now update the amplitude of the other lines
    Float64 amplitudeRef =
        GetFittedAmplitude(line_index) / m_NominalAmplitudes[line_index];
    for (Int32 index = 0; index != size(); ++index) {
      m_FittedAmplitudes[index] = m_NominalAmplitudes[index] * amplitudeRef;
    }
  }
  return limited;
}

Float64 TLineModelElementParam::GetMaxNominalAmplitude() const {
  auto it = std::max_element(m_NominalAmplitudes.cbegin(),
                             m_NominalAmplitudes.cend());
  return (*it);
}

void TLineModelElementParam::SetAllOffsetsEnabled(Float64 val) {
  for (Int32 index = 0; index != size(); ++index) {
    if (GetLines()[index].IsOffsetFitEnabled())
      m_Offsets[index] = val;
  }
}

Float64 TLineModelElementParam::GetLineProfileDerivVel(
    const CLineProfile &profile, Float64 x, Float64 x0, Float64 sigma,
    bool isEmission) const {
  const Float64 c = SPEED_OF_LIGHT_IN_VACCUM;
  const Float64 pfsSimuCompensationFactor = 1.0; // TODO should be removed !!!
  Float64 v = getVelocity(), v_to_sigma = pfsSimuCompensationFactor / c * x0;

  Float64 profile_derivSigma = profile.GetLineProfileDerivSigma(x, x0, sigma);

  // sincs lsf is an instrumental response, then derivative of this latter with
  // respect to velocity is null
  switch (m_LineWidthType) {
  case INSTRUMENTDRIVEN:
    return 0.0;
  case COMBINED:
    return v_to_sigma * v_to_sigma * v / sigma * profile_derivSigma;
  case VELOCITYDRIVEN:
    return v_to_sigma * profile_derivSigma;
  default:
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Invalid LineWidthType : " << m_LineWidthType);
  }
  return 0.0;
}

Float64 TLineModelElementParam::GetLineProfileDerivVel(Int32 index, Float64 x,
                                                       Float64 x0,
                                                       Float64 sigma) const {
  return GetLineProfileDerivVel(*getLineProfile(index), x, x0, sigma,
                                m_isEmission);
}
