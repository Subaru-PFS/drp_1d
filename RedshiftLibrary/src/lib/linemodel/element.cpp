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

#include "RedshiftLibrary/linemodel/element.h"

#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

#include <algorithm>
#include <cfloat>
#include <climits>

using namespace std;
using namespace NSEpic;

TLineModelElementParam::TLineModelElementParam(CLineVector lines,
                                               Float64 velocityEmission,
                                               Float64 velocityAbsorption)
    : m_Lines(std::move(lines)), m_VelocityEmission(velocityEmission),
      m_VelocityAbsorption(velocityAbsorption),
      m_FittedAmplitudes(m_Lines.size(), NAN),
      m_FittedAmplitudeErrorSigmas(m_Lines.size(), NAN),
      m_fittingGroupInfo(undefStr) {
  m_NominalAmplitudes.reserve(m_Lines.size());
  m_Offsets.reserve(m_Lines.size());
  for (Int32 index = 0; index != m_Lines.size(); ++index) {
    auto const &line = m_Lines[index];
    m_NominalAmplitudes.push_back(line.GetNominalAmplitude());
    m_Offsets.push_back(line.GetOffset());
    m_LinesIds[line.GetID()] = index;
  }
  auto const &first_line = m_Lines.front();

  m_type = first_line.GetType();
  m_isEmission = first_line.IsEmission();
}

/**
 * \brief Constructs the object setting members according to arguments and
 *defaults.
 **/
CLineModelElement::CLineModelElement(
    const TLineModelElementParam_ptr elementParam, const std::string &widthType)
    : m_ElementParam(std::move(elementParam)),
      m_OutsideLambdaRangeOverlapThreshold(
          0.33), // 33% overlap minimum in order to keep the line
      m_OutsideLambdaRange(
          true) // example: 0.33 means 66% of the line is allowed to be outside
                // the spectrum with the line still considered inside the
                // lambda range

{

  m_size = m_ElementParam->m_Lines.size();

  if (widthType == "instrumentdriven") {
    m_ElementParam->m_LineWidthType = INSTRUMENTDRIVEN;
  } else if (widthType == "combined") {
    m_ElementParam->m_LineWidthType = COMBINED;
  } else if (widthType == "velocitydriven") {
    m_ElementParam->m_LineWidthType = VELOCITYDRIVEN;
  } else {
    THROWG(INTERNAL_ERROR, Formatter() << "Unknown LineWidthType" << widthType);
  }

  Float64 sign = IsEmission() ? 1.0 : -1.0;
  m_ElementParam->m_SignFactors.assign(GetSize(), sign);
  for (Int32 index = 0; index != GetSize(); ++index) {
    if (getLineProfile(index)->isAsym() || getLineProfile(index)->isSymIgm())
      m_ElementParam->m_asymLineIndices.push_back(index);
  }

  SetElementAmplitude(NAN, NAN);
}

void CLineModelElement::reset() {
  // init the fitted amplitude values and related variables
  m_ElementParam->m_FittedAmplitudes.assign(GetSize(), NAN);
  m_ElementParam->m_FittedAmplitudeErrorSigmas.assign(GetSize(), NAN);

  SetFittingGroupInfo(undefStr);
  m_ElementParam->m_sumGauss = NAN;
  m_ElementParam->m_sumCross = NAN;
  m_ElementParam->m_dtmFree = NAN;
}

Int32 CLineModelElement::getLineIndex(Int32 line_id) const {
  Int32 index = undefIdx;
  auto const &indicesMap = GetLinesIndices();
  auto const &it = indicesMap.find(line_id);
  if (it != indicesMap.cend())
    index = it->second;
  return index;
}

Int32 CLineModelElement::getLineIndex(const std::string &LineTagStr) const {
  Int32 line_index = undefIdx;
  auto const &lines = GetLines();
  auto it =
      std::find_if(lines.cbegin(), lines.cend(), [&LineTagStr](auto &line) {
        auto const &name = line.GetName();
        std::size_t foundstra = name.find(LineTagStr.c_str());
        return foundstra != std::string::npos;
      });
  if (it != lines.cend())
    line_index = it - lines.cbegin();

  return line_index;
}

// redirecting to the lsf method for computing instrument responce
/**
 * Get instrumental response (including source response) from LSF
 * combine quadratically with the instrinsic width of the Line itself. The line
 * width in this case represents to the velocity
 * */
Float64 CLineModelElement::GetLineWidth(Float64 redshiftedlambda,
                                        bool isEmission) const {
  const Float64 c = m_speedOfLightInVacuum;
  Float64 v = isEmission ? getVelocityEmission() : getVelocityAbsorption();
  const Float64 pfsSimuCompensationFactor = 1.0;

  if (!m_LSF)
    THROWG(INTERNAL_ERROR, "LSF object is not initailized.");
  Float64 instrumentSigma = m_LSF->GetWidth(redshiftedlambda);

  Float64 velocitySigma = pfsSimuCompensationFactor * v / c *
                          redshiftedlambda; //, useless /(1+z)*(1+z);
  switch (m_ElementParam->m_LineWidthType) {
  case INSTRUMENTDRIVEN: // only instrumental sigma
    velocitySigma = 0.;
    break;
  case VELOCITYDRIVEN: // only velocity sigma
    instrumentSigma = 0.;
    break;
  case COMBINED: // combination of the two
    break;
  default:
    // TODO this should not happen here, but at parameter setting stage
    THROWG(INTERNAL_ERROR, Formatter() << "Invalid LSF type "
                                       << m_ElementParam->m_LineWidthType);
  }

  Float64 sigma =
      sqrt(instrumentSigma * instrumentSigma + velocitySigma * velocitySigma);
  return sigma;
}

Float64 CLineModelElement::GetLineProfileDerivVel(const CLineProfile &profile,
                                                  Float64 x, Float64 x0,
                                                  Float64 sigma,
                                                  bool isEmission) const {
  const Float64 c = m_speedOfLightInVacuum;
  const Float64 pfsSimuCompensationFactor = 1.0;
  Float64 v = isEmission ? getVelocityEmission() : getVelocityAbsorption(),
          v_to_sigma = pfsSimuCompensationFactor / c * x0;

  Float64 profile_derivSigma = profile.GetLineProfileDerivSigma(x, x0, sigma);

  // sincs lsf is an instrumental response, then derivative of this latter with
  // respect to velocity is null
  switch (m_ElementParam->m_LineWidthType) {
  case INSTRUMENTDRIVEN:
    return 0.0;
  case COMBINED:
    return v_to_sigma * v_to_sigma * v / sigma * profile_derivSigma;
  case VELOCITYDRIVEN:
    return v_to_sigma * profile_derivSigma;
  default:
    THROWG(INTERNAL_ERROR, Formatter() << "Invalid LineWidthType : "
                                       << m_ElementParam->m_LineWidthType);
  }
  return 0.0;
}

Float64 CLineModelElement::getVelocity() const {
  if (!GetSize())
    return NAN;

  return IsEmission() ? getVelocityEmission() : getVelocityAbsorption();
}

/**
 * @brief GetContinuumAtCenterProfile
 * @param line_id
 * @param spectralAxis
 * @param redshift
 * @param lambdaRange
 * @param continuumfluxAxis
 * Add up polynome contribution below lines
 * @return the continuum flux val at the sub element center wavelength. Error
 * returns -999/-9999 if center profile not in range
 *
 */
Float64 CLineModelElement::GetContinuumAtCenterProfile(
    Int32 line_id, const CSpectrumSpectralAxis &spectralAxis, Float64 redshift,
    const CSpectrumFluxAxis &continuumfluxAxis,
    bool enableAmplitudeOffsets) const {
  Float64 mu = GetObservedPosition(line_id, redshift);

  Int32 IdxCenterProfile = spectralAxis.GetIndexAtWaveLength(mu);
  if (IdxCenterProfile < 0 ||
      IdxCenterProfile > continuumfluxAxis.GetSamplesCount() - 1) {
    return NAN;
  }

  Float64 cont = continuumfluxAxis[IdxCenterProfile];

  if (enableAmplitudeOffsets)
    cont += GetPolynomCoeffs().getValue(spectralAxis[IdxCenterProfile]);

  return cont;
}

/**
 * \brief Returns the theoretical support range for the line
 **/
void CLineModelElement::EstimateTheoreticalSupport(
    Int32 line_index, const CSpectrumSpectralAxis &spectralAxis,
    Float64 redshift, const TFloat64Range &lambdaRange) {
  Float64 mu = GetObservedPosition(line_index, redshift);
  if (!m_LSF->checkAvailability(mu)) {
    m_OutsideLambdaRangeList[line_index] = true;
    return;
  }
  Float64 sigma = GetLineWidth(mu, IsEmission());
  Float64 winsize = getLineProfile(line_index)->GetNSigmaSupport() * sigma;
  TInt32Range supportRange =
      EstimateIndexRange(spectralAxis, mu, lambdaRange, winsize);

  m_StartTheoretical[line_index] = supportRange.GetBegin();
  m_EndTheoretical[line_index] = supportRange.GetEnd();
  m_StartNoOverlap[line_index] = supportRange.GetBegin();
  m_EndNoOverlap[line_index] = supportRange.GetEnd();

  if (supportRange.GetBegin() >
      supportRange.GetEnd()) // in this case the line is completely outside the
                             // lambdarange
  {
    m_OutsideLambdaRangeList[line_index] = true;
  } else { // in this case the line is completely inside the lambdarange or with
           // partial overlap

    Float64 minLineOverlap = m_OutsideLambdaRangeOverlapThreshold * winsize;
    Float64 startLbda = spectralAxis[m_StartNoOverlap[line_index]];
    Float64 endLbda = spectralAxis[m_EndNoOverlap[line_index]];

    if (startLbda >= (lambdaRange.GetEnd() - minLineOverlap) ||
        endLbda <= (lambdaRange.GetBegin() + minLineOverlap)) {
      m_OutsideLambdaRangeList[line_index] = true;
    } else {
      m_OutsideLambdaRangeList[line_index] = false;
    }
  }

  return;
}

/**
 * \brief Returns the index range for a given window size (Angstrom)
 **/
TInt32Range CLineModelElement::EstimateIndexRange(
    const CSpectrumSpectralAxis &spectralAxis, Float64 mu,
    const TFloat64Range &lambdaRange, Float64 winsizeAngstrom) {
  TInt32Range supportRange;
  Float64 winsize = winsizeAngstrom;

  Float64 lambda_start = mu - winsize / 2.0;
  if (lambda_start < lambdaRange.GetBegin()) {
    lambda_start = lambdaRange.GetBegin();
  }
  supportRange.SetBegin(spectralAxis.GetIndexAtWaveLength(lambda_start));

  Float64 lambda_end = mu + winsize / 2.0;
  if (lambda_end > lambdaRange.GetEnd()) {
    lambda_end = lambdaRange.GetEnd();
  }
  supportRange.SetEnd(spectralAxis.GetIndexAtWaveLength(lambda_end));

  // correct the end value if higher then lambdaRange end
  // Log.LogDebug( "    multiline: spectralAxis[supportRange.GetEnd()] = %f",
  // spectralAxis[supportRange.GetEnd()]);
  if (spectralAxis[supportRange.GetEnd()] > lambdaRange.GetEnd()) {
    supportRange.SetEnd(supportRange.GetEnd() - 1);
  }
  // correct the end value if not higher or equal to the begin value
  if (supportRange.GetEnd() < supportRange.GetBegin()) {
    supportRange.SetEnd(supportRange.GetBegin() - 1);
  }

  return supportRange;
}

// set the global outside lambda range
void CLineModelElement::SetOutsideLambdaRange() {
  m_OutsideLambdaRange = true;
  for (auto const &outside_lambda_range_id : m_OutsideLambdaRangeList)
    m_OutsideLambdaRange = m_OutsideLambdaRange && outside_lambda_range_id;
}

/**
 * \brief Limits each m_Lines element within the argument lambdaRange, and sets
 *the m_FittedAmplitudes to -1. Sets the global outside lambda range. Inits the
 *fitted amplitude values.
 **/
void CLineModelElement::prepareSupport(
    const CSpectrumSpectralAxis &spectralAxis, Float64 redshift,
    const TFloat64Range &lambdaRange) {

  Int32 nLines = GetSize();
  m_OutsideLambdaRange = true;
  m_StartNoOverlap.assign(nLines, undefIdx);
  m_EndNoOverlap.assign(nLines, undefIdx);
  m_StartTheoretical.assign(nLines, undefIdx);
  m_EndTheoretical.assign(nLines, undefIdx);
  m_OutsideLambdaRangeList.assign(nLines, true);
  m_LineIsActiveOnSupport.assign(nLines, TBoolList(nLines, false));

  for (Int32 index = 0; index != nLines; ++index) {
    EstimateTheoreticalSupport(index, spectralAxis, redshift, lambdaRange);
    // set the lines active on their own support
    m_LineIsActiveOnSupport[index][index] = true;
  }
  SetOutsideLambdaRange();

  bool supportNoOverlap_has_duplicates = true;
  Int32 x1 = 0;
  Int32 y1 = 0;
  Int32 x2 = 0;
  Int32 y2 = 0;
  Int32 icmpt = 0;
  Int32 ncmpt = 20;
  while (supportNoOverlap_has_duplicates && icmpt < ncmpt) {
    icmpt++;
    for (Int32 i = 0; i != nLines; ++i) {
      if (m_OutsideLambdaRangeList[i])
        continue;

      if (m_StartNoOverlap[i] > m_EndNoOverlap[i])
        continue;

      for (Int32 j = 0; j != nLines; ++j) {
        if (m_OutsideLambdaRangeList[j])
          continue;

        if (m_StartNoOverlap[j] > m_EndNoOverlap[j])
          continue;

        if (i == j)
          continue;

        bool lineActiveSupportToBeCorrected = false;
        //
        x1 = m_StartNoOverlap[i];
        x2 = m_EndNoOverlap[i];
        y1 = m_StartNoOverlap[j];
        y2 = m_EndNoOverlap[j];
        // compute overlapping region
        Int32 max = std::max(x1, y1);
        Int32 min = std::min(x2, y2);
        if (max - min < 0) { // case of overlapping
          m_StartNoOverlap[i] = std::min(x1, y1);
          m_EndNoOverlap[i] = std::max(x2, y2);
          m_StartNoOverlap[j] =
              m_EndNoOverlap[i]; // deactivate j when end is start -1
          m_EndNoOverlap[j] = m_EndNoOverlap[i] - 1; // deactivate j

          lineActiveSupportToBeCorrected = true;
        }

        if (lineActiveSupportToBeCorrected) {
          // set the lines active on the overlapping support
          m_LineIsActiveOnSupport[i][j] = true;
          m_LineIsActiveOnSupport[j][i] = true;
          // append all the previously overlapping lines as active on the
          // support
          for (Int32 i2 = 0; i2 != nLines; ++i2) {
            if (m_OutsideLambdaRangeList[i2])
              continue;

            if (m_LineIsActiveOnSupport[i][i2]) {
              m_LineIsActiveOnSupport[i2][j] = true;
              m_LineIsActiveOnSupport[j][i2] = true;
            }
          }
          // append all the previously overlapping lines as active on the
          // support
          for (Int32 j2 = 0; j2 != nLines; ++j2) {
            if (m_OutsideLambdaRangeList[j2])
              continue;

            if (m_LineIsActiveOnSupport[j][j2]) {
              m_LineIsActiveOnSupport[j2][i] = true;
              m_LineIsActiveOnSupport[i][j2] = true;
            }
          }
        }
      }
    }

    supportNoOverlap_has_duplicates = false;

    // check that there are no overlapping sub-supports in the list
    for (Int32 i = 0; i != nLines; ++i) {
      if (supportNoOverlap_has_duplicates) {
        break;
      }
      if (m_OutsideLambdaRangeList[i]) {
        continue;
      }
      for (Int32 j = 0; j != nLines; ++j) {
        if (m_OutsideLambdaRangeList[j]) {
          continue;
        }
        if (i == j) {
          continue;
        }

        x1 = m_StartNoOverlap[i];
        x2 = m_EndNoOverlap[i];
        y1 = m_StartNoOverlap[j];
        y2 = m_EndNoOverlap[j];
        Int32 max = std::max(x1, y1);
        Int32 min = std::min(x2, y2);
        if (max - min < 0) {
          supportNoOverlap_has_duplicates = true;
          break;
        }
      }
    }
  }
}

/**
 * \brief Creates an empty list of ranges as the return value. If not
 *m_OutsideLambdaRange, for each m_Lines element which is also not outside
 *lambda range, add its support to the return value.
 **/
TInt32RangeList CLineModelElement::getSupport() const {
  TInt32RangeList support;
  if (m_OutsideLambdaRange == false) {
    for (Int32 index = 0; index != GetSize(); ++index) {
      if (m_OutsideLambdaRangeList[index])
        continue;

      support.push_back(
          TInt32Range(m_StartNoOverlap[index], m_EndNoOverlap[index]));
    }
  }
  return support;
}

TInt32RangeList CLineModelElement::getTheoreticalSupport() const {
  TInt32RangeList support;

  if (m_OutsideLambdaRange == false) {
    for (Int32 index = 0; index != GetSize(); ++index) {
      if (m_OutsideLambdaRangeList[index])
        continue;

      support.push_back(
          TInt32Range(m_StartTheoretical[index], m_EndTheoretical[index]));
    }
  }
  return support;
}

/**
 * \brief Creates an empty list of ranges as the return value. If not
 *m_OutsideLambdaRange, for each m_Lines element belonging to the argument
 *subeIdx which is also not outside lambda range, add its support to the return
 *value.
 **/
TInt32Range CLineModelElement::getSupportSubElt(Int32 index) const {
  TInt32Range support =
      TInt32Range(m_StartNoOverlap[index], m_EndNoOverlap[index]);
  return support;
}

/**
 * \brief Returns the theoretical support of the line (sub-element).
 **/
TInt32Range CLineModelElement::getTheoreticalSupportSubElt(Int32 index) const {
  TInt32Range support =
      TInt32Range(m_StartTheoretical[index], m_EndTheoretical[index]);
  return support;
}

/**
 * \brief Calls GetLineWidth using the arguments and a calculated argument mu.
 **/
void CLineModelElement::getObservedPositionAndLineWidth(
    Int32 index, Float64 redshift, Float64 &mu, Float64 &sigma,
    bool doAsymfitdelta) const {
  mu = GetObservedPosition(index, redshift, doAsymfitdelta);
  if (!m_LSF->checkAvailability(mu)) {
    THROWG(INTERNAL_ERROR, "Line position does not belong to LSF range");
  } else
    sigma = GetLineWidth(mu, IsEmission());
  return;
}

/**
 * \brief Get the observed position of the sub-element with id line_id for a
 *given redshift
 **/
Float64 CLineModelElement::GetObservedPosition(Int32 index, Float64 redshift,
                                               bool doAsymfitdelta) const {
  Float64 dzOffset = m_ElementParam->m_Offsets[index] / m_speedOfLightInVacuum;

  auto const &line = GetLines()[index];
  Float64 mu = line.GetPosition() * (1 + redshift) * (1 + dzOffset);

  // deals with delta of asym profile
  if (doAsymfitdelta) {
    mu -= line.GetProfile()->GetDelta();
  }
  return mu;
}

/**
 * \brief Returns the line profile of the sub-element with id line_id at
 *wavelength x, for a given redshift.
 **/
Float64 CLineModelElement::GetLineProfileAtRedshift(Int32 index,
                                                    Float64 redshift,
                                                    Float64 x) const {
  Float64 mu = NAN;
  Float64 sigma = NAN;
  getObservedPositionAndLineWidth(index, redshift, mu, sigma,
                                  false); // do not apply Lya asym offset

  auto const &profile = getLineProfile(index);

  return profile->GetLineProfileVal(x, mu, sigma);
}

/**
 * \brief Returns NA  N if m_OutsideLambdaRange, and the fitted amplitude /
 *nominal amplitude of the first element otherwise.
 **/
Float64 CLineModelElement::GetElementAmplitude() const {
  if (m_OutsideLambdaRange) {
    return NAN;
  }
  for (Int32 index = 0; index != GetSize(); ++index) {
    auto const &nominal_amplitude = m_ElementParam->m_NominalAmplitudes[index];
    if (!m_OutsideLambdaRangeList[index] && nominal_amplitude != 0.0) {
      return m_ElementParam->m_FittedAmplitudes[index] / nominal_amplitude;
    }
  }
  return NAN;
}

/**
 * \brief Returns NAN if m_OutsideLambdaRange, and the fitted error / nominal
 *amplitude of the first element otherwise.
 **/
Float64 CLineModelElement::GetElementError() const {
  if (m_OutsideLambdaRange) {
    return NAN;
  }

  for (Int32 index = 0; index != GetSize(); ++index) {
    auto const &nominal_amplitude = m_ElementParam->m_NominalAmplitudes[index];
    if (!m_OutsideLambdaRangeList[index] && nominal_amplitude != 0.0) {
      return m_ElementParam->m_FittedAmplitudeErrorSigmas[index] /
             nominal_amplitude;
    }
  }
  return NAN;
}

Float64 CLineModelElement::GetMaxNominalAmplitude() const {
  auto it = std::max_element(m_ElementParam->m_NominalAmplitudes.cbegin(),
                             m_ElementParam->m_NominalAmplitudes.cend());
  return (*it);
}

void CLineModelElement::SetFittedAmplitude(Int32 index, Float64 fittedAmp,
                                           Float64 SNR) {

  if (m_OutsideLambdaRangeList[index] || std::isnan(fittedAmp)) {
    m_ElementParam->m_FittedAmplitudes[index] = NAN;
    m_ElementParam->m_FittedAmplitudeErrorSigmas[index] = NAN;
    return;
  }

  m_ElementParam->m_FittedAmplitudes[index] = fittedAmp;

  // limit the absorption to 0.0-1.0, so that it's never <0
  //*
  if (m_ElementParam->m_SignFactors[index] == -1 &&
      m_ElementParam->m_absLinesLimit > 0.0 &&
      m_ElementParam->m_FittedAmplitudes[index] >
          m_ElementParam->m_absLinesLimit) {
    m_ElementParam->m_FittedAmplitudes[index] = m_ElementParam->m_absLinesLimit;
  }

  m_ElementParam->m_FittedAmplitudeErrorSigmas[index] = SNR;
}

/**
 * \brief If outside lambda range, sets fitted amplitudes and errors to -1. If
 *inside, sets each line's fitted amplitude and error to -1 if line outside
 *lambda range, or amplitude to A * nominal amplitude and error to SNR * nominal
 *amplitude.
 **/
void CLineModelElement::SetElementAmplitude(Float64 A, Float64 SNR) {
  auto &fa = m_ElementParam->m_FittedAmplitudes;
  auto &faes = m_ElementParam->m_FittedAmplitudeErrorSigmas;
  auto &na = m_ElementParam->m_NominalAmplitudes;

  if (std::isnan(A) || m_OutsideLambdaRange) {
    fa.assign(GetSize(), NAN);
    faes.assign(GetSize(), NAN);
    return;
  }

  for (Int32 index = 0; index != GetSize(); ++index) {
    if (m_OutsideLambdaRangeList[index]) {
      fa[index] = NAN;
      faes[index] = NAN;
      continue;
    }
    fa[index] = A * m_ElementParam->m_NominalAmplitudes[index];
    // limit the absorption to 0.0-1.0, so that it's never <0
    if (m_ElementParam->m_SignFactors[index] == -1 &&
        m_ElementParam->m_absLinesLimit > 0.0 &&
        fa[index] > m_ElementParam->m_absLinesLimit) {
      fa[index] = m_ElementParam->m_absLinesLimit;
    }

    faes[index] = SNR * na[index];
  }
}

/**
 * \brief Adds to the model's flux, at each line not outside lambda range, the
 *value contained in the corresponding lambda for each catalog line.
 **/
void CLineModelElement::addToSpectrumModel(
    const CSpectrumSpectralAxis &modelspectralAxis,
    CSpectrumFluxAxis &modelfluxAxis,
    const CSpectrumFluxAxis &continuumfluxAxis, Float64 redshift,
    Int32 line_index) const {
  if (m_OutsideLambdaRange)
    return;

  for (Int32 index = 0; index != GetSize(); ++index) { // loop on the interval
    if (m_OutsideLambdaRangeList[index])
      continue;

    if (line_index != undefIdx && !(m_LineIsActiveOnSupport[line_index][index]))
      continue;

    for (Int32 i = m_StartNoOverlap[index]; i <= m_EndNoOverlap[index]; i++) {
      Float64 lambda = modelspectralAxis[i];
      Float64 Yi =
          getModelAtLambda(lambda, redshift, continuumfluxAxis[i], index);
      modelfluxAxis[i] += Yi;
      if (std::isnan(modelfluxAxis[i]))
        THROWG(INTERNAL_ERROR,
               Formatter() << "addToSpectrumModel has a NaN flux Line" << index
                           << ": ContinuumFlux " << continuumfluxAxis[i]
                           << ", ModelAtLambda Yi = " << Yi << " for range ["
                           << m_StartNoOverlap[index] << ", "
                           << m_EndNoOverlap[index] << "]");
    }
  }
  return;
}

void CLineModelElement::addToSpectrumModelDerivVel(
    const CSpectrumSpectralAxis &modelspectralAxis,
    CSpectrumFluxAxis &modelfluxAxis,
    const CSpectrumFluxAxis &continuumfluxAxis, Float64 redshift,
    bool emissionLine) const {
  if (m_OutsideLambdaRange)
    return;

  for (Int32 index = 0; index != GetSize(); ++index) {
    if (m_OutsideLambdaRangeList[index])
      continue;

    if ((emissionLine ^ IsEmission()))
      continue;

    Float64 A = m_ElementParam->m_FittedAmplitudes[index];
    if (std::isnan(A))
      continue;
    // THROWG(INTERNAL_ERROR,"FittedAmplitude cannot
    // be NAN");//to be uncommented

    for (Int32 i = m_StartNoOverlap[index]; i <= m_EndNoOverlap[index]; i++) {

      Float64 x = modelspectralAxis[i];
      Float64 mu = NAN;
      Float64 sigma = NAN;
      getObservedPositionAndLineWidth(index, redshift, mu, sigma, false);

      if (m_ElementParam->m_SignFactors[index] == -1)
        modelfluxAxis[i] += m_ElementParam->m_SignFactors[index] * A *
                            continuumfluxAxis[i] *
                            GetLineProfileDerivVel(*getLineProfile(index), x,
                                                   mu, sigma, IsEmission());
      else
        modelfluxAxis[i] += m_ElementParam->m_SignFactors[index] * A *
                            GetLineProfileDerivVel(*getLineProfile(index), x,
                                                   mu, sigma, IsEmission());
    }
  }
  return;
}

/**
 * \brief Returns the sum of the amplitude of each line on redshifted lambda.
 **/
Float64 CLineModelElement::getModelAtLambda(Float64 lambda, Float64 redshift,
                                            Float64 continuumFlux,
                                            Int32 line_index) const {
  if (m_OutsideLambdaRange)
    return 0.0;

  Float64 Yi = 0.0;

  Float64 x = lambda;

  for (Int32 index = 0; index != GetSize(); ++index) {
    if (m_OutsideLambdaRangeList[index])
      continue;

    if (line_index >= 0 && !m_LineIsActiveOnSupport[line_index][index])
      continue;

    Float64 A = m_ElementParam->m_FittedAmplitudes[index];
    if (std::isnan(A))
      continue;
    // THROWG(INTERNAL_ERROR,"FittedAmplitude cannot
    // be NAN");
    if (A < 0.)
      continue;

    Float64 fluxval = m_ElementParam->m_SignFactors[index] * A *
                      GetLineProfileAtRedshift(index, redshift, x);
    Yi += m_ElementParam->m_SignFactors[index] == -1 ? continuumFlux * fluxval
                                                     : fluxval;

    if (std::isnan(Yi))
      THROWG(INTERNAL_ERROR,
             Formatter() << "NaN fluxval for Line nb: " << index
                         << " and GetLineProfileAtRedshift: "
                         << GetLineProfileAtRedshift(index, redshift, x));
  }
  return Yi;
}

Float64 CLineModelElement::GetModelDerivAmplitudeAtLambda(
    Float64 lambda, Float64 redshift, Float64 continuumFlux) const {
  if (m_OutsideLambdaRange)
    return 0.0;

  Float64 Yi = 0.0;

  Float64 x = lambda;

  for (Int32 index = 0; index != GetSize(); ++index) { // loop on lines
    if (m_OutsideLambdaRangeList[index])
      continue;

    Float64 fluxval = m_ElementParam->m_SignFactors[index] *
                      m_ElementParam->m_NominalAmplitudes[index] *
                      GetLineProfileAtRedshift(index, redshift, x);
    Yi += m_ElementParam->m_SignFactors[index] == -1 ? continuumFlux * fluxval
                                                     : fluxval;
  }
  return Yi;
}

Float64
CLineModelElement::GetModelDerivVelAtLambda(Float64 lambda, Float64 redshift,
                                            Float64 continuumFlux) const {
  if (m_OutsideLambdaRange)
    return 0.0;

  Float64 Yi = 0.0;

  Float64 x = lambda;

  for (Int32 index = 0; index != GetSize(); ++index) { // loop on lines
    if (m_OutsideLambdaRangeList[index])
      continue;

    Float64 A = m_ElementParam->m_FittedAmplitudes[index];
    if (std::isnan(A))
      continue;
    // THROWG(INTERNAL_ERROR,"FittedAmplitude cannot
    // be NAN");
    if (A < 0.)
      continue;

    Float64 mu = NAN;
    Float64 sigma = NAN;
    getObservedPositionAndLineWidth(index, redshift, mu, sigma,
                                    false); // do not apply Lya asym offset

    const auto &profile = getLineProfile(index);
    bool isEmission = GetElementType() == CLine::EType::nType_Emission;
    Float64 lineprofile_derivVel =
        GetLineProfileDerivVel(*profile, x, mu, sigma, isEmission);
    Float64 fluxval =
        m_ElementParam->m_SignFactors[index] * A * lineprofile_derivVel;

    Yi += m_ElementParam->m_SignFactors[index] == -1 ? continuumFlux * fluxval
                                                     : fluxval;
  }
  return Yi;
}

Float64 CLineModelElement::GetModelDerivContinuumAmpAtLambda(
    Float64 lambda, Float64 redshift, Float64 continuumFluxUnscale) const {
  if (m_OutsideLambdaRange)
    return 0.0;

  Float64 Yi = 0.0;

  Float64 x = lambda;

  for (Int32 index = 0; index != GetSize(); ++index) { // loop on lines

    if (m_OutsideLambdaRangeList[index])
      continue;

    if (m_ElementParam->m_SignFactors[index] == 1)
      continue;

    Float64 A = m_ElementParam->m_FittedAmplitudes[index];
    if (std::isnan(A))
      continue;
    // THROWG(INTERNAL_ERROR,"FittedAmplitude cannot
    // be NAN");

    Yi += m_ElementParam->m_SignFactors[index] * continuumFluxUnscale * A *
          GetLineProfileAtRedshift(index, redshift, x);
  }
  return Yi;
}

/* Given the value of the partial deriv of the flux of this multiline at the
 * given lamda when The continuum is a variable of z
 */
Float64
CLineModelElement::GetModelDerivZAtLambda(Float64 lambda, Float64 redshift,
                                          Float64 continuumFlux,
                                          Float64 continuumFluxDerivZ) const {
  if (m_OutsideLambdaRange)
    return 0.0;

  Float64 Yi = 0.0;

  Float64 x = lambda;

  for (Int32 index = 0; index != GetSize(); ++index) { // loop on lines
    if (m_OutsideLambdaRangeList[index])
      continue;

    Float64 A = m_ElementParam->m_FittedAmplitudes[index];
    if (std::isnan(A))
      continue;
    // THROWG(INTERNAL_ERROR,"FittedAmplitude cannot
    // be NAN");

    Float64 mu = NAN;
    Float64 sigma = NAN;
    getObservedPositionAndLineWidth(index, redshift, mu, sigma,
                                    false); // do not apply Lya asym offset
    Float64 lambda_rest =
        GetObservedPosition(index, 0.0, false); // get restframe wavelentgh
    const auto &profile = getLineProfile(index);

    Float64 profile_derivz_val =
        lambda_rest * profile->GetLineProfileDerivX0(x, mu, sigma);

    Float64 fluxval =
        m_ElementParam->m_SignFactors[index] * A * profile_derivz_val;

    Yi += m_ElementParam->m_SignFactors[index] == -1
              ? continuumFlux * fluxval -
                    A * continuumFluxDerivZ *
                        profile->GetLineProfileVal(x, mu, sigma)
              : fluxval;
  }
  return Yi;
}

/**
 * \brief For lines inside lambda range, sets the flux to the continuum flux.
 **/
void CLineModelElement::initSpectrumModel(
    CSpectrumFluxAxis &modelfluxAxis,
    const CSpectrumFluxAxis &continuumfluxAxis, Int32 line_index) const {

  if (m_OutsideLambdaRange)
    return;

  for (Int32 index = 0; index != GetSize(); ++index) { // loop on the interval
    if (m_OutsideLambdaRangeList[index])
      continue;

    if (line_index != undefIdx && !(m_LineIsActiveOnSupport[line_index][index]))
      continue;

    for (Int32 i = m_StartNoOverlap[index]; i <= m_EndNoOverlap[index]; i++)
      modelfluxAxis[i] = continuumfluxAxis[i];
  }
  return;
}

/**
 * \brief For lines inside lambda range, sets the flux to the polynomial.
 **/
void CLineModelElement::initSpectrumModelPolynomial(
    CSpectrumFluxAxis &modelfluxAxis, const CSpectrumSpectralAxis &spcAxis,
    Int32 line_index) const {

  if (m_OutsideLambdaRange)
    return;

  for (Int32 index = 0; index != GetSize(); ++index) { // loop on the interval
    if (m_OutsideLambdaRangeList[index])
      continue;

    if (line_index != undefIdx && !(m_LineIsActiveOnSupport[line_index][index]))
      continue;

    for (Int32 i = m_StartNoOverlap[index]; i <= m_EndNoOverlap[index]; i++)
      modelfluxAxis[i] =
          m_ElementParam->m_ampOffsetsCoeffs.getValue(spcAxis[i]);
  }
  return;
}

/**
 * \brief If the fitted amplitude of line with id line_id is above the limit,
 *sets it to either that limit or zero, whichever is greater.
 **/
bool CLineModelElement::LimitFittedAmplitude(Int32 line_index, Float64 limit) {

  bool limited = false;
  if (GetFittedAmplitude(line_index) > limit) {
    limited = true;
    m_ElementParam->m_FittedAmplitudes[line_index] = std::max(0.0, limit);

    // now update the amplitude of the other lines
    Float64 amplitudeRef = GetFittedAmplitude(line_index) /
                           m_ElementParam->m_NominalAmplitudes[line_index];
    for (Int32 index = 0; index != GetSize(); ++index) {
      m_ElementParam->m_FittedAmplitudes[index] =
          m_ElementParam->m_NominalAmplitudes[index] * amplitudeRef;
    }
  }
  return limited;
}

void CLineModelElement::debug(std::ostream &os) const {
  os << m_ElementParam->m_fittingGroupInfo << "\t";
  os << getVelocityAbsorption() << "\t" << getVelocityEmission() << "\n";

/*   for (Int32 line_index = 0; line_index != GetSize(); ++line_index) {
    os << GetLineName(line_index) << "\t" << GetNominalAmplitude(line_index)
       << "\t" << getLineProfile(line_index)->GetName() << "\t"
       << GetSignFactor(line_index) << "\n";
  }
 */}

void CLineModelElement::dumpElement(std::ostream &os) const {
  //  debug(os); // to dump lines info
  os << "m_OutsideLambdaRange\t" << m_OutsideLambdaRange << "\n";
  os << "m_fittingGroupInfo\t" << m_ElementParam->m_fittingGroupInfo << "\n";
  os << "m_ElementType\t" << CLine::ETypeString.at(GetElementType()) << "\n";

  os << "m_VelocityEmission\t" << getVelocityEmission() << "\n";
  os << "m_VelocityAbsorption\t" << getVelocityAbsorption() << "\n";
  os << "m_OutsideLambdaRangeOverlapThreshold\t"
     << m_OutsideLambdaRangeOverlapThreshold << "\n";
  os << "m_sumCross\t" << m_ElementParam->m_sumCross << "\n";
  os << "m_sumGauss\t" << m_ElementParam->m_sumGauss << "\n";
  os << "m_dtmFree\t" << m_ElementParam->m_dtmFree << "\n";

  os << "m_absLinesLimit\t" << m_ElementParam->m_absLinesLimit << "\n";

  os << "m_LineIsActiveOnSupport \n";

  for (Int32 i = 0; i != GetSize(); ++i)
    for (Int32 j = 0; j != GetSize(); ++j)
      os << i << "\t" << j << "\t" << m_LineIsActiveOnSupport[i][j] << "\n";

  os << "\n";
  os << "Line \t OutsideLR\t Sign \t Amp \t AmpErrSigma \t NomAmp \t StartNO "
        "\t EndNO \t "
        "StartTheo \t EndTheo\n";

  for (Int32 i = 0; i != GetSize(); ++i) {
    os << i << "\t " << m_OutsideLambdaRangeList[i] << "\t "
       << m_ElementParam->m_SignFactors[i] << "\t "
       << m_ElementParam->m_FittedAmplitudes[i] << "\t"
       << m_ElementParam->m_FittedAmplitudeErrorSigmas[i] << "\t"
       << m_ElementParam->m_NominalAmplitudes[i] << "\t" << m_StartNoOverlap[i]
       << "\t" << m_EndNoOverlap[i] << "\t" << m_StartTheoretical[i] << "\t"
       << m_EndTheoretical[i] << "\n";
  }

  os << "\n";
  os << "m_asymLineIndices \n";
  for (Int32 i = 0; i < m_ElementParam->m_asymLineIndices.size(); i++)
    os << i << "\t" << m_ElementParam->m_asymLineIndices[i] << "\n";
}

Int32 CLineModelElement::computeCrossProducts(
    Float64 redshift, const CSpectrumSpectralAxis &spectralAxis,
    const CSpectrumFluxAxis &noContinuumfluxAxis,
    const CSpectrumFluxAxis &continuumfluxAxis, Int32 line_index) {

  const CSpectrumNoiseAxis &error = noContinuumfluxAxis.GetError();
  auto &nominalAmplitudes = m_ElementParam->m_NominalAmplitudes;
  Float64 y = 0.0;
  Float64 x = 0.0;
  Float64 yg = 0.0;
  Float64 c = 1.0;

  Float64 err2 = 0.0;
  Int32 num = 0;

  for (Int32 index = 0; index != GetSize(); ++index) { // loop on lines
    if (m_OutsideLambdaRangeList[index])
      continue;

    if (line_index != undefIdx && !isLineActiveOnSupport(line_index, index))
      continue;

    for (Int32 i = getStartNoOverlap(index); i <= getEndNoOverlap(index); i++) {
      c = continuumfluxAxis[i];
      y = noContinuumfluxAxis[i];
      x = spectralAxis[i];

      yg = 0.0;

      for (Int32 index2 = 0; index2 != GetSize();
           ++index2) { // loop for the signal synthesis
        if (m_OutsideLambdaRangeList[index2] ||
            !m_LineIsActiveOnSupport[index][index2])
          continue;

        Int32 sf = getSignFactor(index2);
        if (sf == -1) {
          yg += sf * c * nominalAmplitudes[index2] *
                GetLineProfileAtRedshift(index2, redshift, x);
        } else {
          yg += sf * nominalAmplitudes[index2] *
                GetLineProfileAtRedshift(index2, redshift, x);
        }
      }
      num++;
      err2 = 1.0 / (error[i] * error[i]);
      m_ElementParam->m_dtmFree += yg * y * err2;
      m_ElementParam->m_sumGauss += yg * yg * err2;
    }
  }

  return num;
}

/*
void CLineModelElement::fitAmplitude(
    Float64 redshift, const CSpectrumSpectralAxis &spectralAxis,
    const CSpectrumFluxAxis &noContinuumfluxAxis,
    const CSpectrumFluxAxis &continuumfluxAxis, Int32 line_index) {


  computeCrossProducts(redshift, spectralAxis, noContinuumfluxAxis,
                       continuumfluxAxis, line_index);
}
*/
