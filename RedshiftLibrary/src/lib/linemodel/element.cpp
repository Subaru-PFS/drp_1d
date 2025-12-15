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

#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/linemodel/element.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace std;
using namespace NSEpic;

/**
 * \brief Constructs the object setting members according to arguments and
 *defaults.
 **/
CLineModelElement::CLineModelElement(
    const TLineModelElementParam_ptr elementParam, Float64 maxDistanceToLine,
    Int32 minSamplesNumberForLineFit)
    : m_ElementParam(std::move(elementParam)),
      m_maxDistanceToLine(maxDistanceToLine),
      m_minSamplesNumberForLineFit(minSamplesNumberForLineFit),
      m_OutsideLambdaRange(true), m_size(m_ElementParam->size()){};

void TLineModelElementParam::resetFittingParams() {
  // init the fitted amplitude values and related variables
  m_FittedAmplitudes.assign(size(), NAN);
  m_FittedAmplitudesStd.assign(size(), NAN);
  m_fittingGroupInfo = undefStr;
  m_ampOffsetsCoeffs = CPolynomCoeffs();

  m_sumGauss = NAN;
  m_sumCross = NAN;
  m_dtmFree = NAN;
}

Int32 TLineModelElementParam::getLineIndex(Int32 line_id) const {
  Int32 index = undefIdx;

  auto const &it = m_LinesIds.find(line_id);
  if (it != m_LinesIds.cend())
    index = it->second;
  return index;
}

Int32 TLineModelElementParam::getLineIndex(
    const std::string &LineTagStr) const {
  Int32 line_index = undefIdx;

  auto it =
      std::find_if(m_Lines.cbegin(), m_Lines.cend(), [&LineTagStr](auto &line) {
        auto const &name = line.GetName();
        std::size_t foundstra = name.find(LineTagStr.c_str());
        return foundstra != std::string::npos;
      });
  if (it != m_Lines.cend())
    line_index = it - m_Lines.cbegin();

  return line_index;
}

// redirecting to the lsf method for computing instrument responce
/**
 * Get instrumental response (including source response) from LSF
 * combine quadratically with the instrinsic width of the Line itself. The line
 * width in this case represents to the velocity
 * */
Float64 CLineModelElement::GetLineWidth(Float64 redshiftedlambda) const {
  const Float64 c = SPEED_OF_LIGHT_IN_VACCUM;
  Float64 v = m_ElementParam->getVelocity();

  if (!m_LSF)
    THROWG(ErrorCode::INTERNAL_ERROR, "lsf object is not initailized.");
  Float64 instrumentSigma = m_LSF->GetWidth(redshiftedlambda);

  Float64 velocitySigma = v / c * redshiftedlambda;
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
    THROWG(ErrorCode::INTERNAL_ERROR, Formatter()
                                          << "Invalid LSF type "
                                          << m_ElementParam->m_LineWidthType);
  }

  Float64 sigma =
      sqrt(instrumentSigma * instrumentSigma + velocitySigma * velocitySigma);
  return sigma;
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
std::pair<Float64, Float64> CLineModelElement::GetContinuumAtCenterProfile(
    Int32 line_id, const CSpectrumSpectralAxis &spectralAxis, Float64 redshift,
    const CSpectrumFluxAxis &continuumfluxAxis,
    bool enableAmplitudeOffsets) const {
  Float64 mu = GetObservedPosition(line_id, redshift);

  Int32 IdxCenterProfile = spectralAxis.GetIndexAtWaveLength(mu);
  if (IdxCenterProfile < 0 ||
      IdxCenterProfile > continuumfluxAxis.GetSamplesCount() - 1) {
    return std::make_pair(NAN, NAN);
  }

  Float64 cont = continuumfluxAxis[IdxCenterProfile];
  Float64 contStd = NAN;
  if (enableAmplitudeOffsets) {
    auto const &polyCoeffs = getElementParam()->GetPolynomCoeffs();
    contStd = std::sqrt(polyCoeffs.getVariance(spectralAxis[IdxCenterProfile]));
  }
  return std::make_pair(cont, contStd);
}

void CLineModelElement::EstimateSupport(
    Int32 line_index, const CSpectrumSpectralAxis &spectralAxis,
    Float64 redshift, const TFloat64Range &lambdaRange, Float64 max_offset) {
  Float64 mu = GetObservedPosition(line_index, redshift);
  if (!m_LSF->checkAvailability(mu)) {
    m_OutsideLambdaRangeList[line_index] = true;
    return;
  }
  Float64 const sigma = GetLineWidth(mu);
  Float64 const max_offset_angstrom =
      (max_offset / SPEED_OF_LIGHT_IN_VACCUM) * mu;
  Float64 winsize =
      getElementParam()->getLineProfile(line_index)->GetNSigmaSupport() * sigma;
  winsize += 2 * max_offset_angstrom;
  TInt32Range supportRange =
      EstimateIndexRange(spectralAxis, mu, lambdaRange, winsize);

  m_range[line_index] = supportRange;
  m_rangeNoOverlap[line_index] = supportRange;

  EstimateLineVisbility(line_index, spectralAxis, mu, sigma,
                        max_offset_angstrom);
}

void CLineModelElement::EstimateLineVisbility(
    Int32 line_index, const CSpectrumSpectralAxis &spectralAxis,
    Float64 line_lambda, Float64 sigma, Float64 max_offset) {

  auto const &supportRange = m_range[line_index];
  if (supportRange.GetLength() < 0) {
    // in this case the line is completely outside the
    // lambdarange
    m_OutsideLambdaRangeList[line_index] = true;
    return;
  }

  // in this case the line is completely inside the lambdarange or with
  // partial overlap
  Int32 const nsupport = supportRange.GetLength() + 1;
  TFloat64List distance(nsupport);
  auto const &wave_iter =
      spectralAxis.GetSamplesVector().begin() + supportRange.GetBegin();
  std::transform(
      wave_iter, wave_iter + nsupport, distance.begin(),
      [line_lambda](Float64 lambda) { return std::abs(lambda - line_lambda); });
  auto const &min_distance =
      *std::min_element(distance.begin(), distance.end());

  if (min_distance > (m_maxDistanceToLine * sigma + max_offset) ||
      (nsupport < m_minSamplesNumberForLineFit)) {
    m_OutsideLambdaRangeList[line_index] = true;
  } else {
    m_OutsideLambdaRangeList[line_index] = false;
  }

  return;
}

/**
 * \brief Returns the index range for a given window size (Angstrom)
 **/
TInt32Range CLineModelElement::EstimateIndexRange(
    const CSpectrumSpectralAxis &spectralAxis, Float64 mu,
    const TFloat64Range &lambdaRange, Float64 winsizeAngstrom) {

  Float64 const winsize = winsizeAngstrom;
  Float64 const lambda_start = mu - winsize / 2.0;
  Float64 const lambda_end = mu + winsize / 2.0;

  TFloat64Range elementLambdaRange{lambda_start, lambda_end};
  auto has_intersection = elementLambdaRange.IntersectWith(lambdaRange);
  TInt32Range supportRange;
  if (!has_intersection || elementLambdaRange.GetIsEmpty()) {
    supportRange = {0, -1};
    return supportRange;
  }

  try {
    supportRange =
        spectralAxis.GetIndexRangeAtWaveLengthRange(elementLambdaRange);
  } catch (const AmzException &exception) {
    if (exception.getErrorCode() == ErrorCode::IE_CRANGE_NO_INTERSECTION) {
      Int32 imin =
          spectralAxis.GetIndexAtWaveLength(elementLambdaRange.GetBegin());
      TInt32Range supportRange{imin, imin - 1};
      return supportRange;
    } else {
      throw exception;
    }
  }

  return supportRange;
}

void CLineModelElement::computeOutsideLambdaRange() {
  m_OutsideLambdaRange = true;
  for (auto const &outside_lambda_range_id : m_OutsideLambdaRangeList)
    if (!outside_lambda_range_id) {
      m_OutsideLambdaRange = false;
      break;
    }
}

/**
 * \brief Limits each m_Lines element within the argument lambdaRange, and sets
 *the m_FittedAmplitudes to -1. Sets the global outside lambda range. Inits the
 *fitted amplitude values.
 **/

void CLineModelElement::initSupport(const CSpectrumSpectralAxis &spectralAxis,
                                    Float64 redshift,
                                    const TFloat64Range &lambdaRange,
                                    Float64 max_offset) {

  Int32 nLines = GetSize();
  m_OutsideLambdaRange = true;
  m_rangeNoOverlap.assign(nLines, TInt32Range{undefIdx, undefIdx});
  m_range.assign(nLines, TInt32Range{undefIdx, undefIdx});
  m_OutsideLambdaRangeList.assign(nLines, true);
  m_LineIsActiveOnSupport.assign(nLines, TBoolList(nLines, false));
  m_sortedLineIndices.clear();

  for (Int32 index = 0; index != nLines; ++index) {
    EstimateSupport(index, spectralAxis, redshift, lambdaRange, max_offset);
    // set the lines active on their own support
    m_LineIsActiveOnSupport[index][index] = true;
  }
  computeOutsideLambdaRange();
}

bool CLineModelElement::mergeIfOverlapping(Int32 i, Int32 j) {

  if (!m_rangeNoOverlap[i].unionWith(m_rangeNoOverlap[j]))
    return false;

  // deactivate j
  m_rangeNoOverlap[j] = {m_rangeNoOverlap[i].GetEnd(),
                         m_rangeNoOverlap[i].GetEnd() - 1};
  return true;
}

void CLineModelElement::sortLinesByLeftIndex() {
  Int32 nLines = GetSize();
  m_sortedLineIndices.reserve(nLines);

  // list visible lines
  for (Int32 index = 0; index != nLines; ++index)
    if (!m_OutsideLambdaRangeList[index])
      m_sortedLineIndices.push_back(index);

  // then sort by left most index
  std::sort(m_sortedLineIndices.begin(), m_sortedLineIndices.end(),
            [this](Int32 l, Int32 r) {
              return m_rangeNoOverlap[l] < m_rangeNoOverlap[r];
            });
}

void CLineModelElement::mergeOverlapingLines() {
  auto alreadyMerged = [this](Int32 index) {
    return m_rangeNoOverlap[index].GetLength() < 0;
  };

  // assumes lines are sorted by left boundary
  for (auto iter_left = m_sortedLineIndices.begin();
       iter_left != m_sortedLineIndices.end(); ++iter_left) {
    if (alreadyMerged(*iter_left))
      continue;

    for (auto iter_right = iter_left + 1;
         iter_right != m_sortedLineIndices.end(); ++iter_right) {
      if (alreadyMerged(*iter_right))
        continue;

      if (mergeIfOverlapping(*iter_left, *iter_right)) {
        propagateOverlap(*iter_left, *iter_right);
      }
    }
  }
}

void CLineModelElement::prepareSupport(
    const CSpectrumSpectralAxis &spectralAxis, Float64 redshift,
    const TFloat64Range &lambdaRange, Float64 max_offset) {

  initSupport(spectralAxis, redshift, lambdaRange, max_offset);

  if (!IsOutsideLambdaRange()) {
    sortLinesByLeftIndex();
    mergeOverlapingLines();
    if (detectRemainingOverlaps())
      THROWG(ErrorCode::INTERNAL_ERROR, "Remaining overlaps");
  }
}

bool CLineModelElement::detectRemainingOverlaps() {
  Int32 nLines = GetSize();

  for (Int32 i = 0; i != nLines; ++i) {
    if (m_OutsideLambdaRangeList[i])
      continue;

    for (Int32 j = i + 1; j != nLines; ++j) {
      if (m_OutsideLambdaRangeList[j])
        continue;

      if (m_rangeNoOverlap[i].HasIntersectionWith(m_rangeNoOverlap[j]))
        return true;
    }
  }
  return false;
}

void CLineModelElement::propagateOverlap(Int32 i, Int32 j) {
  m_LineIsActiveOnSupport[i][j] = true;
  m_LineIsActiveOnSupport[j][i] = true;

  Int32 nLines = GetSize();

  // Propagate through i
  for (Int32 i2 = 0; i2 != nLines; ++i2) {
    if (m_OutsideLambdaRangeList[i2])
      continue;
    if (m_LineIsActiveOnSupport[i][i2]) {
      m_LineIsActiveOnSupport[i2][j] = true;
      m_LineIsActiveOnSupport[j][i2] = true;
    }
  }

  // Propagate through j
  for (Int32 j2 = 0; j2 != nLines; ++j2) {
    if (m_OutsideLambdaRangeList[j2])
      continue;
    if (m_LineIsActiveOnSupport[j][j2]) {
      m_LineIsActiveOnSupport[j2][i] = true;
      m_LineIsActiveOnSupport[i][j2] = true;
    }
  }
}

/**
 * \brief Creates an empty list of ranges as the return value. If not
 *m_OutsideLambdaRange, for each m_Lines element which is also not outside
 *lambda range, add its support to the return value.
 **/
TInt32RangeList CLineModelElement::getSupportNoOverlap() const {
  TInt32RangeList support;
  if (m_OutsideLambdaRange)
    return support;

  for (Int32 index = 0; index != GetSize(); ++index) {
    if (m_OutsideLambdaRangeList[index])
      continue;

    support.push_back(m_rangeNoOverlap[index]);
  }
  return support;
}

TInt32RangeList CLineModelElement::getSupport() const {
  TInt32RangeList support;

  if (m_OutsideLambdaRange)
    return support;

  for (Int32 index = 0; index != GetSize(); ++index) {
    if (m_OutsideLambdaRangeList[index])
      continue;

    support.push_back(m_range[index]);
  }

  return support;
}

/**
 * \brief Calls GetLineWidth using the arguments and a calculated argument mu.
 **/
std::pair<Float64, Float64> CLineModelElement::getObservedPositionAndLineWidth(
    Float64 redshift, Int32 index, bool doAsymfitdelta) const {

  Float64 const mu = (index == undefIdx)
                         ? GetObservedPosition(redshift, doAsymfitdelta)
                         : GetObservedPosition(index, redshift, doAsymfitdelta);

  if (!m_LSF->checkAvailability(mu))
    THROWG(ErrorCode::INTERNAL_ERROR,
           "Line position does not belong to LSF range");
  Float64 const sigma = GetLineWidth(mu);
  return std::make_pair(mu, sigma);
}

/**
 * \brief Get the observed position of the sub-element with id line_id for a
 *given redshift
 **/
Float64 CLineModelElement::GetObservedPosition(Int32 index, Float64 redshift,
                                               bool doAsymfitdelta) const {
  Float64 dzOffset =
      m_ElementParam->getLambdaOffset(index) / SPEED_OF_LIGHT_IN_VACCUM;

  auto const &line = m_ElementParam->GetLines()[index];
  Float64 mu = line.GetPosition() * (1 + redshift) * (1 + dzOffset);

  // deals with delta of asym profile
  if (doAsymfitdelta) {
    mu -= line.GetProfile()->GetDelta();
  }
  return mu;
}

Float64 CLineModelElement::GetObservedPosition(Float64 redshift,
                                               bool doAsymfitdelta) const {
  Float64 mu = 0.0;
  Int32 nvalid = 0;
  for (Int32 index = 0; index < GetSize(); ++index) {
    if (IsOutsideLambdaRangeLine(index))
      continue;
    mu += GetObservedPosition(index, redshift, doAsymfitdelta);
    ++nvalid;
  }
  if (nvalid == 0)
    THROWG(ErrorCode::INTERNAL_ERROR,
           "Trying to get mean position of an element with "
           "all lines outside range");
  mu /= nvalid;
  return mu;
}

/**
 * \brief Returns the line profile of the sub-element with id line_id at
 *wavelength x, for a given redshift.
 **/
Float64 CLineModelElement::GetLineProfileAtRedshift(Int32 index,
                                                    Float64 redshift,
                                                    Float64 x) const {
  auto const &[mu, sigma] =
      getObservedPositionAndLineWidth(redshift, index,
                                      false); // do not apply Lya asym offset
  if (sigma == 0.0) {
    Flag.warning(WarningCode::NULL_LINES_PROFILE,
                 Formatter() << "null line width sigma, for line "
                             << getElementParam()->GetLineName(index));
    return 0.0;
  }

  auto const &profile = getElementParam()->getLineProfile(index);
  return profile->GetLineProfileVal(x, mu, sigma);
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
  if (m_OutsideLambdaRange || getElementParam()->isNotFittable())
    return;

  for (Int32 index = 0; index != GetSize(); ++index) { // loop on the interval
    if (m_OutsideLambdaRangeList[index])
      continue;

    if (line_index != undefIdx && !(m_LineIsActiveOnSupport[line_index][index]))
      continue;

    for (Int32 i : m_rangeNoOverlap[index]) {
      Float64 lambda = modelspectralAxis[i];
      Float64 Yi =
          getModelAtLambda(lambda, redshift, continuumfluxAxis[i], index);
      modelfluxAxis[i] += Yi;
      if (std::isnan(modelfluxAxis[i]))
        THROWG(ErrorCode::INTERNAL_ERROR,
               Formatter() << "NaN flux Line: "
                           << getElementParam()->GetLineName(index)
                           << ", ContinuumFlux " << continuumfluxAxis[i]
                           << ", ModelAtLambda Yi = " << Yi << "for range "
                           << m_rangeNoOverlap[index]);
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

    if ((emissionLine != getElementParam()->IsEmission()))
      continue;

    Float64 A = m_ElementParam->m_FittedAmplitudes[index];
    if (std::isnan(A))
      THROWG(ErrorCode::INTERNAL_ERROR, "FittedAmplitude cannot be NAN");

    for (Int32 i : m_rangeNoOverlap[index]) {

      Float64 const x = modelspectralAxis[i];
      auto const &[mu, sigma] =
          getObservedPositionAndLineWidth(redshift, index, false);

      if (m_ElementParam->m_SignFactors[index] == -1)
        modelfluxAxis[i] +=
            m_ElementParam->m_SignFactors[index] * A * continuumfluxAxis[i] *
            m_ElementParam->GetLineProfileDerivVel(index, x, mu, sigma);
      else
        modelfluxAxis[i] +=
            m_ElementParam->m_SignFactors[index] * A *
            m_ElementParam->GetLineProfileDerivVel(index, x, mu, sigma);
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
      THROWG(ErrorCode::INTERNAL_ERROR, "FittedAmplitude cannot be NAN");
    if (A <= 0.)
      continue;

    Float64 fluxval = m_ElementParam->m_SignFactors[index] * A *
                      GetLineProfileAtRedshift(index, redshift, x);
    Yi += m_ElementParam->m_SignFactors[index] == -1 ? continuumFlux * fluxval
                                                     : fluxval;

    if (std::isnan(Yi))
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "NaN fluxval for Line: "
                         << getElementParam()->GetLineName(index)
                         << ", amplitude: " << A << ", line profile:"
                         << GetLineProfileAtRedshift(index, redshift, x)
                         << ", continnum: " << continuumFlux);
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
    if (m_ElementParam->m_NominalAmplitudes[index] == 0.0)
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

  for (Int32 index = 0; index != GetSize(); ++index) { // loop on lines
    if (m_OutsideLambdaRangeList[index])
      continue;

    Float64 const A = m_ElementParam->m_FittedAmplitudes[index];
    if (std::isnan(A))
      THROWG(ErrorCode::INTERNAL_ERROR, "FittedAmplitude cannot be NAN");
    if (A < 0.)
      continue;

    auto const &[mu, sigma] =
        getObservedPositionAndLineWidth(redshift, index,
                                        false); // do not apply Lya asym offset

    Float64 const lineprofile_derivVel =
        m_ElementParam->GetLineProfileDerivVel(index, lambda, mu, sigma);
    Float64 const fluxval =
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
      THROWG(ErrorCode::INTERNAL_ERROR, "FittedAmplitude cannot be NAN");
    if (A <= 0.0)
      continue;
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

  for (Int32 index = 0; index != GetSize(); ++index) { // loop on lines
    if (m_OutsideLambdaRangeList[index])
      continue;

    Float64 const A = m_ElementParam->m_FittedAmplitudes[index];
    if (std::isnan(A))
      THROWG(ErrorCode::INTERNAL_ERROR, "FittedAmplitude cannot be NAN");

    auto const &[mu, sigma] =
        getObservedPositionAndLineWidth(redshift, index,
                                        false); // do not apply Lya asym offset
    Float64 const lambda_rest =
        GetObservedPosition(index, 0.0, false); // get restframe wavelentgh
    auto const &profile = getElementParam()->getLineProfile(index);

    Float64 const profile_derivz_val =
        lambda_rest * profile->GetLineProfileDerivX0(lambda, mu, sigma);

    Float64 const fluxval =
        m_ElementParam->m_SignFactors[index] * A * profile_derivz_val;

    Yi += m_ElementParam->m_SignFactors[index] == -1
              ? continuumFlux * fluxval -
                    A * continuumFluxDerivZ *
                        profile->GetLineProfileVal(lambda, mu, sigma)
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

    for (Int32 i : m_rangeNoOverlap[index])
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

    for (Int32 i : m_rangeNoOverlap[index])
      modelfluxAxis[i] =
          m_ElementParam->m_ampOffsetsCoeffs.getValue(spcAxis[i]);
  }
  return;
}

void CLineModelElement::debug(std::ostream &os) const {}

void CLineModelElement::dumpElement(std::ostream &os) const {
  //  debug(os); // to dump lines info
  os << "m_OutsideLambdaRange\t" << m_OutsideLambdaRange << "\n";
  os << "m_fittingGroupInfo\t" << m_ElementParam->m_fittingGroupInfo << "\n";
  os << "m_ElementType\t"
     << CLine::ETypeString.at(getElementParam()->GetElementType()) << "\n";

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
       << m_ElementParam->m_FittedAmplitudesStd[i] << "\t"
       << m_ElementParam->m_NominalAmplitudes[i] << "\t"
       << m_rangeNoOverlap[i].GetBegin() << "\t" << m_rangeNoOverlap[i].GetEnd()
       << "\t" << m_range[i].GetBegin() << "\t" << m_range[i].GetEnd() << "\n";
  }

  os << "\n";
  os << "m_asymLineIndices \n";
  for (Int32 i = 0; i < ssize(m_ElementParam->m_asymLineIndices); i++)
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

    for (Int32 i : m_rangeNoOverlap[index]) {
      c = continuumfluxAxis[i];
      y = noContinuumfluxAxis[i];
      x = spectralAxis[i];

      yg = 0.0;

      for (Int32 index2 = 0; index2 != GetSize();
           ++index2) { // loop for the signal synthesis
        if (m_OutsideLambdaRangeList[index2] ||
            !m_LineIsActiveOnSupport[index][index2])
          continue;
        Int32 sf = m_ElementParam->getSignFactor(index2);
        Float64 amp = nominalAmplitudes[index2];
        if (sf == -1)
          amp *= -c;
        if (amp == 0.0)
          continue;
        yg += amp * GetLineProfileAtRedshift(index2, redshift, x);
      }
      num++;
      err2 = 1.0 / (error[i] * error[i]);
      m_ElementParam->m_dtmFree += yg * y * err2;
      m_ElementParam->m_sumGauss += yg * yg * err2;
    }
  }

  return num;
}
