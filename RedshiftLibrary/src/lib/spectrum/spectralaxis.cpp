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
#include <cmath>

#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/line/airvacuum.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"

using namespace NSEpic;
using namespace std;

/**
 * Constructor, flags log scale when set.
 */
CSpectrumSpectralAxis::CSpectrumSpectralAxis(Int32 n) : CSpectrumAxis(n) {}

/**
 * Constructor, flags log scale when set.
 */
CSpectrumSpectralAxis::CSpectrumSpectralAxis(const TFloat64List &samples,
                                             std::string const &AirVacuum)
    : CSpectrumAxis(samples) {
  convertToVacuum(AirVacuum);
}

CSpectrumSpectralAxis::CSpectrumSpectralAxis(TFloat64List &&samples,
                                             std::string const &AirVacuum)
    : CSpectrumAxis(std::move(samples)) {
  convertToVacuum(AirVacuum);
}

// only used by client
CSpectrumSpectralAxis::CSpectrumSpectralAxis(const Float64 *samples, Int32 n,
                                             std::string const &AirVacuum)
    : CSpectrumAxis(samples, n) {
  convertToVacuum(AirVacuum);
}

void CSpectrumSpectralAxis::convertToVacuum(std::string const &AirVacuum) {
  if (AirVacuum != "") {
    m_Samples = CAirVacuumConverter::Get(AirVacuum)->AirToVac(m_Samples);
    Log.LogInfo(
        Formatter()
        << "SpectralAxis converted from air to vacuum using translation from: "
        << AirVacuum);
  }
}

/**
 * Constructor, shifts origin along direction an offset distance.
 */
CSpectrumSpectralAxis::CSpectrumSpectralAxis(
    const CSpectrumSpectralAxis &origin, Float64 wavelengthOffset,
    EShiftDirection direction)
    : CSpectrumAxis(origin.GetSamplesCount()) {
  ShiftByWaveLength(origin, wavelengthOffset, direction);
}

CSpectrumSpectralAxis::CSpectrumSpectralAxis(Int32 n, Float64 value)
    : CSpectrumAxis(n, value) {}

CSpectrumSpectralAxis &CSpectrumSpectralAxis::operator*=(const Float64 op) {
  CSpectrumAxis::operator*=(op);
  if (op < 0 && !indeterminate(m_isSorted))
    m_isSorted = m_isSorted ? static_cast<tribool>(false) : indeterminate;
  if (!op)
    m_isSorted = (GetSamplesCount() < 2);
  return *this;
}
CSpectrumSpectralAxis &CSpectrumSpectralAxis::operator/=(const Float64 op) {
  CSpectrumAxis::operator/=(op);
  if (op < 0 && !indeterminate(m_isSorted))
    m_isSorted = m_isSorted ? static_cast<tribool>(false) : indeterminate;
  if (!op)
    m_isSorted = (GetSamplesCount() < 2);
  return *this;
}
CSpectrumSpectralAxis
CSpectrumSpectralAxis::MaskAxis(const TFloat64List &mask) const {
  CSpectrumSpectralAxis spc_axis = CSpectrumAxis::MaskAxis(mask);
  spc_axis.m_isSorted = m_isSorted;
  if (spc_axis.GetSamplesCount() < 2)
    spc_axis.m_isSorted = true;
  spc_axis.m_isLogSampled = indeterminate;
  return spc_axis;
}

void CSpectrumSpectralAxis::SetSize(Int32 s) {
  Float64 sample_size = GetSamplesCount();
  CSpectrumAxis::SetSize(s);
  if (s < 2)
    m_isSorted = true;
  else {
    if (sample_size < s)
      resetAxisProperties();
  }
}
/**
 * Shift current axis the input offset in the input direction.
 */
void CSpectrumSpectralAxis::ShiftByWaveLength(Float64 wavelengthOffset,
                                              EShiftDirection direction) {
  ShiftByWaveLength(*this, wavelengthOffset, direction);
}

/**
 * Copy the input axis samples and shift the axis the specified offset in the
 * specidifed direction.
 */
void CSpectrumSpectralAxis::ShiftByWaveLength(
    const CSpectrumSpectralAxis &origin, Float64 wavelengthOffset,
    EShiftDirection direction) {
  if (wavelengthOffset < 0.)
    THROWG(INTERNAL_ERROR, "wavelengthOffset can not be negative");
  if (!(direction == nShiftForward || direction == nShiftBackward))
    THROWG(INTERNAL_ERROR, "Unknown shift direction");

  *this = origin;
  if (wavelengthOffset == 0.0)
    return;

  if (direction == nShiftForward) {
    operator*=(wavelengthOffset);
  } else if (direction == nShiftBackward) {
    operator/=(wavelengthOffset);
  }
}

void CSpectrumSpectralAxis::ApplyOffset(Float64 wavelengthOffset) {
  Int32 nSamples = m_Samples.size();
  for (Int32 i = 0; i < nSamples; i++) {
    m_Samples[i] += wavelengthOffset;
  }
  if (!wavelengthOffset)
    m_isLogSampled = false;
}

/**
 * Return the wavelength interval between two consecutive samples.
 */
Float64 CSpectrumSpectralAxis::GetResolution(Float64 atWavelength) const {
  if (GetSamplesCount() < 2)
    return 0.0;

  Int32 i = 0;

  if (atWavelength >= 0.0) {
    i = GetIndexAtWaveLength(atWavelength);

    if (i > m_Samples.size() - 1)
      i = m_Samples.size() - 1;
    if (i < 1)
      i = 1;

    return m_Samples[i] - m_Samples[i - 1];
  } else {
    return m_Samples[1] - m_Samples[0];
  }

  return 0.0;
}

/**
 *
 */
Float64 CSpectrumSpectralAxis::GetMeanResolution() const {
  if (GetSamplesCount() < 2)
    return 0.0;

  Float64 resolution = 0.0;
  Int32 nsum = 0;
  Int32 nSamples = m_Samples.size() - 1;
  for (Int32 i = 0; i < nSamples; i++) {
    resolution += (m_Samples[i + 1] - m_Samples[i]);
    nsum++;
  }

  if (nsum > 0) {
    resolution /= nsum;
  }

  return resolution;
}

/**
 *
 */
TLambdaRange CSpectrumSpectralAxis::GetLambdaRange() const {
  if (m_Samples.size() < 2)
    return TLambdaRange(0.0, 0.0);
  return TLambdaRange(m_Samples[0], m_Samples[m_Samples.size() - 1]);
}

/**
 *
 */
void CSpectrumSpectralAxis::GetMask(const TFloat64Range &lambdaRange,
                                    CMask &mask) const {
  TFloat64Range range = lambdaRange;

  mask.SetSize(m_Samples.size());

  // weight = Spectrum over lambdarange flag
  for (Int32 i = 0; i < m_Samples.size(); i++) {
    mask[i] = (Mask)0;
    // If this sample is somewhere in a valid lambdaRande, tag weight with 1

    if (m_Samples[i] >= range.GetBegin() && m_Samples[i] <= range.GetEnd()) {
      mask[i] = (Mask)1;
    }
  }
}

/**
 *
 */
Float64 CSpectrumSpectralAxis::IntersectMaskAndComputeOverlapFraction(
    const TFloat64Range &lambdaRange, const CMask &omask) const {

  Int32 selfRate = 0;
  Int32 otherRate = 0;

  Int32 nSamples = m_Samples.size();
  // weight = Spectrum over lambdarange flag
  for (Int32 i = 0; i < nSamples; i++) {
    // If this sample is somewhere in a valid lambdaRande, tag weight with 1
    if (m_Samples[i] >= lambdaRange.GetBegin() &&
        m_Samples[i] <= lambdaRange.GetEnd()) {
      if (omask[i]) {
        otherRate++;
      }
      selfRate++;
    }
  }

  if (selfRate == 0.0)
    return 0;

  return (Float64)otherRate / (Float64)selfRate;
}

/**
 *
 */
void CSpectrumSpectralAxis::ClampLambdaRange(
    const TFloat64Range &range, TFloat64Range &clampedRange) const {
  TFloat64Range effectiveRange = GetLambdaRange();

  clampedRange = TFloat64Range();

  if (range.GetIsEmpty())
    THROWG(INVALID_WAVELENGTH_RANGE, "lambda range is empty");

  if (effectiveRange.GetIsEmpty())
    THROWG(INVALID_SPECTRUM_WAVELENGTH,
           "spectral axis has less than 2 samples");

  // Clamp lambda start
  Float64 start = range.GetBegin();
  if (start < effectiveRange.GetBegin())
    start = effectiveRange.GetBegin();

  // CLamp lambda end
  Float64 end = range.GetEnd();
  if (end > effectiveRange.GetEnd())
    end = effectiveRange.GetEnd();

  clampedRange = TFloat64Range(start, end);

  if (clampedRange.GetLength() <= 0.0)
    THROWG(INVALID_WAVELENGTH_RANGE,
           "lambda range clamped to spectral axis is empty");
}

/**
 *
 */
TInt32Range CSpectrumSpectralAxis::GetIndexesAtWaveLengthRange(
    const TFloat64Range &waveLengthRange) const {
  TInt32Range r;

  r.SetBegin(GetIndexAtWaveLength(waveLengthRange.GetBegin()));
  r.SetEnd(GetIndexAtWaveLength(waveLengthRange.GetEnd()));

  return r;
}

/**
 *
 */
Int32 CSpectrumSpectralAxis::GetIndexAtWaveLength(Float64 waveLength) const {
  Int32 m;
  Int32 lo = 0;
  Int32 hi = GetSamplesCount() - 1;

  if (waveLength <= m_Samples[lo])
    return lo;

  if (waveLength >= m_Samples[hi])
    return hi;

  auto it = std::lower_bound(m_Samples.begin(), m_Samples.end(), waveLength);

  return (it - m_Samples.begin());
}

/**
 * Check if spectralAxis is well rebinned in log
 */
bool CSpectrumSpectralAxis::CheckLoglambdaSampling() const {
  Float64 logGridStep = log(m_Samples[1] / m_Samples[0]);

  Float64 relativelogGridStepTol =
      1e-1; // only 10% precision (not better to keep input spectra with
            // truncatd decimals)
  Float64 maxAbsRelativeError = 0.0;
  Float64 lbda1 = m_Samples[0];
  for (Int32 t = 1; t < m_Samples.size(); t++) {
    Float64 lbda2 = m_Samples[t];
    Float64 _logGridStep = log(lbda2 / lbda1);

    Float64 relativeErrAbs =
        std::abs((_logGridStep - logGridStep) / logGridStep);
    maxAbsRelativeError = max(relativeErrAbs, maxAbsRelativeError);

    if (relativeErrAbs >
        relativelogGridStepTol) { // return without setting anything
      m_isLogSampled = false;     // unsetting
      Log.LogDebug("   CSpectrumSpectralAxis::CheckLoglambdaSampling: "
                   "Log-regular lambda check FAILED");

      return false;
    }
    lbda1 = lbda2;
  }
  //  recompute log step with more precision:
  m_regularLogSamplingStep =
      log(m_Samples.back() / m_Samples.front()) / (GetSamplesCount() - 1);
  m_isLogSampled = true;
  Log.LogDetail(Formatter()
                << "   CSpectrumSpectralAxis::CheckLoglambdaSampling: max Abs "
                   "Relative Error (log lbda step)= "
                << maxAbsRelativeError);

  return true;
}

/**
 * Brief:
 * In the actual version we consider that a spectral axis can be log sampled
 * while having its values in Angstrom (i.e., non-log)
 *
 */
bool CSpectrumSpectralAxis::IsLogSampled(Float64 logGridstep) const {
  if (!IsLogSampled())
    return false;

  if (std::abs(m_regularLogSamplingStep - logGridstep) > 1E-7) {
    Log.LogDetail("   CSpectrumSpectralAxis::IsLogSampled: Log-regular "
                  "sampling with bad step");
    return false;
  }

  return true;
}

bool CSpectrumSpectralAxis::IsLogSampled() const {
  if (indeterminate(m_isLogSampled))
    CheckLoglambdaSampling();

  return bool(m_isLogSampled);
}

Float64 CSpectrumSpectralAxis::GetlogGridStep() const {
  if (!IsLogSampled()) {
    THROWG(INTERNAL_ERROR, " axis is not logsampled");
  }

  return m_regularLogSamplingStep;
}

// still TODO: check end-to-end redshift coverage
TFloat64List CSpectrumSpectralAxis::GetSubSamplingMask(Int32 ssratio) const {
  return GetSubSamplingMask(ssratio, TInt32Range(0, GetSamplesCount() - 1));
}

TFloat64List CSpectrumSpectralAxis::GetSubSamplingMask(
    Int32 ssratio, TFloat64Range const &lambdarange) const {
  Int32 imin = -1, imax = m_Samples.size();
  lambdarange.getClosedIntervalIndices(m_Samples, imin, imax);
  return GetSubSamplingMask(ssratio, TInt32Range(imin, imax));
}

/*@ssratio stands for sub-samplingRatio*/
TFloat64List
CSpectrumSpectralAxis::GetSubSamplingMask(Int32 ssratio,
                                          const TInt32Range &ilbda) const {
  if (!IsLogSampled()) {
    THROWG(INTERNAL_ERROR, "Cannot subsample spectrum");
  }
  if (ilbda.GetBegin() < 0)
    THROWG(INTERNAL_ERROR, "range's lower bound < 0");
  if (ilbda.GetEnd() > m_Samples.size() - 1)
    THROWG(INTERNAL_ERROR, "range's upper bound > samples size");

  Int32 s = GetSamplesCount();
  if (ssratio == 1)
    return TFloat64List(s, 1.);
  TFloat64List mask(s, 0.);
  for (Int32 i = ilbda.GetEnd(); i >= ilbda.GetBegin();
       i -= ssratio) { // ensure that z[0] remains the same
    mask[i] = 1;
  }
  return mask;
}

/**
 * Brief:
 * check that division of two floating values gives an int, and return modulo
 value
 *  Int32 logstep_int = Int32(logstep*1E12);
    Int32 rebinlogstep_int = Int32(rebinlogstep*1E12);
    auto lambda_redshift_modulo = logstep_int %rebinlogstep_int;
    auto modulo_2 = logstep_int -
 trunc(logstep_int/rebinlogstep_int)*rebinlogstep_int;
*/
Int32 CSpectrumSpectralAxis::GetLogSamplingIntegerRatio(Float64 logstep,
                                                        Float64 &modulo) const {
  if (!IsLogSampled()) {
    THROWG(INTERNAL_ERROR,
           "axis is not logsampled, thus cannot get integer ratio");
  }

  Int32 ratio = std::round(logstep / m_regularLogSamplingStep);
  modulo = logstep - ratio * m_regularLogSamplingStep;
  return ratio;
}

void CSpectrumSpectralAxis::RecomputePreciseLoglambda() {
  if (!IsLogSampled()) {
    THROWG(INTERNAL_ERROR, "axis is not logsampled");
  }
  TFloat64Range lrange = GetLambdaRange();
  TFloat64List new_Samples =
      lrange.SpreadOverLogEpsilon(m_regularLogSamplingStep);

  // gain one more decimal
  Int32 bs;
  const Int32 nm1 = m_Samples.size() - 1;
  bs = min(100, nm1 / 10);
  if (bs < 30) {
    THROWG(INTERNAL_ERROR, "not enough points to recompute logLambda");
  }

  Float64 bias_start = 0., bias_end = 0.;
  // take the mean value (assuming rounding to even),
  //  should take the max value if truncation
  for (Int32 k = 0; k < bs; k++) {
    bias_start += log(new_Samples[k] / m_Samples[k]);
    bias_end += log(new_Samples[nm1 - k] / m_Samples[nm1 - k]);
  }
  bias_start /= bs;
  bias_end /= bs;
  Float64 lstart = log(lrange.GetBegin());
  Float64 lend = log(lrange.GetEnd());
  Float64 new_lstart = lstart - bias_start;
  Float64 new_lend = lend - bias_end;
  TFloat64Range new_lrange;
  new_lrange.Set(exp(new_lstart), exp(new_lend));

  Float64 new_regularLogSamplingStep = (new_lend - new_lstart) / nm1;
  m_Samples = new_lrange.SpreadOverLogEpsilon(new_regularLogSamplingStep);
  m_regularLogSamplingStep = new_regularLogSamplingStep;
}

/**
 * @brief Check if spectral axis is sorted in the increasing order
 * empty or constant vectors are considered as non-sorted
 */
bool CSpectrumSpectralAxis::isSorted() const {
  if (!indeterminate(m_isSorted))
    return bool(m_isSorted);
  // first time checking
  m_isSorted = false;
  if (std::is_sorted(std::begin(m_Samples), std::end(m_Samples))) {
    auto it = std::adjacent_find(m_Samples.begin(),
                                 m_Samples.end()); // look for duplicates
    if (it == m_Samples.end())                     // no duplicates
      m_isSorted = true;
  }
  return bool(m_isSorted);
}
void CSpectrumSpectralAxis::resetAxisProperties() {
  // reset states since m_Samples is going to change
  m_isSorted = indeterminate;
  m_isLogSampled = indeterminate;
}
