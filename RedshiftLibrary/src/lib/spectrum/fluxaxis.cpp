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
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/common/mean.h"
#include "RedshiftLibrary/common/median.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/spectrum/fluxaxis.h"

using namespace NSEpic;
using namespace std;

CSpectrumFluxAxis::CSpectrumFluxAxis(Int32 n, Float64 value)
    : CSpectrumAxis(n, value), m_StdError(n) {}

CSpectrumFluxAxis::CSpectrumFluxAxis(CSpectrumAxis otherFlux,
                                     CSpectrumNoiseAxis otherError)
    : CSpectrumAxis(std::move(otherFlux)), m_StdError(std::move(otherError)) {}

CSpectrumFluxAxis::CSpectrumFluxAxis(const CSpectrumAxis otherFlux)
    : CSpectrumAxis(std::move(otherFlux)), m_StdError(GetSamplesCount()) {}

CSpectrumFluxAxis::CSpectrumFluxAxis(const Float64 *samples, Int32 n)
    : CSpectrumAxis(samples, n), m_StdError(n) {}

CSpectrumFluxAxis::CSpectrumFluxAxis(const TFloat64List &samples)
    : CSpectrumAxis(samples), m_StdError(samples.size()) // default to 1
{}

CSpectrumFluxAxis::CSpectrumFluxAxis(TFloat64List &&samples)
    : CSpectrumAxis(std::move(samples)),
      m_StdError(GetSamplesCount()) // default to 1
{}

CSpectrumFluxAxis::CSpectrumFluxAxis(const Float64 *samples, Int32 n,
                                     const Float64 *error, const Int32 m)
    : CSpectrumAxis(samples, n), m_StdError(error, m) {
  if (m != n) {
    THROWG(ErrorCode::INTERNAL_ERROR,
           "FluxAxis and NoiseAxis sizes do not match");
  }
}

void CSpectrumFluxAxis::setError(CSpectrumNoiseAxis otherError) {
  if (otherError.GetSamplesCount() != m_StdError.GetSamplesCount())
    THROWG(ErrorCode::INTERNAL_ERROR,
           "FluxAxis and NoiseAxis sizes do not match");
  m_StdError = std::move(otherError);
}

void CSpectrumFluxAxis::SetSize(Int32 s) {
  CSpectrumAxis::SetSize(s);
  m_StdError.SetSize(s);
}
void CSpectrumFluxAxis::clear() {
  CSpectrumAxis::clear();
  m_StdError.clear();
}

bool CSpectrumFluxAxis::ApplyMedianSmooth(Int32 kernelHalfWidth) {
  if (kernelHalfWidth == 0)
    return false;

  if (GetSamplesCount() < (kernelHalfWidth) + 1)
    return false;

  TAxisSampleList tmp(m_Samples.size());

  Int32 left = 0;
  Int32 right = 0;
  CMedian<Float64> median;
  Int32 N = GetSamplesCount();

  for (Int32 i = 0; i < N; i++) {
    left = max(0, i - kernelHalfWidth);
    right = min(N - 1, i + kernelHalfWidth);
    tmp[i] =
        median.Find(m_Samples.begin() + left, m_Samples.begin() + right + 1);
  }
  m_Samples = std::move(tmp);

  return true;
}

bool CSpectrumFluxAxis::ApplyMeanSmooth(Int32 kernelHalfWidth) {
  if (kernelHalfWidth == 0)
    return false;

  if (GetSamplesCount() < (kernelHalfWidth) + 1)
    return false;

  TAxisSampleList tmp(m_Samples.size());

  Int32 left = 0;
  Int32 right = 0;
  CMean<Float64> mean;
  Int32 N = GetSamplesCount();

  for (Int32 i = 0; i < N; i++) {
    left = max(0, i - kernelHalfWidth);
    right = min(N - 1, i + kernelHalfWidth);

    tmp[i] = mean.Find(m_Samples.begin() + left, m_Samples.begin() + right + 1);
  }

  m_Samples = std::move(tmp);

  return true;
}

Float64 CSpectrumFluxAxis::computeMaxAbsValue(Int32 imin, Int32 imax) const {

  Float64 maxabsval = std::abs(*std::max_element(
      m_Samples.begin() + imin, m_Samples.begin() + imax + 1,
      [](Float64 a, Float64 b) { return std::abs(a) < std::abs(b); }));

  return maxabsval;
}

bool CSpectrumFluxAxis::ComputeMeanAndSDev(const CMask &mask, Float64 &mean,
                                           Float64 &sdev) const {
  if (mask.GetMasksCount() != GetSamplesCount())
    THROWG(ErrorCode::INTERNAL_ERROR,
           "mask.GetMasksCount() != GetSamplesCount()");

  const CSpectrumNoiseAxis &error = GetError();

  Int32 j;

  Float64 sum = 0.0, sum2 = 0.0, weigthSum = 0.0, weigthSum2 = 0.0, weight;

  for (j = 0; j < GetSamplesCount(); j++) {
#ifdef DEBUG_BUILD
    if (!(mask[j] == 1 || mask[j] == 0))
      THROWG(ErrorCode::INTERNAL_ERROR, "bad mask");
#endif

    weight = 1.0 / (error[j] * error[j]);

    sum += mask[j] * m_Samples[j] * weight;
    sum2 += mask[j] * m_Samples[j] * m_Samples[j] * weight;
    weigthSum += mask[j] * weight;
    weigthSum2 += mask[j] * weight * weight;
  }

  if (weigthSum > 0.0) {
    mean = sum / weigthSum;
    sdev = sqrt((sum2 - mean * mean * weigthSum) /
                (weigthSum - weigthSum2 / weigthSum));
  } else {
    mean = NAN;
    sdev = NAN;
    return false;
  }

  return true;
}

Float64 CSpectrumFluxAxis::ComputeRMSDiff(const CSpectrumFluxAxis &other) {
  Float64 er2 = 0.f;
  Float64 er = 0.f;

  if (other.GetSamplesCount() != GetSamplesCount())
    THROWG(ErrorCode::INTERNAL_ERROR,
           "other.GetSamplesCount() != GetSamplesCount()");

  int n = GetSamplesCount();
  Float64 weight = (Float64)n;
  for (int j = 0; j < n; j++) {
    er2 += (m_Samples[j] - other[j]) * (m_Samples[j] - other[j]) / weight;
  }

  er = sqrt(er2);
  return er;
}

const TBoolList CSpectrumFluxAxis::checkFlux() const {
  TBoolList isValid(m_Samples.size(), true);
  for (std::size_t i = 0; i < m_Samples.size(); i++) {
    if (std::isnan(m_Samples[i]) || std::isinf(m_Samples[i]) ||
        (m_Samples[i] != m_Samples[i]))
      isValid[i] = false;
  }
  return isValid;
}

bool CSpectrumFluxAxis::correctFluxAndNoiseAxis(Int32 iMin, Int32 iMax,
                                                Float64 coeffCorr) {
  bool corrected = false;
  Int32 nCorrected = 0;
  CSpectrumNoiseAxis error = GetError();
  Float64 maxNoise = -DBL_MAX;
  Float64 minFlux = DBL_MAX;

  // check noise & flux
  TBoolList isNoiseValid = error.checkNoise();
  TBoolList isFluxValid = checkFlux();

  for (Int32 i = iMin; i < iMax; i++) {
    if (!isNoiseValid[i] || !isFluxValid[i])
      continue;

    maxNoise = std::max(maxNoise, error[i]);
    minFlux = std::min(std::abs(m_Samples[i]), std::abs(minFlux));
  }
  if (minFlux == DBL_MAX) {
    THROWG(ErrorCode::SPECTRUM_CORRECTION_ERROR,
           Formatter() << "Unable to find a min flux value.");
  }
  if (maxNoise == -DBL_MAX)
    THROWG(ErrorCode::SPECTRUM_CORRECTION_ERROR,
           "Unable to find a max noise value");

  for (Int32 i = iMin; i < iMax; i++) {
    // check noise & flux
    bool validSample = isNoiseValid[i] && isFluxValid[i];

    if (validSample)
      continue;
    error[i] = maxNoise * coeffCorr;
    m_Samples[i] = minFlux / coeffCorr;
    corrected = true;
    nCorrected++;
  }

  if (corrected) {
    Log.LogInfo(Formatter()
                << "    CSpectrumFluxAxis::" << __func__ << "- Corrected "
                << nCorrected << " invalid samples with coeff (=" << coeffCorr
                << "), minFlux=" << minFlux << ", maxNoise=" << maxNoise);
    setError(std::move(error));
  }

  return corrected;
}

bool CSpectrumFluxAxis::Subtract(const CSpectrumFluxAxis &other) {
  if (other.GetSamplesCount() != GetSamplesCount())
    THROWG(ErrorCode::INTERNAL_ERROR,
           "other.GetSamplesCount() != GetSamplesCount()");

  Int32 N = GetSamplesCount();
  for (Int32 i = 0; i < N; i++) {
    m_Samples[i] = m_Samples[i] - other[i];
  }
  return true;
}

bool CSpectrumFluxAxis::Invert() {
  Int32 N = GetSamplesCount();
  for (Int32 i = 0; i < N; i++) {
    m_Samples[i] = -m_Samples[i];
  }
  return true;
}
