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
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/continuum/irregularsamplingmedian.h"
#include "RedshiftLibrary/log/log.h"
#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

using namespace NSEpic;
using namespace std;

CSpectrum::CSpectrum()
    : m_estimationMethod(""), m_medianWindowSize(-1),
      m_medianEvenReflection(true), m_Name("") {}

CSpectrum::CSpectrum(const std::string &name) : m_Name(name){};

CSpectrum::CSpectrum(const CSpectrum &other, const TFloat64List &mask)
    : m_estimationMethod(other.m_estimationMethod),
      m_medianWindowSize(other.m_medianWindowSize),
      m_medianEvenReflection(other.m_medianEvenReflection),
      m_Name(other.m_Name), m_spcType(other.m_spcType), m_LSF(other.m_LSF),
      alreadyRemoved(other.alreadyRemoved),
      m_SpectralAxis(Int32(0), other.m_SpectralAxis.IsInLogScale()) {
  const CSpectrumNoiseAxis &otherRawError = other.m_RawFluxAxis.GetError(),
                           &otherContinuumError =
                               other.m_ContinuumFluxAxis.GetError(),
                           &otherWithoutContinuumError =
                               other.m_WithoutContinuumFluxAxis.GetError();

  CSpectrumNoiseAxis &RawError = m_RawFluxAxis.GetError(),
                     &ContinuumError = m_ContinuumFluxAxis.GetError(),
                     &WithoutContinuumError =
                         m_WithoutContinuumFluxAxis.GetError();

  other.m_SpectralAxis.MaskAxis(mask, m_SpectralAxis);
  other.m_RawFluxAxis.MaskAxis(mask, m_RawFluxAxis);

  if (!otherRawError.isEmpty())
    otherRawError.MaskAxis(mask, RawError);

  if (other.alreadyRemoved) {
    other.m_ContinuumFluxAxis.MaskAxis(mask, m_ContinuumFluxAxis);
    other.m_WithoutContinuumFluxAxis.MaskAxis(mask, m_WithoutContinuumFluxAxis);

    if (!otherContinuumError.isEmpty())
      otherContinuumError.MaskAxis(mask, ContinuumError);

    if (!otherWithoutContinuumError.isEmpty())
      otherWithoutContinuumError.MaskAxis(mask, WithoutContinuumError);
  }
}

CSpectrum::CSpectrum(CSpectrumSpectralAxis spectralAxis,
                     CSpectrumFluxAxis fluxAxis)
    : CSpectrum(std::move(spectralAxis), std::move(fluxAxis), nullptr) {}

CSpectrum::CSpectrum(CSpectrumSpectralAxis spectralAxis,
                     CSpectrumFluxAxis fluxAxis,
                     const std::shared_ptr<const CLSF> &lsf)
    : m_SpectralAxis(std::move(spectralAxis)),
      m_RawFluxAxis(std::move(fluxAxis)), m_estimationMethod(""),
      m_medianWindowSize(-1), m_medianEvenReflection(true), m_Name(""),
      m_LSF(lsf) {
  if (!IsValid()) {
    THROWG(INVALID_SPECTRUM,
           "Invalid spectrum with empty axes, non-matching size "
           "or unsorted spectral axis");
  }
}

// copy constructor
//  copy everything exept fullpath, rebined buffer (m_FineGridInterpolated,
//  m_pfgFlux) and const members (m_dLambdaFineGrid &m_method2baseline)
CSpectrum::CSpectrum(const CSpectrum &other)
    : m_estimationMethod(other.m_estimationMethod),
      m_medianWindowSize(other.m_medianWindowSize),
      m_medianEvenReflection(other.m_medianEvenReflection),
      m_SpectralAxis(other.m_SpectralAxis), m_RawFluxAxis(other.m_RawFluxAxis),
      m_ContinuumFluxAxis(other.m_ContinuumFluxAxis),
      m_WithoutContinuumFluxAxis(other.m_WithoutContinuumFluxAxis),
      m_spcType(other.m_spcType), m_LSF(other.m_LSF), m_Name(other.m_Name),
      alreadyRemoved(other.alreadyRemoved) {
  if (!IsValid()) {
    THROWG(INVALID_SPECTRUM,
           "Invalid spectrum with empty axes, non-matching size "
           "or unsorted spectral axis");
  }
}

CSpectrum::CSpectrum(CSpectrum &&other)
    : m_estimationMethod(std::move(other.m_estimationMethod)),
      m_medianWindowSize(other.m_medianWindowSize),
      m_medianEvenReflection(other.m_medianEvenReflection),
      m_SpectralAxis(std::move(other.m_SpectralAxis)),
      m_RawFluxAxis(std::move(other.m_RawFluxAxis)),
      m_ContinuumFluxAxis(std::move(other.m_ContinuumFluxAxis)),
      m_WithoutContinuumFluxAxis(std::move(other.m_WithoutContinuumFluxAxis)),
      m_spcType(other.m_spcType), m_LSF(std::move(other.m_LSF)),
      m_Name(std::move(other.m_Name)), alreadyRemoved(other.alreadyRemoved) {
  if (!IsValid()) {
    THROWG(INVALID_SPECTRUM,
           "Invalid spectrum with empty axes, non-matching size "
           "or unsorted spectral axis");
  }
}

CSpectrum::~CSpectrum() {}

// copy assignment operator
//  same logic as copy constructor
CSpectrum &CSpectrum::operator=(const CSpectrum &other) {
  m_SpectralAxis = other.m_SpectralAxis;
  m_RawFluxAxis = other.m_RawFluxAxis;
  m_ContinuumFluxAxis = other.m_ContinuumFluxAxis;
  m_WithoutContinuumFluxAxis = other.m_WithoutContinuumFluxAxis;
  m_spcType = other.m_spcType;
  SetType(m_spcType);

  m_LSF = other.m_LSF;

  m_estimationMethod = other.m_estimationMethod;
  m_medianWindowSize = other.m_medianWindowSize;
  m_medianEvenReflection = other.m_medianEvenReflection;
  m_Name = other.m_Name;
  alreadyRemoved = other.alreadyRemoved;
  return *this;
}

CSpectrum &CSpectrum::operator=(CSpectrum &&other) {
  m_SpectralAxis = std::move(other.m_SpectralAxis);
  m_RawFluxAxis = std::move(other.m_RawFluxAxis);
  m_ContinuumFluxAxis = std::move(other.m_ContinuumFluxAxis);
  m_WithoutContinuumFluxAxis = std::move(other.m_WithoutContinuumFluxAxis);
  m_spcType = other.m_spcType;
  SetType(m_spcType);

  m_LSF = std::move(other.m_LSF);

  m_estimationMethod = std::move(other.m_estimationMethod);
  m_medianWindowSize = other.m_medianWindowSize;
  m_medianEvenReflection = other.m_medianEvenReflection;
  m_Name = std::move(other.m_Name);
  alreadyRemoved = other.alreadyRemoved;

  return *this;
}

void CSpectrum::SetSpectralAxis(const CSpectrumSpectralAxis &spectralaxis) {
  if (spectralaxis.GetSamplesCount() != GetFluxAxis().GetSamplesCount()) {
    THROWG(INVALID_SPECTRUM,
           "CSpectrum::SetSpectralAxis: new spectral axis has "
           "not the same size than flux axis");
  }
  m_SpectralAxis = spectralaxis;
}

void CSpectrum::SetSpectralAxis(CSpectrumSpectralAxis &&spectralaxis) {
  if (spectralaxis.GetSamplesCount() != GetFluxAxis().GetSamplesCount()) {
    THROWG(INVALID_SPECTRUM,
           "CSpectrum::SetSpectralAxis: new spectral axis has "
           "not the same size than flux axis");
  }
  m_SpectralAxis = std::move(spectralaxis);
}

void CSpectrum::SetFluxAxis(const CSpectrumFluxAxis &fluxaxis) {
  if (fluxaxis.GetSamplesCount() != m_SpectralAxis.GetSamplesCount()) {
    THROWG(INVALID_SPECTRUM,
           "CSpectrum::SetFluxAxis: new flux axis has not the "
           "same size than spectral axis");
  }

  ResetContinuum();
  SetType(EType::nType_raw);
  ClearFineGrid();
  m_RawFluxAxis = fluxaxis;
}

void CSpectrum::SetFluxAxis(CSpectrumFluxAxis &&fluxaxis) {
  if (fluxaxis.GetSamplesCount() != m_SpectralAxis.GetSamplesCount() ||
      m_SpectralAxis.isEmpty()) {
    THROWG(INVALID_SPECTRUM,
           "CSpectrum::SetFluxAxis: new flux axis has not the "
           "same size than spectral axis");
  }

  ResetContinuum();
  SetType(EType::nType_raw);
  ClearFineGrid();
  m_RawFluxAxis = std::move(fluxaxis);
}

void CSpectrum::SetSpectralAndFluxAxes(CSpectrumSpectralAxis spcaxis,
                                       CSpectrumFluxAxis fluxaxis) {
  if (fluxaxis.GetSamplesCount() != spcaxis.GetSamplesCount() ||
      spcaxis.isEmpty()) {
    THROWG(INVALID_SPECTRUM, "CSpectrum::SetSpectralAndFluxAxes: new flux axis "
                             "has not the same size than new spectral axis");
  }

  ResetContinuum();
  SetType(EType::nType_raw);
  ClearFineGrid();

  m_SpectralAxis = std::move(spcaxis);
  SetFluxAxis(std::move(fluxaxis));
}

void CSpectrum::InitSpectrum(CParameterStore &parameterStore) {
  if (!IsValid())
    THROWG(INVALID_SPECTRUM,
           "Invalid spectrum with empty axes, non-matching size "
           "or unsorted spectral axis");

  Float64 smoothWidth = parameterStore.Get<Float64>("smoothWidth");
  std::string medianRemovalMethod =
      parameterStore.Get<std::string>("continuumRemoval.method");
  Float64 medianKernelWidth =
      parameterStore.Get<Float64>("continuumRemoval.medianKernelWidth");
  bool medianEvenReflection =
      parameterStore.Get<bool>("continuumRemoval.medianEvenReflection");
  SetType(EType::nType_raw);
  if (smoothWidth > 0) {
    m_RawFluxAxis.ApplyMeanSmooth(smoothWidth);
  }
  ResetContinuum();
  SetContinuumEstimationMethod(medianRemovalMethod);
  SetMedianWinsize(medianKernelWidth);
  SetMedianEvenReflection(medianEvenReflection);
}

/**
 * below should be calculated in the case of precomputedfinegrid
 */
bool CSpectrum::RebinFineGrid() const {
  // Precalculate a fine grid template to be used for the 'closest value' rebin
  // method
  Int32 n = GetSampleCount();
  if (!n)
    return false;

  Float64 lmin = m_SpectralAxis[0]; // template wavelength never starts at 0
  Float64 lmax = m_SpectralAxis[n - 1];
  Int32 nTgt = (lmax - lmin) / m_dLambdaFineGrid + 2.0 / m_dLambdaFineGrid;

  m_pfgFlux.resize(nTgt);

  const TAxisSampleList &Ysrc = GetFluxAxis().GetSamplesVector();
  const TAxisSampleList &Xsrc = m_SpectralAxis.GetSamplesVector();

  // Initialize and allocate the gsl objects
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline_init(spline, Xsrc.data(), Ysrc.data(), n);
  gsl_interp_accel *accelerator = gsl_interp_accel_alloc();

  Int32 k = 0;
  Float64 x = 0.0;
  for (k = 0; k < nTgt; k++) {
    x = lmin + k * m_dLambdaFineGrid;
    if (x < m_SpectralAxis[0] || x > m_SpectralAxis[n - 1]) {
      m_pfgFlux[k] = 0.0;
    } else {
      m_pfgFlux[k] = gsl_spline_eval(spline, x, accelerator);
    }
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(accelerator);
  m_FineGridInterpolated = true;
  return true;
}

void CSpectrum::ResetContinuum() const {
  alreadyRemoved = false;

  m_ContinuumFluxAxis.clear();
  m_WithoutContinuumFluxAxis.clear();
}

bool CSpectrum::RemoveContinuum(CContinuum &remover) const {
  ResetContinuum();

  return remover.RemoveContinuum(*this, m_WithoutContinuumFluxAxis);
}

/**
 * Estimate the continuum component
 */
void CSpectrum::EstimateContinuum() const {
  Log.LogDetail("Continuum estimation on input spectrum: using %s",
                m_estimationMethod.c_str());

  if (m_estimationMethod == "IrregularSamplingMedian") {
    CContinuumIrregularSamplingMedian continuum;
    continuum.SetMedianKernelWidth(m_medianWindowSize);
    continuum.SetMeanKernelWidth(m_medianWindowSize);
    continuum.SetMedianEvenReflection(m_medianEvenReflection);
    RemoveContinuum(continuum);
    Log.LogDetail("Continuum estimation - medianKernelWidth = %.2f",
                  m_medianWindowSize);
  } else if (m_estimationMethod == "raw") {
    Int32 nbSamples = this->GetSampleCount();
    m_WithoutContinuumFluxAxis.SetSize(nbSamples);
    for (Int32 k = 0; k < nbSamples; k++) {
      m_WithoutContinuumFluxAxis[k] = 0.0;
    }
  } else if (m_estimationMethod == "zero") {
    m_WithoutContinuumFluxAxis = m_RawFluxAxis;
  } else if (m_estimationMethod == "manual") {
    m_WithoutContinuumFluxAxis = m_RawFluxAxis;
    m_WithoutContinuumFluxAxis.Subtract(m_ContinuumFluxAxis);
  } else {
    THROWG(INTERNAL_ERROR, "Unknown Estimation method");
  }

  Log.LogDetail("===============================================");

  // Fill m_ContinuumFluxAxis
  if (m_estimationMethod != "manual") {
    m_ContinuumFluxAxis = m_RawFluxAxis;
    m_ContinuumFluxAxis.Subtract(m_WithoutContinuumFluxAxis);
  }

  alreadyRemoved = true;
}

/**
 * Invert the flux axis
 */
bool CSpectrum::InvertFlux() {
  m_RawFluxAxis.Invert();
  if (alreadyRemoved) {
    m_ContinuumFluxAxis.Invert();
    m_WithoutContinuumFluxAxis.Invert();
  }
  return true;
}

/**
 * Convert the spectral axis to a neperian logarithm scale
 */
bool CSpectrum::ConvertToLogScale() {
  return m_SpectralAxis.ConvertToLogScale();
}

/**
 * Convert the spectral axis to a linear scale
 */
bool CSpectrum::ConvertToLinearScale() {
  return m_SpectralAxis.ConvertToLinearScale();
}

Float64 CSpectrum::GetResolution() const {
  return m_SpectralAxis.GetResolution();
}

Float64 CSpectrum::GetMeanResolution() const {
  return m_SpectralAxis.GetMeanResolution();
}

/**
 * Return the lambda range of the entire spectrum.
 * Range is always expressed in linear scale NOT in log scale even if the
 * underlying spcetrum is in log scale
 */
TLambdaRange CSpectrum::GetLambdaRange() const {
  return m_SpectralAxis.GetLambdaRange();
}

bool CSpectrum::GetMeanAndStdFluxInRange(TFloat64Range wlRange, Float64 &mean,
                                         Float64 &std) const {
  // wlrange should be totally included in the spectrum lambdarange
  if (wlRange.GetBegin() < m_SpectralAxis.GetLambdaRange().GetBegin()) {
    return false;
  }
  if (wlRange.GetEnd() > m_SpectralAxis.GetLambdaRange().GetEnd()) {
    return false;
  }

  CMask mask;
  m_SpectralAxis.GetMask(wlRange, mask);
  const CSpectrumNoiseAxis &error = GetFluxAxis().GetError();
  Float64 _Mean = 0.0;
  Float64 _SDev = 0.0;
  GetFluxAxis().ComputeMeanAndSDev(mask, _Mean, _SDev);

  mean = _Mean;
  std = _SDev;
  return true;
}

bool CSpectrum::GetLinearRegInRange(TFloat64Range wlRange, Float64 &a,
                                    Float64 &b) const {
  // wlrange should be totally included in the spectrum lambdarange
  if (wlRange.GetBegin() < m_SpectralAxis.GetLambdaRange().GetBegin()) {
    return false;
  }
  if (wlRange.GetEnd() > m_SpectralAxis.GetLambdaRange().GetEnd()) {
    return false;
  }

  const CSpectrumNoiseAxis &error = GetFluxAxis().GetError();

  TInt32Range iRange = m_SpectralAxis.GetIndexesAtWaveLengthRange(wlRange);
  Int32 n = iRange.GetEnd() - iRange.GetBegin() + 1;
  Float64 *x = (Float64 *)malloc(n * sizeof(Float64));
  Float64 *y = (Float64 *)malloc(n * sizeof(Float64));
  Float64 *w = (Float64 *)malloc(n * sizeof(Float64));

  for (Int32 k = 0; k < n; k++) {
    Int32 ik = k + iRange.GetBegin();
    w[k] = 1.0 / (error[ik] * error[ik]);
    x[k] = m_SpectralAxis[ik];
    y[k] = GetFluxAxis().GetSamples()[ik];
  }

  double c0, c1, cov00, cov01, cov11, chisq;
  gsl_fit_wlinear(x, 1, w, 1, y, 1, n, &c0, &c1, &cov00, &cov01, &cov11,
                  &chisq);

  a = c1;
  b = c0;
  // todo: use chisq (weighted sum of squares of the residuals) for fitting
  // quality,
  free(x);
  free(y);
  free(w);
  return true;
}

const std::string &CSpectrum::GetName() const { return m_Name; }

void CSpectrum::SetName(std::string name) { m_Name = std::move(name); }

const CSpectrum::EType CSpectrum::GetType() const { return m_spcType; }

void CSpectrum::SetType(const CSpectrum::EType type) const {
  if (m_spcType != type) {
    ClearFineGrid();
    m_spcType = type;
  }
}

const bool CSpectrum::checkFlux(Float64 flux, Int32 index) const {
  // Log.LogDebug("    CSpectrum::checkFlux - Found flux value (=%e) at
  // index=%d", flux, index);
  bool validValue = true;
  if (std::isnan(flux)) {
    validValue = false;
    // Log.LogDebug("    CSpectrum::checkFlux - Found nan flux value (=%e) at
    // index=%d", flux, index);
  }
  if (std::isinf(flux)) {
    validValue = false;
    // Log.LogDebug("    CSpectrum::checkFlux - Found inf flux value (=%e) at
    // index=%d", flux, index);
  }
  if (flux != flux) {
    validValue = false;
    // Log.LogDebug("    CSpectrum::checkFlux - Found invalid flux value (=%e)
    // at index=%d", flux, index);
  }
  return validValue;
}

const bool CSpectrum::checkNoise(Float64 error, Int32 index) const {
  // Log.LogDebug("    CSpectrum::checkNoise - Found noise value (=%e) at
  // index=%d", error, index);
  bool validValue = true;
  if (error < DBL_MIN) {
    // check if noise is below minimum normalized positive value of double
    validValue = false;
    // Log.LogDebug("    CSpectrum::checkNoise - Found subnormal noise value
    // (=%e) at index=%d", error, index);
  }
  if (std::isnan(error)) {
    validValue = false;
    // Log.LogDebug("    CSpectrum::checkNoise - Found nan noise value (=%e) at
    // index=%d", error, index);
  }
  if (std::isinf(error)) {
    validValue = false;
    // Log.LogDebug("    CSpectrum::checkNoise - Found inf noise value (=%e) at
    // index=%d", error, index);
  }
  if (error != error) {
    validValue = false;
    // Log.LogDebug("    CSpectrum::checkNoise - Found invalid noise value (=%e)
    // at index=%d", error, index);
  }
  return validValue;
}

const bool CSpectrum::IsFluxValid(Float64 LambdaMin, Float64 LambdaMax) const {
  bool allzero = true;
  bool invalidValue = false;
  Int32 nInvalid = 0;

  if (!IsValid())
    THROWG(INVALID_SPECTRUM,
           "Invalid spectrum with empty axes, non-matching size "
           "or unsorted spectral axis");

  const Float64 *flux = GetFluxAxis().GetSamples();
  if (LambdaMin < m_SpectralAxis[0] ||
      LambdaMax > m_SpectralAxis[m_SpectralAxis.GetSamplesCount() - 1]) {
    return false;
  } else {
    Int32 iMin = m_SpectralAxis.GetIndexAtWaveLength(LambdaMin);
    Int32 iMax = m_SpectralAxis.GetIndexAtWaveLength(LambdaMax);
    Log.LogDetail("CSpectrum::IsFluxValid - checking on the configured "
                  "lambdarange = (%f, %f)",
                  LambdaMin, LambdaMax);
    Log.LogDetail("CSpectrum::IsFluxValid - checking on the true observed "
                  "spectral axis lambdarange = (%f, %f)",
                  m_SpectralAxis[iMin], m_SpectralAxis[iMax]);
    for (Int32 i = iMin; i < iMax; i++) {
      // check flux
      bool validSample = checkFlux(flux[i], i);

      if (!validSample) {
        invalidValue = true;
        nInvalid++;
      }

      // all zero check
      if (flux[i] != 0.0) {
        allzero = false;
        // Log.LogDebug("    CSpectrum::IsFluxValid - Found non zero and valid
        // flux value (=%e) at index=%d", i, flux[i]);
      }
    }
    bool valid = !invalidValue && !allzero;
    if (nInvalid > 0) {
      Log.LogDetail(
          "    CSpectrum::IsFluxValid - Found %d invalid flux samples",
          nInvalid);
    }
    return valid;
  }
}

const bool CSpectrum::IsNoiseValid(Float64 LambdaMin, Float64 LambdaMax) const {
  bool valid = true;
  Int32 nInvalid = 0;

  if (IsNoiseEmpty())
    return false;

  const TFloat64List &error = GetFluxAxis().GetError().GetSamplesVector();
  if (LambdaMin < m_SpectralAxis[0] ||
      LambdaMax > m_SpectralAxis[m_SpectralAxis.GetSamplesCount() - 1]) {
    return false;
  } else {
    Int32 iMin = m_SpectralAxis.GetIndexAtWaveLength(LambdaMin);
    Int32 iMax = m_SpectralAxis.GetIndexAtWaveLength(LambdaMax);
    Log.LogDetail("CSpectrum::IsNoiseValid - checking on the configured "
                  "lambdarange = (%f, %f)",
                  LambdaMin, LambdaMax);
    Log.LogDetail("CSpectrum::IsNoiseValid - checking on the true observed "
                  "spectral axis lambdarange = (%f, %f)",
                  m_SpectralAxis[iMin], m_SpectralAxis[iMax]);
    for (Int32 i = iMin; i < iMax; i++) {
      // check noise
      bool validSample =
          checkNoise(error[i], i); // checkNoise can be moved to CSpectrumNoise

      if (!validSample) {
        valid = false;
        nInvalid++;
      }
    }
    if (nInvalid > 0) {
      Log.LogDetail(
          "    CSpectrum::IsNoiseValid - Found %d invalid noise samples",
          nInvalid);
    }
    return valid;
  }
}

bool CSpectrum::correctSpectrum(Float64 LambdaMin, Float64 LambdaMax,
                                Float64 coeffCorr) {
  bool corrected = false;
  Int32 nCorrected = 0;
  if (!IsValid())
    THROWG(INVALID_SPECTRUM,
           "Invalid spectrum with empty axes, non-matching size "
           "or unsorted spectral axis");

  CSpectrumFluxAxis fluxaxis = std::move(GetFluxAxis_());
  TFloat64List &error = fluxaxis.GetError().GetSamplesVector();
  TFloat64List &flux = fluxaxis.GetSamplesVector();

  Int32 iMin = m_SpectralAxis.GetIndexAtWaveLength(LambdaMin);
  Int32 iMax = m_SpectralAxis.GetIndexAtWaveLength(LambdaMax);
  // Log.LogDebug("    CSpectrum::correctSpectrum - debug - iMin=%d and wmin=%f,
  // iMax=%d and wmax=%f", iMin, m_SpectralAxis[iMin], iMax,
  // m_SpectralAxis[iMax]);

  Float64 maxNoise = -DBL_MAX;
  Float64 minFlux = DBL_MAX;
  for (Int32 i = iMin; i < iMax; i++) {
    // check noise & flux
    if (!checkNoise(error[i], i) || !checkFlux(flux[i], i))
      continue;

    maxNoise = std::max(maxNoise, error[i]);
    minFlux = std::min(std::abs(flux[i]), std::abs(minFlux));
  }
  if (minFlux == DBL_MAX) {
    Log.LogError("    CSpectrum::correctSpectrum - unable to set minFlux value "
                 "(=%e). Setting it to 0.",
                 minFlux);
    minFlux = 0.0;
  }
  if (maxNoise == -DBL_MAX)
    THROWG(INTERNAL_ERROR, "Unable to set maxNoise value");

  for (Int32 i = iMin; i < iMax; i++) {
    // check noise & flux
    bool validSample = checkNoise(error[i], i) && checkFlux(flux[i], i);

    if (validSample)
      continue;
    error[i] = maxNoise * coeffCorr;
    flux[i] = minFlux / coeffCorr;
    corrected = true;
    nCorrected++;
  }

  SetFluxAxis(std::move(fluxaxis));

  if (corrected) {
    Log.LogInfo("    CSpectrum::correctSpectrum - Corrected %d invalid samples "
                "with coeff (=%f), minFlux=%e, maxNoise=%e",
                nCorrected, coeffCorr, minFlux, maxNoise);
    ResetContinuum();
  }

  return corrected;
}

const std::string &CSpectrum::GetFullPath() const { return m_FullPath; }

const Float64 CSpectrum::GetMedianWinsize() const { return m_medianWindowSize; }

const bool CSpectrum::GetMedianEvenReflection() const {
  return m_medianEvenReflection;
}

const std::string &CSpectrum::GetContinuumEstimationMethod() const {
  return m_estimationMethod;
}

void CSpectrum::SetFullPath(const char *nameP) { m_FullPath = nameP; }

void CSpectrum::SetMedianWinsize(Float64 winsize) {
  if (m_medianWindowSize != winsize &&
      m_estimationMethod == "IrregularSamplingMedian") {
    ResetContinuum();
  }
  m_medianWindowSize = winsize;
}

void CSpectrum::SetMedianEvenReflection(bool medianEvenReflection) {
  if (m_medianEvenReflection != medianEvenReflection &&
      m_estimationMethod == "IrregularSamplingMedian") {
    ResetContinuum();
  }
  m_medianEvenReflection = medianEvenReflection;
}

void CSpectrum::SetContinuumEstimationMethod(std::string method) const {
  if (m_estimationMethod != method) {
    m_estimationMethod = method;
    ResetContinuum();
  }
}

/*
 *  force manual setting of the continuum
 */
void CSpectrum::SetContinuumEstimationMethod(
    const CSpectrumFluxAxis &ContinuumFluxAxis) {
  m_estimationMethod = "manual";
  ResetContinuum();

  if (ContinuumFluxAxis.GetSamplesCount() != GetSampleCount()) {
    THROWG(INTERNAL_ERROR,
           "Cannot set continuum flux: non-matching size with spectral axis");
  }

  m_ContinuumFluxAxis = ContinuumFluxAxis;
}

void CSpectrum::ClearFineGrid() const {
  m_FineGridInterpolated = false;
  m_pfgFlux.clear();
}

///
/// * This rebin method targets processing speed:
/// - it uses already allocated rebinedFluxAxis, rebinedSpectralAxis and
/// rebinedMask
/// - opt_interp = 'lin' : linear interpolation is performed by default
/// - opt_interp = 'precomputedfinegrid' : nearest grid point interpolation is
/// performed using m_pfgFlux which is the precomputed fine grid
/// - opt_interp = 'spline' : GSL/spline interpolation is performed (TODO - not
/// tested)
/// - opt_interp = 'ngp' : nearest grid point is performed (TODO - not tested)
/**
 * targetSpectralAxis should be expressed in same frame as source SpetralAxis
 */
bool CSpectrum::Rebin(const TFloat64Range &range,
                      const CSpectrumSpectralAxis &targetSpectralAxis,
                      CSpectrum &rebinedSpectrum, CMask &rebinedMask,
                      const std::string &opt_interp,
                      const std::string &opt_error_interp) const {

  if (!IsValid())
    THROWG(INVALID_SPECTRUM,
           "Invalid spectrum with empty axes, non-matching size "
           "or unsorted spectral axis");

  Int32 s = targetSpectralAxis.GetSamplesCount();
  TFloat64Range logIntersectedLambdaRange(log(range.GetBegin()),
                                          log(range.GetEnd()));
  TFloat64Range currentRange = logIntersectedLambdaRange;
  if (m_SpectralAxis.IsInLinearScale() !=
      targetSpectralAxis.IsInLinearScale()) {
    Log.LogError("Problem spectral axis and target spectral axis are not in "
                 "the same scale\n");
    return false;
  }
  if (m_SpectralAxis.IsInLinearScale()) {
    currentRange = range;
  }
  // find start/end indexs for both axes
  if (m_SpectralAxis[0] > currentRange.GetBegin() ||
      m_SpectralAxis[m_SpectralAxis.GetSamplesCount() - 1] <
          currentRange.GetEnd()) {
    THROWG(INTERNAL_ERROR, "input spectral range is not "
                           "included in spectral axis");
  }

  if (opt_interp == "precomputedfinegrid" && m_FineGridInterpolated == false) {
    RebinFineGrid();
  }

  if (m_pfgFlux.size() == 0 && opt_interp == "precomputedfinegrid" &&
      m_FineGridInterpolated == true) {
    Log.LogError("Problem buffer couldnt be computed\n");
    return false;
  }

  CSpectrumFluxAxis rebinedFluxAxis = std::move(rebinedSpectrum.m_RawFluxAxis);
  rebinedFluxAxis.SetSize(s); // does not re-allocate if already allocated

  rebinedMask.SetSize(s);

  const TAxisSampleList &Xsrc = m_SpectralAxis.GetSamplesVector();
  const TAxisSampleList &Ysrc = GetFluxAxis().GetSamplesVector();
  const TAxisSampleList &Xtgt = targetSpectralAxis.GetSamplesVector();
  TAxisSampleList &Yrebin = rebinedFluxAxis.GetSamplesVector();
  const TFloat64List &Error = GetFluxAxis().GetError().GetSamplesVector();
  TFloat64List &ErrorRebin = rebinedFluxAxis.GetError().GetSamplesVector();

  // Move cursors up to lambda range start
  Int32 j = 0;
  while (j < targetSpectralAxis.GetSamplesCount() &&
         Xtgt[j] < currentRange.GetBegin()) {
    rebinedMask[j] = 0;
    Yrebin[j] = 0.0;
    if (opt_error_interp == "rebin" || opt_error_interp == "rebinVariance")
      ErrorRebin[j] = INFINITY;
    j++;
  }

  if (opt_interp == "lin") {
    // Default linear interp.
    Int32 k = 0;
    // For each sample in the valid lambda range interval.
    while (k <= m_SpectralAxis.GetSamplesCount() - 1 &&
           Xsrc[k] <= currentRange.GetEnd()) {
      // For each sample in the target spectrum that are in between two
      // continous source sample
      while (j < targetSpectralAxis.GetSamplesCount() &&
             Xtgt[j] <= Xsrc[k + 1]) {
        // perform linear interpolation of the flux
        Float64 xSrcStep = (Xsrc[k + 1] - Xsrc[k]);
        Float64 t = (Xtgt[j] - Xsrc[k]) / xSrcStep;
        Yrebin[j] = Ysrc[k] + (Ysrc[k + 1] - Ysrc[k]) * t;
        rebinedMask[j] = 1;

        if (opt_error_interp != "no") {
          if (opt_error_interp == "rebin")
            ErrorRebin[j] = Error[k] + (Error[k + 1] - Error[k]) * t;
          else if (opt_error_interp == "rebinVariance") {
            ErrorRebin[j] = sqrt(Error[k] * Error[k] * (1 - t) * (1 - t) +
                                 Error[k + 1] * Error[k + 1] * t * t);
            //*
            Float64 xDestStep = NAN;
            Float64 xStepCompensation = 1.;
            if (j < targetSpectralAxis.GetSamplesCount() - 1) {
              xDestStep = Xtgt[j + 1] - Xtgt[j];
              xStepCompensation = xSrcStep / xDestStep;
            } else if (j > 0) {
              xDestStep = Xtgt[j] - Xtgt[j - 1];
              xStepCompensation = xSrcStep / xDestStep;
            }
            ErrorRebin[j] *= sqrt(xStepCompensation);
          }
        }
        j++;
      }

      k++;
    }
  } else if (opt_interp == "precomputedfinegrid") {
    // Precomputed FINE GRID nearest sample, 20150801
    Int32 k = 0;
    Float64 lmin = m_SpectralAxis[0];
    // For each sample in the target spectrum
    while (j < targetSpectralAxis.GetSamplesCount() &&
           Xtgt[j] <= currentRange.GetEnd()) {
      k = int((Xtgt[j] - lmin) / m_dLambdaFineGrid + 0.5);
      Yrebin[j] = m_pfgFlux[k];
      rebinedMask[j] = 1;

      // note: error rebin not implemented for precomputedfinegrid
      if (opt_error_interp != "no")
        return false;

      j++;
    }
  } else if (opt_interp == "spline") {
    // GSL method spline
    // Initialize and allocate the gsl objects
    Int32 n = m_SpectralAxis.GetSamplesCount();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, n);
    gsl_spline_init(spline, Xsrc.data(), Ysrc.data(), n);
    gsl_interp_accel *accelerator = gsl_interp_accel_alloc();

    // For each sample in the valid lambda range interval.
    while (j < targetSpectralAxis.GetSamplesCount() &&
           Xtgt[j] <= currentRange.GetEnd()) {
      Yrebin[j] = gsl_spline_eval(spline, Xtgt[j], accelerator);
      rebinedMask[j] = 1;

      // note: error rebin not implemented for spline interp
      if (opt_error_interp != "no")
        return false;

      j++;
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(accelerator);

  } else if (opt_interp == "ngp") {
    // nearest sample, lookup
    Int32 k = 0;
    Int32 n = m_SpectralAxis.GetSamplesCount();
    while (j < targetSpectralAxis.GetSamplesCount() &&
           Xtgt[j] <= currentRange.GetEnd()) {
      // k = gsl_interp_bsearch (Xsrc.data(), Xtgt[j], kprev, n);
      k = CIndexing<Float64>::getCloserIndex(Xsrc, Xtgt[j]);
      Float64 xSrcStep = NAN;
      if (k == Xsrc.size() - 1)
        xSrcStep = Xsrc[k] - Xsrc[k - 1];
      else
        xSrcStep = Xsrc[k + 1] - Xsrc[k];

      // closest value
      Yrebin[j] = Ysrc[k];

      if (opt_error_interp != "no") {
        ErrorRebin[j] = Error[k];

        if (opt_error_interp == "rebinVariance") {
          Float64 xDestStep = NAN;
          Float64 xStepCompensation = 1.;
          if (j < targetSpectralAxis.GetSamplesCount() - 1) {
            xDestStep = Xtgt[j + 1] - Xtgt[j];
            xStepCompensation = xSrcStep / xDestStep;
          } else if (j > 0) {
            xDestStep = Xtgt[j] - Xtgt[j - 1];
            xStepCompensation = xSrcStep / xDestStep;
          }
          ErrorRebin[j] *= sqrt(xStepCompensation);
        }
      }

      rebinedMask[j] = 1;
      j++;
    }
  }

  while (j < targetSpectralAxis.GetSamplesCount()) {
    rebinedMask[j] = 0;
    Yrebin[j] = 0.0;
    if (opt_error_interp == "rebin" || opt_error_interp == "rebinVariance")
      ErrorRebin[j] = INFINITY;
    j++;
  }

  rebinedSpectrum.ResetContinuum();
  rebinedSpectrum.SetType(EType::nType_raw);
  rebinedSpectrum.SetSpectralAndFluxAxes(targetSpectralAxis,
                                         std::move(rebinedFluxAxis));

  return true;
}

void CSpectrum::ScaleFluxAxis(Float64 scale) {
  m_RawFluxAxis *= scale;
  if (alreadyRemoved) {
    m_ContinuumFluxAxis *= scale;
    m_WithoutContinuumFluxAxis *= scale;
  }
  if (m_FineGridInterpolated)
    ClearFineGrid();
}

void CSpectrum::ValidateSpectrum(TFloat64Range lambdaRange,
                                 bool enableInputSpcCorrect) {
  if (!IsValid())
    THROWG(INVALID_SPECTRUM, "Invalid spectrum with empty axes or non-matching "
                             "size or unsorted spectral axis");

  TFloat64Range clampedlambdaRange;
  m_SpectralAxis.ClampLambdaRange(lambdaRange, clampedlambdaRange);
  Log.LogInfo("Validate spectrum: (CLambdaRange: %f-%f:%f)",
              clampedlambdaRange.GetBegin(), clampedlambdaRange.GetEnd(),
              GetResolution());

  Float64 lmin = clampedlambdaRange.GetBegin();
  Float64 lmax = clampedlambdaRange.GetEnd();

  // Check if the Spectrum is valid on the lambdarange
  if (enableInputSpcCorrect)
    if (correctSpectrum(lmin, lmax))
      Log.LogInfo(
          "Successfully corrected noise on wavelength range (%.1f ; %.1f)",
          lmin, lmax);

  if (!IsFluxValid(lmin, lmax))
    THROWG(
        INTERNAL_ERROR,
        Formatter() << "Failed to validate spectrum flux on wavelength range ["
                    << lmin << ";" << lmax << "]");
  else
    Log.LogDetail("Successfully validated spectrum flux, on wavelength range "
                  "(%.1f ; %.1f)",
                  lmin, lmax);

  // Check if the noise is valid in the clampedlambdaRange
  if (!IsNoiseValid(lmin, lmax))
    THROWG(INTERNAL_ERROR,
           Formatter() << "Failed to validate noise on wavelength range ["
                       << lmin << ";" << lmax << "]");
  else
    Log.LogDetail(
        "Successfully validated noise on wavelength range (%.1f ; %.1f)", lmin,
        lmax);

  if (!m_LSF)
    return;
  // check if spectrum LSF spectralAxis covers lambdaRange
  if (!m_LSF->checkAvailability(lambdaRange.GetBegin()) ||
      !m_LSF->checkAvailability(lambdaRange.GetEnd())) {
    THROWG(INTERNAL_ERROR, Formatter()
                               << "Failed to validate lsf on wavelength range ["
                               << lmin << ";" << lmax << "]");
  }
}
