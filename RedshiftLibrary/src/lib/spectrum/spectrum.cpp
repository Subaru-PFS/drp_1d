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
#include <map>

#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/continuum/irregularsamplingmedian.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/spectrum/rebin/rebinLinear.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

using namespace NSEpic;
using namespace std;

CSpectrum::CSpectrum()
    : m_estimationMethod(""), m_medianWindowSize(-1),
      m_medianEvenReflection(true), m_Name(""),
      m_rebin(std::unique_ptr<CRebin>(new CRebinLinear(*this))) {}

CSpectrum::CSpectrum(const std::string &name)
    : m_Name(name), m_rebin(std::unique_ptr<CRebin>(new CRebinLinear(*this))){};

CSpectrum::CSpectrum(const CSpectrum &other, const TFloat64List &mask)
    : m_estimationMethod(other.m_estimationMethod),
      m_medianWindowSize(other.m_medianWindowSize),
      m_medianEvenReflection(other.m_medianEvenReflection),
      m_Name(other.m_Name), m_spcType(other.m_spcType), m_LSF(other.m_LSF),
      alreadyRemoved(other.alreadyRemoved), m_SpectralAxis(Int32(0)),
      m_rebin(CRebin::create(other.m_rebin->getType(), *this)),
      m_photData(other.m_photData), m_obsId(other.m_obsId) {
  const CSpectrumNoiseAxis &otherRawError = other.m_RawFluxAxis.GetError(),
                           &otherContinuumError =
                               other.m_ContinuumFluxAxis.GetError(),
                           &otherWithoutContinuumError =
                               other.m_WithoutContinuumFluxAxis.GetError();

  m_SpectralAxis = other.m_SpectralAxis.MaskAxis(mask);
  m_RawFluxAxis = other.m_RawFluxAxis.MaskAxis(mask);

  if (!otherRawError.isEmpty())
    m_RawFluxAxis.setError(otherRawError.MaskAxis(mask));

  if (other.alreadyRemoved) {
    m_ContinuumFluxAxis =
        CSpectrumFluxAxis(other.m_ContinuumFluxAxis.MaskAxis(mask));
    m_WithoutContinuumFluxAxis =
        CSpectrumFluxAxis(other.m_WithoutContinuumFluxAxis.MaskAxis(mask));

    if (!otherContinuumError.isEmpty())
      m_ContinuumFluxAxis.setError(otherContinuumError.MaskAxis(mask));

    if (!otherWithoutContinuumError.isEmpty())
      m_WithoutContinuumFluxAxis.setError(
          otherWithoutContinuumError.MaskAxis(mask));
  }
}

CSpectrum::CSpectrum(CSpectrumSpectralAxis spectralAxis,
                     CSpectrumFluxAxis fluxAxis)
    : CSpectrum(std::move(spectralAxis), std::move(fluxAxis), nullptr) {
  m_rebin = std::unique_ptr<CRebin>(new CRebinLinear(*this));
}

CSpectrum::CSpectrum(CSpectrumSpectralAxis spectralAxis,
                     CSpectrumFluxAxis fluxAxis,
                     const std::shared_ptr<const CLSF> &lsf)
    : m_SpectralAxis(std::move(spectralAxis)),
      m_RawFluxAxis(std::move(fluxAxis)), m_estimationMethod(""),
      m_medianWindowSize(-1), m_medianEvenReflection(true), m_Name(""),
      m_LSF(lsf), m_rebin(std::unique_ptr<CRebin>(new CRebinLinear(*this))) {
  if (!IsValid()) {
    THROWG(ErrorCode::INVALID_SPECTRUM,
           "Invalid spectrum with empty axes, non-matching size "
           "or unsorted spectral axis");
  }
}

// copy constructor
//  copy everything exept fullpath and const members &m_method2baseline
CSpectrum::CSpectrum(const CSpectrum &other)
    : m_estimationMethod(other.m_estimationMethod),
      m_medianWindowSize(other.m_medianWindowSize),
      m_medianEvenReflection(other.m_medianEvenReflection),
      m_SpectralAxis(other.m_SpectralAxis), m_RawFluxAxis(other.m_RawFluxAxis),
      m_ContinuumFluxAxis(other.m_ContinuumFluxAxis),
      m_WithoutContinuumFluxAxis(other.m_WithoutContinuumFluxAxis),
      m_spcType(other.m_spcType), m_LSF(other.m_LSF), m_Name(other.m_Name),
      alreadyRemoved(other.alreadyRemoved),
      m_rebin(CRebin::create(other.m_rebin->getType(), *this)),
      m_photData(other.m_photData), m_obsId(other.m_obsId) {
  if (!IsValid()) {
    THROWG(ErrorCode::INVALID_SPECTRUM,
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
      m_Name(std::move(other.m_Name)), alreadyRemoved(other.alreadyRemoved),
      m_rebin(CRebin::create(other.m_rebin->getType(), *this)),
      m_photData(std::move(other.m_photData)),
      m_obsId(std::move(other.m_obsId)) {
  if (!IsValid()) {
    THROWG(ErrorCode::INVALID_SPECTRUM,
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
  m_rebin = CRebin::create(other.m_rebin->getType(), *this);
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
  m_rebin = CRebin::create(other.m_rebin->getType(), *this);

  return *this;
}

void CSpectrum::SetSpectralAxis(const CSpectrumSpectralAxis &spectralaxis) {
  if (spectralaxis.GetSamplesCount() != GetFluxAxis().GetSamplesCount()) {
    THROWG(ErrorCode::INVALID_SPECTRUM,
           "CSpectrum::SetSpectralAxis: new spectral axis has "
           "not the same size than flux axis");
  }
  m_SpectralAxis = spectralaxis;
}

void CSpectrum::SetSpectralAxis(CSpectrumSpectralAxis &&spectralaxis) {
  if (spectralaxis.GetSamplesCount() != GetFluxAxis().GetSamplesCount()) {
    THROWG(ErrorCode::INVALID_SPECTRUM,
           "CSpectrum::SetSpectralAxis: new spectral axis has "
           "not the same size than flux axis");
  }
  m_SpectralAxis = std::move(spectralaxis);
}

void CSpectrum::SetFluxAxis(const CSpectrumFluxAxis &fluxaxis) {
  if (fluxaxis.GetSamplesCount() != m_SpectralAxis.GetSamplesCount()) {
    THROWG(ErrorCode::INVALID_SPECTRUM,
           "CSpectrum::SetFluxAxis: new flux axis has not the "
           "same size than spectral axis");
  }

  ResetContinuum();
  SetType(EType::nType_raw);
  m_rebin->reset();
  m_RawFluxAxis = fluxaxis;
}

void CSpectrum::SetFluxAxis(CSpectrumFluxAxis &&fluxaxis) {
  if (fluxaxis.GetSamplesCount() != m_SpectralAxis.GetSamplesCount() ||
      m_SpectralAxis.isEmpty()) {
    THROWG(ErrorCode::INVALID_SPECTRUM,
           "CSpectrum::SetFluxAxis: new flux axis has not the "
           "same size than spectral axis");
  }

  ResetContinuum();
  SetType(EType::nType_raw);
  m_rebin->reset();
  m_RawFluxAxis = std::move(fluxaxis);
}

void CSpectrum::SetSpectralAndFluxAxes(CSpectrumSpectralAxis spcaxis,
                                       CSpectrumFluxAxis fluxaxis) {
  if (fluxaxis.GetSamplesCount() != spcaxis.GetSamplesCount() ||
      spcaxis.isEmpty()) {
    THROWG(ErrorCode::INVALID_SPECTRUM,
           "CSpectrum::SetSpectralAndFluxAxes: new flux axis "
           "has not the same size than new spectral axis");
  }

  ResetContinuum();
  SetType(EType::nType_raw);
  m_rebin->reset();

  m_SpectralAxis = std::move(spcaxis);
  SetFluxAxis(std::move(fluxaxis));
}

void CSpectrum::InitSpectrumContinuum(CParameterStore &parameterStore) {
  if (!IsValid())
    THROWG(ErrorCode::INVALID_SPECTRUM,
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
  Log.LogDetail(Formatter() << "Continuum estimation on input spectrum: using "
                            << m_estimationMethod);

  if (m_estimationMethod == "irregularSamplingMedian") {
    CContinuumIrregularSamplingMedian continuum;
    continuum.SetMedianKernelWidth(m_medianWindowSize);
    continuum.SetMeanKernelWidth(m_medianWindowSize);
    continuum.SetMedianEvenReflection(m_medianEvenReflection);
    RemoveContinuum(continuum);
    Log.LogDetail(Formatter() << "Continuum estimation - medianKernelWidth ="
                              << m_medianWindowSize);
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
    THROWG(ErrorCode::INTERNAL_ERROR, "Unknown Estimation method");
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
  if (wlRange.GetBegin() < m_SpectralAxis.GetLambdaRange().GetBegin() ||
      wlRange.GetEnd() > m_SpectralAxis.GetLambdaRange().GetEnd())
    return false;

  const CSpectrumNoiseAxis &error = GetErrorAxis();
  const CSpectrumFluxAxis &flux = GetFluxAxis();

  TInt32Range iRange = m_SpectralAxis.GetIndexesAtWaveLengthRange(wlRange);
  Int32 n = iRange.GetLength() + 1;
  Float64 x[n];
  Float64 y[n];
  Float64 w[n];

  for (Int32 k = 0; k < n; k++) {
    Int32 ik = k + iRange.GetBegin();
    w[k] = 1.0 / (error[ik] * error[ik]);
    x[k] = m_SpectralAxis[ik];
    y[k] = flux[ik];
  }

  double c0, c1, cov00, cov01, cov11, chisq;
  gsl_fit_wlinear(x, 1, w, 1, y, 1, n, &c0, &c1, &cov00, &cov01, &cov11,
                  &chisq);

  a = c1;
  b = c0;
  return true;
}

const std::string &CSpectrum::GetName() const { return m_Name; }

const std::string &CSpectrum::getObsID() const { return m_obsId; }

void CSpectrum::setObsID(const std::string &obsID) { m_obsId = obsID; }

void CSpectrum::SetName(std::string name) { m_Name = std::move(name); }

const CSpectrum::EType CSpectrum::GetType() const { return m_spcType; }

void CSpectrum::SetType(const CSpectrum::EType type) const {
  if (m_spcType != type) {
    m_rebin->reset();
    m_spcType = type;
  }
}

void CSpectrum::ValidateFlux(Float64 LambdaMin, Float64 LambdaMax) const {
  if (IsFluxEmpty())
    THROWG(ErrorCode::INVALID_SPECTRUM_FLUX, "Invalid spectrum: empty flux");

  bool allzero = true;
  Int32 nInvalid = 0;

  Int32 iMin = m_SpectralAxis.GetIndexAtWaveLength(LambdaMin);
  Int32 iMax = m_SpectralAxis.GetIndexAtWaveLength(LambdaMax);
  Log.LogDetail(Formatter()
                << "CSpectrum::ValidateFlux - checking on the configured "
                   "lambdaRange = ("
                << LambdaMin << ", " << LambdaMax << ")");
  Log.LogDetail(Formatter()
                << "CSpectrum::ValidateFlux - checking on the true observed "
                   "spectral axis lambdarange = ("
                << m_SpectralAxis[iMin] << ", " << m_SpectralAxis[iMax] << ")");

  const CSpectrumFluxAxis &flux = GetFluxAxis();
  std::map<string, int> invalidElements = {};

  // check flux
  TBoolList validSamples = flux.checkFlux();
  for (Int32 i = iMin; i < iMax; i++) {
    // collect invalid values
    if (!validSamples[i]) {
      ++nInvalid;
      ++invalidElements[to_string(flux[i])];
    }

    // all zeroes check
    if (flux[i] != 0.0)
      allzero = false;
  }

  if (nInvalid > 0) {
    string errorMessage = "Failed to validate spectrum flux. " +
                          to_string(nInvalid) + " invalid values:";
    for (auto &itr : invalidElements) {
      errorMessage += itr.first + " : " + to_string(itr.second);
    }
    THROWG(ErrorCode::INVALID_SPECTRUM_FLUX, errorMessage);
  }

  if (allzero) {
    THROWG(ErrorCode::INVALID_SPECTRUM_FLUX,
           "Failed to validate spectrum flux: all values are zeroes.");
  }
}

void CSpectrum::ValidateNoise(Float64 LambdaMin, Float64 LambdaMax) const {
  Int32 nInvalid = 0;

  if (IsNoiseEmpty())
    THROWG(ErrorCode::INVALID_NOISE, "Invalid spectrum: empty noise.");

  const TFloat64List &error = GetFluxAxis().GetError().GetSamplesVector();
  Int32 iMin = m_SpectralAxis.GetIndexAtWaveLength(LambdaMin);
  Int32 iMax = m_SpectralAxis.GetIndexAtWaveLength(LambdaMax);
  Log.LogDetail(Formatter()
                << "CSpectrum::ValidateNoise - checking on the configured "
                   "lambdaRange = ("
                << LambdaMin << ", " << LambdaMax << ")");
  Log.LogDetail(Formatter()
                << "CSpectrum::ValidateNoise - checking on the true observed "
                   "spectral axis lambdarange = ("
                << m_SpectralAxis[iMin] << ", " << m_SpectralAxis[iMax] << ")");

  std::map<string, int> invalidElements = {};

  // check noise
  TBoolList validSamples = GetFluxAxis().GetError().checkNoise();
  for (Int32 i = iMin; i < iMax; i++) {
    if (!validSamples[i]) {
      ++nInvalid;
      invalidElements[to_string(error[i])]++;
    }
  }
  if (nInvalid > 0) {
    string errorMessage = "Failed to validate spectrum noise. " +
                          to_string(nInvalid) + " invalid values:";
    for (auto &itr : invalidElements) {
      errorMessage += itr.first + " : " + to_string(itr.second);
    }
    Log.LogDetail(Formatter() << "    CSpectrum::ValidateNoise - Found "
                              << nInvalid << " invalid noise samples");
    THROWG(ErrorCode::INVALID_NOISE, errorMessage);
  }
}

bool CSpectrum::correctSpectrum(Float64 LambdaMin, Float64 LambdaMax,
                                Float64 coeffCorr) {
  if (!IsValid())
    THROWG(ErrorCode::INVALID_SPECTRUM,
           "Invalid spectrum with empty axes, non-matching size "
           "or unsorted spectral axis");

  Int32 iMin = m_SpectralAxis.GetIndexAtWaveLength(LambdaMin);
  Int32 iMax = m_SpectralAxis.GetIndexAtWaveLength(LambdaMax);

  bool corrected =
      GetFluxAxis_().correctFluxAndNoiseAxis(iMin, iMax, coeffCorr);

  if (corrected)
    ResetContinuum();

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
      m_estimationMethod == "irregularSamplingMedian") {
    ResetContinuum();
  }
  m_medianWindowSize = winsize;
}

void CSpectrum::SetMedianEvenReflection(bool medianEvenReflection) {
  if (m_medianEvenReflection != medianEvenReflection &&
      m_estimationMethod == "irregularSamplingMedian") {
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
    THROWG(ErrorCode::INTERNAL_ERROR,
           "Cannot set continuum flux: non-matching size with spectral axis");
  }

  m_ContinuumFluxAxis = ContinuumFluxAxis;
}

void CSpectrum::setRebinInterpMethod(const std::string &opt_interp) const {
  m_rebin = std::move(*m_rebin).convert(opt_interp);
}

// Test methode Rebin

///
/// * This rebin method targets processing speed:
/// - it uses already allocated rebinedFluxAxis, rebinedSpectralAxis and
/// rebinedMask
/// - opt_interp = 'lin' : linear interpolation is performed by default
/// - opt_interp = 'preComputedFineGrid' : nearest grid point interpolation
/// is performed using m_pfgFlux which is the precomputed fine grid
/// - opt_interp = 'spline' : GSL/spline interpolation is performed (TODO -
/// not tested)
/// - opt_interp = 'ngp' : nearest grid point is performed (TODO - not
/// tested)
/**
 * targetSpectralAxis should be expressed in same frame as source
 * SpetralAxis
 */
void CSpectrum::Rebin(const TFloat64Range &range,
                      const CSpectrumSpectralAxis &targetSpectralAxis,
                      CSpectrum &rebinedSpectrum, CMask &rebinedMask,
                      const std::string &opt_error_interp) const {
  if (!IsValid())
    THROWG(ErrorCode::INVALID_SPECTRUM,
           "Invalid spectrum with empty axes, non-matching size "
           "or unsorted spectral axis");

  m_rebin->compute(range, targetSpectralAxis, rebinedSpectrum, rebinedMask,
                   opt_error_interp);
}

void CSpectrum::ApplyAmplitude(Float64 amplitude) {
  m_RawFluxAxis *= amplitude;
  if (alreadyRemoved) {
    m_ContinuumFluxAxis *= amplitude;
    m_WithoutContinuumFluxAxis *= amplitude;
  }
  m_rebin->reset();
}

void CSpectrum::ValidateSpectrum(TFloat64Range lambdaRange,
                                 bool enableInputSpcCorrect) {
  if (!IsValid())
    THROWG(ErrorCode::INVALID_SPECTRUM,
           "Invalid spectrum with empty axes or non-matching "
           "size or unsorted spectral axis");

  TFloat64Range clampedlambdaRange;
  m_SpectralAxis.ClampLambdaRange(lambdaRange, clampedlambdaRange);
  Log.LogInfo(Formatter() << "Validate spectrum: (clamped lambda range: "
                          << clampedlambdaRange.GetBegin() << "-"
                          << clampedlambdaRange.GetEnd() << ":"
                          << GetResolution() << ")");

  Float64 lmin = clampedlambdaRange.GetBegin();
  Float64 lmax = clampedlambdaRange.GetEnd();

  // Check if the Spectrum is valid on the clamped lambdarange
  if (enableInputSpcCorrect)
    if (correctSpectrum(lmin, lmax))
      Log.LogInfo(Formatter()
                  << "Successfully corrected noise on wavelength range (%."
                  << lmin << " ; %." << lmax << ")");

  ValidateFlux(lmin, lmax);
  Log.LogDetail(Formatter()
                << "Successfully validated spectrum flux, on wavelength range "
                   "(%."
                << lmin << " ; %." << lmax << ")");

  // Check if the noise is valid in the clamped lambdaRange
  ValidateNoise(lmin, lmax);
  Log.LogDetail(Formatter()
                << "Successfully validated noise on wavelength range (%."
                << lmin << " ; %." << lmax << ")");

  if (!m_LSF)
    return;
  // check if spectrum LSF spectralAxis covers clamped lambdaRange
  if (!m_LSF->checkAvailability(lmin) || !m_LSF->checkAvailability(lmax)) {
    THROWG(ErrorCode::INVALID_LSF,
           Formatter() << "Failed to validate lsf on wavelength range [" << lmin
                       << ";" << lmax << "]");
  }
}

std::pair<Float64, Float64> CSpectrum::integrateFluxes_usingTrapez(
    CSpectrumSpectralAxis const &spectralAxis,
    CSpectrumFluxAxis const &fluxAxis, TInt32RangeList const &indexRangeList) {

  Float64 sumFlux = 0.0;
  Float64 sumErr = 0.0;

  if (spectralAxis.GetSamplesCount() < 2)
    THROWG(ErrorCode::INTERNAL_ERROR, "Not enough samples to integrate flux");

  if (spectralAxis.GetSamplesCount() != fluxAxis.GetSamplesCount())
    THROWG(ErrorCode::INTERNAL_ERROR,
           "spectral axis and flux axis have different samples number");

  const auto &Error = fluxAxis.GetError();
  for (auto &r : indexRangeList) {
    for (Int32 t = r.GetBegin(), e = r.GetEnd(); t < e; t++) {
      Float64 trapweight = (spectralAxis[t + 1] - spectralAxis[t]) * 0.5;
      sumFlux += trapweight * (fluxAxis[t + 1] + fluxAxis[t]);

      Float64 ea = Error[t] * Error[t];
      Float64 eb = Error[t + 1] * Error[t + 1];
      sumErr += trapweight * trapweight * (eb + ea);
    }
  }
  return std::make_pair(sumFlux, sumErr);
}
