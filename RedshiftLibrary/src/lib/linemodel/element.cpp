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

using namespace NSEpic;

/**
 * \brief Constructs the object setting members according to arguments and
 *defaults.
 **/
CLineModelElement::CLineModelElement(
    std::vector<CLine> rs, const std::string &widthType,
    const std::shared_ptr<TFittedData> &fittedData, TInt32List catalogIndexes)
    : m_Lines(rs), m_fittingGroupInfo("-1"), m_fittedData(fittedData),
      m_OutsideLambdaRangeOverlapThreshold(
          0.33), // 33% overlap minimum in order to keep the line
      m_OutsideLambdaRange(
          true), // example: 0.33 means 66% of the line is allowed to be outside
                 // the spectrum with the line still considered inside the
                 // lambda range

      m_SignFactors(rs.size()),

      m_absLinesLimit(
          1.0) //-1: disable the ABS lines amplitude cut, any other value is
               // used as a limit for the abs line coeff (typically: 1.0)

{
  if (widthType == "instrumentdriven") {
    m_LineWidthType = INSTRUMENTDRIVEN;
  } else if (widthType == "combined") {
    m_LineWidthType = COMBINED;
  } else if (widthType == "velocitydriven") {
    m_LineWidthType = VELOCITYDRIVEN;
  } else {
    THROWG(INTERNAL_ERROR, Formatter() << "Unknown LineWidthType" << widthType);
  }

  Int32 nLines = m_Lines.size();
  for (Int32 i = 0; i < nLines; i++) {
    if (m_Lines[i].GetType() == CLine::nType_Emission)
      m_SignFactors[i] = 1.0;
    else
      m_SignFactors[i] = -1.0;
  }

  m_LineCatalogIndexes = catalogIndexes;

  for (Int32 k2 = 0; k2 < nLines; k2++) {
    if (getLineProfile(k2).isAsym() || getLineProfile(k2).isSymIgm())
      m_asymLineIndices.push_back(k2);
  }

  SetFittedAmplitude(NAN, NAN);
}

Int32 CLineModelElement::findElementIndex(Int32 LineCatalogIndex) const {
  Int32 idx = undefIdx;
  for (Int32 iElts = 0; iElts < m_LineCatalogIndexes.size(); iElts++) {
    if (m_LineCatalogIndexes[iElts] == LineCatalogIndex) {
      idx = iElts;
      break;
    }
  }

  return idx;
}

Int32 CLineModelElement::getLineIndexInCatalog(
    Int32 idxLine, const CLineCatalog::TLineVector &catalog) const {
  Int32 lineIndex = undefIdx;
  lineIndex = m_LineCatalogIndexes[idxLine];
  if (lineIndex < 0 || lineIndex >= catalog.size())
    THROWG(INTERNAL_ERROR, "Line idx out-of-bound");

  return lineIndex;
}

Int32 CLineModelElement::GetSize() const { return m_LineCatalogIndexes.size(); }

bool CLineModelElement::IsOutsideLambdaRange() const {
  return m_OutsideLambdaRange;
}

// redirecting to the lsf method for computing instrument responce
/**
 * Get instrumental response (including source response) from LSF
 * combine quadratically with the instrinsic width of the Line itself. The line
 * width in this case represents to the velocity
 * */
Float64 CLineModelElement::GetLineWidth(Float64 redshiftedlambda, Float64 z,
                                        bool isEmission) const {
  const Float64 c = m_speedOfLightInVacuum;
  Float64 v = isEmission ? getVelocityEmission() : getVelocityAbsorption();
  const Float64 pfsSimuCompensationFactor = 1.0;

  if (!m_LSF)
    THROWG(INTERNAL_ERROR, "LSF object is not initailized.");
  Float64 instrumentSigma = m_LSF->GetWidth(redshiftedlambda);

  Float64 velocitySigma = pfsSimuCompensationFactor * v / c *
                          redshiftedlambda; //, useless /(1+z)*(1+z);
  switch (m_LineWidthType) {
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
    THROWG(INTERNAL_ERROR, Formatter()
                               << "Invalid LSF type " << m_LineWidthType);
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
  switch (m_LineWidthType) {
  case INSTRUMENTDRIVEN:
    return 0.0;
  case COMBINED:
    return v_to_sigma * v_to_sigma * v / sigma * profile_derivSigma;
  case VELOCITYDRIVEN:
    return v_to_sigma * profile_derivSigma;
  default:
    THROWG(INTERNAL_ERROR,
           Formatter() << "Invalid LineWidthType : " << m_LineWidthType);
  }
  return 0.0;
}

void CLineModelElement::SetLSF(const std::shared_ptr<const CLSF> &lsf) {
  m_LSF = lsf;
}

void CLineModelElement::SetVelocityEmission(Float64 vel) {
  m_fittedData->m_VelocityEmission = vel;
}

void CLineModelElement::SetVelocityAbsorption(Float64 vel) {
  m_fittedData->m_VelocityAbsorption = vel;
}

Float64 CLineModelElement::getVelocityEmission() const {
  return m_fittedData->m_VelocityEmission;
}

Float64 CLineModelElement::getVelocityAbsorption() const {
  return m_fittedData->m_VelocityAbsorption;
}

Float64 CLineModelElement::getVelocity() const {
  if (!m_Lines.size())
    return NAN;

  bool isEmissionLine = m_Lines[0].GetIsEmission();
  return isEmissionLine ? getVelocityEmission() : getVelocityAbsorption();
}

void CLineModelElement::setVelocity(Float64 vel) {
  if (!m_Lines.size())
    THROWG(INTERNAL_ERROR, "Empty line model element, could not set velocity");

  if (m_Lines[0].GetIsEmission()) {
    m_fittedData->m_VelocityEmission = vel;
  } else {
    m_fittedData->m_VelocityAbsorption = vel;
  }
}

const CLineProfile &CLineModelElement::getLineProfile(Int32 lineIdx) const {
  if (lineIdx > m_Lines.size() - 1)
    THROWG(INTERNAL_ERROR, "out-of-bound index");
  return m_Lines[lineIdx].GetProfile();
}

// wrapper function
void CLineModelElement::SetAsymfitParams(const TAsymParams &params, Int32 idx) {
  if (!m_asymLineIndices.size())
    return;
  if (idx >= 0)
    m_Lines[idx].SetAsymParams(params);
  else
    for (auto i : m_asymLineIndices)
      m_Lines[i].SetAsymParams(params);

  return;
}

// wrapper function
void CLineModelElement::SetSymIgmParams(const TSymIgmParams &params,
                                        Int32 idx) {
  if (!m_asymLineIndices.size())
    return;
  if (idx >= 0)
    m_Lines[idx].SetSymIgmParams(params);
  else
    for (auto i : m_asymLineIndices)
      m_Lines[i].SetSymIgmParams(params);

  return;
}

// wrapper function
void CLineModelElement::resetAsymfitParams() {
  for (auto i : m_asymLineIndices)
    m_Lines[i].resetAsymFitParams();
}

// wrapper function
TAsymParams CLineModelElement::GetAsymfitParams(Int32 idx) const {
  if (!m_asymLineIndices.size())
    return TAsymParams(); // case where no asymprofile in linecatalog
  return m_Lines[m_asymLineIndices[idx]].GetAsymParams();
}

// wrapper function
TSymIgmParams CLineModelElement::GetSymIgmParams(Int32 idx) const {
  if (!m_asymLineIndices.size())
    return TSymIgmParams(); // case where no asymprofile in linecatalog
  return m_Lines[m_asymLineIndices[idx]].GetSymIgmParams();
}

Float64 CLineModelElement::GetSumCross() const { return m_sumCross; }

void CLineModelElement::SetSumCross(Float64 val) { m_sumCross = val; }

Float64 CLineModelElement::GetDtmFree() const { return m_dtmFree; }

void CLineModelElement::SetDtmFree(Float64 val) { m_dtmFree = val; }

Float64 CLineModelElement::GetSumGauss() const { return m_sumGauss; }

void CLineModelElement::SetSumGauss(Float64 val) { m_sumGauss = val; }

const std::string &CLineModelElement::GetFittingGroupInfo() const {
  return m_fittingGroupInfo;
}

void CLineModelElement::SetFittingGroupInfo(const std::string &val) {
  m_fittingGroupInfo = val;
}

void CLineModelElement::SetOffset(Float64 val) {
  for (auto &line : m_Lines)
    line.SetOffset(val);
}

void CLineModelElement::SetLineProfile(Int32 lineIdx,
                                       CLineProfile_ptr &&profile) {
  m_Lines[lineIdx].SetProfile(std::move(profile));
}

Float64 CLineModelElement::GetFitAmplitude() const {
  return m_fittedData->m_fitAmplitude;
}

/**
 * \brief If the argument is greater than or equal to the size of m_Lines,
 *returns the string "-1". Otherwise returns a call to the m_Lines GetName.
 **/
const std::string &CLineModelElement::GetLineName(Int32 subeIdx) const {
  if (subeIdx >= m_Lines.size())
    THROWG(INTERNAL_ERROR, "invalid index");

  return m_Lines[subeIdx].GetName();
}

/**
 * \brief Returns the content of the m_SignFactors with index equal to the
 *argument.
 **/
Float64 CLineModelElement::GetSignFactor(Int32 subeIdx) const {
  return m_SignFactors[subeIdx];
}

bool CLineModelElement::SetAbsLinesLimit(Float64 limit) {
  m_absLinesLimit = limit;
  return true;
}

/**
 * @brief GetContinuumAtCenterProfile
 * @param subeIdx
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
    Int32 subeIdx, const CSpectrumSpectralAxis &spectralAxis, Float64 redshift,
    const CSpectrumFluxAxis &continuumfluxAxis,
    const TPolynomCoeffs &line_polynomCoeffs) const {
  Float64 mu = GetObservedPosition(subeIdx, redshift);

  Int32 IdxCenterProfile = spectralAxis.GetIndexAtWaveLength(mu);
  if (IdxCenterProfile < 0 ||
      IdxCenterProfile > continuumfluxAxis.GetSamplesCount() - 1) {
    return NAN;
  }

  Float64 cont = continuumfluxAxis[IdxCenterProfile] + line_polynomCoeffs.x0 +
                 line_polynomCoeffs.x1 * spectralAxis[IdxCenterProfile] +
                 line_polynomCoeffs.x2 * spectralAxis[IdxCenterProfile] *
                     spectralAxis[IdxCenterProfile];
  return cont;
}

/**
 * \brief Returns the theoretical support range for the line
 **/
void CLineModelElement::EstimateTheoreticalSupport(
    Int32 subeIdx, const CSpectrumSpectralAxis &spectralAxis, Float64 redshift,
    const TFloat64Range &lambdaRange) {
  Float64 mu = GetObservedPosition(subeIdx, redshift);
  if (!m_LSF->checkAvailability(mu)) {
    m_OutsideLambdaRangeList[subeIdx] = true;
    return;
  }
  Float64 sigma = GetLineWidth(mu, redshift, m_Lines[subeIdx].GetIsEmission());
  Float64 winsize = getLineProfile(subeIdx).GetNSigmaSupport() * sigma;
  TInt32Range supportRange =
      EstimateIndexRange(spectralAxis, mu, lambdaRange, winsize);

  m_StartTheoretical[subeIdx] = supportRange.GetBegin();
  m_EndTheoretical[subeIdx] = supportRange.GetEnd();
  m_StartNoOverlap[subeIdx] = supportRange.GetBegin();
  m_EndNoOverlap[subeIdx] = supportRange.GetEnd();

  if (supportRange.GetBegin() >
      supportRange.GetEnd()) // in this case the line is completely outside the
                             // lambdarange
  {
    m_OutsideLambdaRangeList[subeIdx] = true;
  } else { // in this case the line is completely inside the lambdarange or with
           // partial overlap

    Float64 minLineOverlap = m_OutsideLambdaRangeOverlapThreshold * winsize;
    Float64 startLbda = spectralAxis[m_StartNoOverlap[subeIdx]];
    Float64 endLbda = spectralAxis[m_EndNoOverlap[subeIdx]];

    if (startLbda >= (lambdaRange.GetEnd() - minLineOverlap) ||
        endLbda <= (lambdaRange.GetBegin() + minLineOverlap)) {
      m_OutsideLambdaRangeList[subeIdx] = true;
    } else {
      m_OutsideLambdaRangeList[subeIdx] = false;
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

/**
 * \brief Limits each m_Lines element within the argument lambdaRange, and sets
 *the m_FittedAmplitudes to -1. Sets the global outside lambda range. Inits the
 *fitted amplitude values.
 **/
void CLineModelElement::prepareSupport(
    const CSpectrumSpectralAxis &spectralAxis, Float64 redshift,
    const TFloat64Range &lambdaRange) {
  Int32 nLines = m_Lines.size();
  m_OutsideLambdaRange = true;
  m_LineIsActiveOnSupport.assign(nLines, TInt32List(nLines, 0));
  m_StartNoOverlap.assign(nLines, -1);
  m_EndNoOverlap.assign(nLines, -1);
  m_StartTheoretical.assign(nLines, -1);
  m_EndTheoretical.assign(nLines, -1);
  m_OutsideLambdaRangeList.assign(nLines, true);
  for (Int32 i = 0; i < nLines; i++) {
    EstimateTheoreticalSupport(i, spectralAxis, redshift, lambdaRange);

    // set the global outside lambda range
    m_OutsideLambdaRange = m_OutsideLambdaRange && m_OutsideLambdaRangeList[i];

    // set the lines active on their own support
    m_LineIsActiveOnSupport[i][i] = 1;
  }

  bool supportNoOverlap_has_duplicates = true;
  Int32 x1 = 0;
  Int32 y1 = 0;
  Int32 x2 = 0;
  Int32 y2 = 0;
  Int32 icmpt = 0;
  Int32 ncmpt = 20;
  while (supportNoOverlap_has_duplicates && icmpt < ncmpt) {
    icmpt++;
    for (Int32 i = 0; i < nLines; i++) {
      if (m_OutsideLambdaRangeList[i]) {
        continue;
      }
      if (m_StartNoOverlap[i] > m_EndNoOverlap[i]) {
        continue;
      }
      for (Int32 j = 0; j < nLines; j++) {
        if (m_OutsideLambdaRangeList[j]) {
          continue;
        }
        if (m_StartNoOverlap[j] > m_EndNoOverlap[j]) {
          continue;
        }
        if (i == j) {
          continue;
        }
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
          m_LineIsActiveOnSupport[i][j] = 1;
          m_LineIsActiveOnSupport[j][i] = 1;
          // append all the previously overlapping lines as active on the
          // support
          for (Int32 i2 = 0; i2 < nLines; i2++) {
            if (m_OutsideLambdaRangeList[i2]) {
              continue;
            }
            if (m_LineIsActiveOnSupport[i][i2] == 1) {
              m_LineIsActiveOnSupport[i2][j] = 1;
              m_LineIsActiveOnSupport[j][i2] = 1;
            }
          }
          // append all the previously overlapping lines as active on the
          // support
          for (Int32 j2 = 0; j2 < nLines; j2++) {
            if (m_OutsideLambdaRangeList[j2]) {
              continue;
            }
            if (m_LineIsActiveOnSupport[j][j2] == 1) {
              m_LineIsActiveOnSupport[j2][i] = 1;
              m_LineIsActiveOnSupport[i][j2] = 1;
            }
          }
        }
      }
    }

    supportNoOverlap_has_duplicates = false;

    // check that there are no overlapping sub-supports in the list
    for (Int32 i = 0; i < nLines; i++) {
      if (supportNoOverlap_has_duplicates) {
        break;
      }
      if (m_OutsideLambdaRangeList[i]) {
        continue;
      }
      for (Int32 j = 0; j < nLines; j++) {
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

  // init the fitted amplitude values and related variables
  m_fittedData->m_FittedAmplitudes.assign(nLines, NAN);
  m_fittedData->m_FittedAmplitudeErrorSigmas.assign(nLines, NAN);
  m_fittedData->m_fitAmplitude = NAN;
  m_sumGauss = NAN;
  m_sumCross = NAN;
  m_dtmFree = NAN;
}

/**
 * \brief Creates an empty list of ranges as the return value. If not
 *m_OutsideLambdaRange, for each m_Lines element which is also not outside
 *lambda range, add its support to the return value.
 **/
TInt32RangeList CLineModelElement::getSupport() const {
  TInt32RangeList support;
  Int32 nLines = m_Lines.size();
  if (m_OutsideLambdaRange == false) {
    for (Int32 i = 0; i < nLines; i++) {
      if (m_OutsideLambdaRangeList[i]) {
        continue;
      }
      support.push_back(TInt32Range(m_StartNoOverlap[i], m_EndNoOverlap[i]));
    }
  }
  return support;
}

TInt32RangeList CLineModelElement::getTheoreticalSupport() const {
  TInt32RangeList support;
  Int32 nLines = m_Lines.size();

  if (m_OutsideLambdaRange == false) {
    for (Int32 i = 0; i < nLines; i++) {
      if (m_OutsideLambdaRangeList[i]) {
        continue;
      }
      support.push_back(
          TInt32Range(m_StartTheoretical[i], m_EndTheoretical[i]));
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
TInt32Range CLineModelElement::getSupportSubElt(Int32 subeIdx) const {
  /*
  TInt32Range support;
  if(m_OutsideLambdaRange==false){
      if( !m_OutsideLambdaRangeList[subeIdx] ){
          support = TInt32Range(m_StartNoOverlap[subeIdx],
  m_EndNoOverlap[subeIdx]);
      }
  }
  */
  TInt32Range support =
      TInt32Range(m_StartNoOverlap[subeIdx], m_EndNoOverlap[subeIdx]);
  return support;
}

/**
 * \brief Returns the theoretical support of the line (sub-element).
 **/
TInt32Range
CLineModelElement::getTheoreticalSupportSubElt(Int32 subeIdx) const {
  /*
  TInt32Range support;
  if(m_OutsideLambdaRange==false){
      if( !m_OutsideLambdaRangeList[subeIdx] ){
          support = TInt32Range(m_StartTheoretical[subeIdx],
  m_EndTheoretical[subeIdx]);
      }
  }
  */
  TInt32Range support =
      TInt32Range(m_StartTheoretical[subeIdx], m_EndTheoretical[subeIdx]);
  return support;
}

/**
 * \brief Calls GetLineWidth using the arguments and a calculated argument mu.
 **/
void CLineModelElement::getObservedPositionAndLineWidth(
    Int32 subeIdx, Float64 redshift, Float64 &mu, Float64 &sigma,
    bool doAsymfitdelta) const {
  mu = GetObservedPosition(subeIdx, redshift, doAsymfitdelta);
  if (!m_LSF->checkAvailability(mu)) {
    THROWG(INTERNAL_ERROR, "Line position does not belong to LSF range");
  } else
    sigma = GetLineWidth(mu, redshift, m_Lines[subeIdx].GetIsEmission());
  return;
}

/**
 * \brief Get the observed position of the sub-element subeIdx for a given
 *redshift
 **/
Float64 CLineModelElement::GetObservedPosition(Int32 subeIdx, Float64 redshift,
                                               bool doAsymfitdelta) const {
  Float64 dzOffset = m_Lines[subeIdx].GetOffset() / m_speedOfLightInVacuum;

  Float64 mu = m_Lines[subeIdx].GetPosition() * (1 + redshift) * (1 + dzOffset);

  // deals with delta of asym profile
  if (doAsymfitdelta) {
    mu -= m_Lines[subeIdx].GetProfile().GetDelta();
  }
  return mu;
}

/**
 * \brief Returns the line profile of the sub-element subIdx at wavelength x,
 *for a given redshift.
 **/
Float64 CLineModelElement::GetLineProfileAtRedshift(Int32 subeIdx,
                                                    Float64 redshift,
                                                    Float64 x) const {
  Float64 mu = NAN;
  Float64 sigma = NAN;
  getObservedPositionAndLineWidth(subeIdx, redshift, mu, sigma,
                                  false); // do not apply Lya asym offset
  // check if lineprofile is SYMIGM
  const CLineProfile &profile = getLineProfile(subeIdx);

  return profile.GetLineProfileVal(x, mu, sigma);
}

/**
 * \brief Returns the fitted amplitude for the argument index.
 **/
Float64 CLineModelElement::GetFittedAmplitude(Int32 subeIdx) const {
  return m_fittedData->m_FittedAmplitudes[subeIdx];
}

/**
 * \brief Returns the fitted amplitude error for the argument index.
 **/
Float64 CLineModelElement::GetFittedAmplitudeErrorSigma(Int32 subeIdx) const {
  return m_fittedData->m_FittedAmplitudeErrorSigmas[subeIdx];
}

/**
 * \brief Returns NAN if m_OutsideLambdaRange, and the fitted amplitude /
 *nominal amplitude of the first element otherwise.
 **/
Float64 CLineModelElement::GetElementAmplitude() const {
  if (m_OutsideLambdaRange) {
    return NAN;
  }
  for (Int32 k = 0; k < m_Lines.size(); k++) {
    if (!m_OutsideLambdaRangeList[k] &&
        m_fittedData->m_NominalAmplitudes[k] != 0.0) {
      return m_fittedData->m_FittedAmplitudes[k] /
             m_fittedData->m_NominalAmplitudes[k];
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
  for (Int32 k = 0; k < m_Lines.size(); k++) {
    if (!m_OutsideLambdaRangeList[k] &&
        m_fittedData->m_NominalAmplitudes[k] != 0.0) {
      return m_fittedData->m_FittedAmplitudeErrorSigmas[k] /
             m_fittedData->m_NominalAmplitudes[k];
    }
  }
  return NAN;
}

/**
 * \brief Returns the nominal amplitude with index subeIdx.
 **/
Float64 CLineModelElement::GetNominalAmplitude(Int32 subeIdx) const {
  return m_fittedData->m_NominalAmplitudes[subeIdx];
}

/**
 * \brief Set the nominal amplitude with index subeIdx.
 **/
bool CLineModelElement::SetNominalAmplitude(Int32 subeIdx, Float64 nominalamp) {
  m_fittedData->m_NominalAmplitudes[subeIdx] = nominalamp;
  return true;
}

// WARNING: setfittedamplitude should not be applied to sub element individually
// ?
void CLineModelElement::SetFittedAmplitude(Int32 subeIdx, Float64 A,
                                           Float64 SNR) {
  if (subeIdx == undefIdx || subeIdx >= m_fittedData->m_FittedAmplitudes.size())
    THROWG(INTERNAL_ERROR, "out-of-bound index");

  if (m_OutsideLambdaRangeList[subeIdx] || std::isnan(A)) {
    m_fittedData->m_FittedAmplitudes[subeIdx] = NAN;
    m_fittedData->m_FittedAmplitudeErrorSigmas[subeIdx] = NAN;
    return;
  }
  m_fittedData->m_FittedAmplitudes[subeIdx] =
      A * m_fittedData->m_NominalAmplitudes[subeIdx];
  // limit the absorption to 0.0-1.0, so that it's never <0
  //*
  if (m_SignFactors[subeIdx] == -1 && m_absLinesLimit > 0.0 &&
      m_fittedData->m_FittedAmplitudes[subeIdx] > m_absLinesLimit) {
    m_fittedData->m_FittedAmplitudes[subeIdx] = m_absLinesLimit;
  }

  m_fittedData->m_FittedAmplitudeErrorSigmas[subeIdx] =
      SNR * m_fittedData->m_NominalAmplitudes[subeIdx]; // todo: check correct
                                                        // formulation for Error
}

/**
 * \brief If outside lambda range, sets fitted amplitudes and errors to -1. If
 *inside, sets each line's fitted amplitude and error to -1 if line outside
 *lambda range, or amplitude to A * nominal amplitude and error to SNR * nominal
 *amplitude.
 **/
void CLineModelElement::SetFittedAmplitude(Float64 A, Float64 SNR) {
  if (std::isnan(A) || m_OutsideLambdaRange) {
    m_fittedData->m_fitAmplitude = NAN;
    m_fittedData->m_FittedAmplitudes.assign(m_Lines.size(), NAN);
    m_fittedData->m_FittedAmplitudeErrorSigmas.assign(m_Lines.size(), NAN);
    return;
  }

  m_fittedData->m_fitAmplitude =
      std::max(0., A); // force amplitudes to zero if negative
  for (Int32 k = 0; k < m_Lines.size(); k++) {
    if (m_OutsideLambdaRangeList[k]) {
      m_fittedData->m_FittedAmplitudes[k] = NAN;
      m_fittedData->m_FittedAmplitudeErrorSigmas[k] = NAN;
      continue;
    }
    m_fittedData->m_FittedAmplitudes[k] =
        A * m_fittedData->m_NominalAmplitudes[k];
    // limit the absorption to 0.0-1.0, so that it's never <0

    if (m_SignFactors[k] == -1 && m_absLinesLimit > 0.0 &&
        m_fittedData->m_FittedAmplitudes[k] > m_absLinesLimit) {
      m_fittedData->m_FittedAmplitudes[k] = m_absLinesLimit;
    }

    m_fittedData->m_FittedAmplitudeErrorSigmas[k] =
        SNR * m_fittedData->m_NominalAmplitudes[k]; // todo: check correct
                                                    // formulation for Error
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
    Int32 lineIdx) const {
  if (m_OutsideLambdaRange)
    return;

  Int32 nLines = m_Lines.size();
  for (Int32 k = 0; k < nLines; k++) { // loop on the interval
    if (m_OutsideLambdaRangeList[k])
      continue;

    if (lineIdx != undefIdx && !(m_LineIsActiveOnSupport[k][lineIdx]))
      continue;

    for (Int32 i = m_StartNoOverlap[k]; i <= m_EndNoOverlap[k]; i++) {
      Float64 lambda = modelspectralAxis[i];
      Float64 Yi = getModelAtLambda(lambda, redshift, continuumfluxAxis[i], k);
      modelfluxAxis[i] += Yi;
      if (std::isnan(modelfluxAxis[i]))
        THROWG(INTERNAL_ERROR,
               Formatter() << "addToSpectrumModel has a NaN flux Line" << k
                           << ": ContinuumFlux " << continuumfluxAxis[i]
                           << ", ModelAtLambda Yi = " << Yi << " for range ["
                           << m_StartNoOverlap[k] << ", " << m_EndNoOverlap[k]
                           << "]");
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

  for (Int32 k = 0; k < m_Lines.size(); k++) { // loop on the interval
    if (m_OutsideLambdaRangeList[k])
      continue;

    if ((emissionLine ^ m_Lines[k].GetIsEmission()))
      continue;

    Float64 A = m_fittedData->m_FittedAmplitudes[k];
    if (std::isnan(A))
      continue;
    // THROWG(INTERNAL_ERROR,"FittedAmplitude cannot
    // be NAN");//to be uncommented

    for (Int32 i = m_StartNoOverlap[k]; i <= m_EndNoOverlap[k]; i++) {

      Float64 x = modelspectralAxis[i];
      Float64 mu = NAN;
      Float64 sigma = NAN;
      getObservedPositionAndLineWidth(k, redshift, mu, sigma, false);

      if (m_SignFactors[k] == -1)
        modelfluxAxis[i] +=
            m_SignFactors[k] * A * continuumfluxAxis[i] *
            GetLineProfileDerivVel(getLineProfile(k), x, mu, sigma,
                                   m_Lines[k].GetIsEmission());
      else
        modelfluxAxis[i] +=
            m_SignFactors[k] * A *
            GetLineProfileDerivVel(getLineProfile(k), x, mu, sigma,
                                   m_Lines[k].GetIsEmission());
    }
  }
  return;
}

/**
 * \brief Returns the sum of the amplitude of each line on redshifted lambda.
 **/
Float64 CLineModelElement::getModelAtLambda(Float64 lambda, Float64 redshift,
                                            Float64 continuumFlux,
                                            Int32 kLineSupport) const {
  if (m_OutsideLambdaRange) {
    return 0.0;
  }
  Float64 Yi = 0.0;

  Float64 x = lambda;
  Int32 nLines = m_Lines.size();

  for (Int32 k2 = 0; k2 < nLines; k2++) // loop on lines
  {
    if (m_OutsideLambdaRangeList[k2]) {
      continue;
    }
    if (kLineSupport >= 0 && m_LineIsActiveOnSupport[k2][kLineSupport] == 0) {
      continue;
    }

    Float64 A = m_fittedData->m_FittedAmplitudes[k2];
    if (std::isnan(A))
      continue;
    // THROWG(INTERNAL_ERROR,"FittedAmplitude cannot
    // be NAN");
    if (A < 0.)
      continue;

    Float64 fluxval =
        m_SignFactors[k2] * A * GetLineProfileAtRedshift(k2, redshift, x);
    Yi += m_SignFactors[k2] == -1 ? continuumFlux * fluxval : fluxval;

    if (std::isnan(Yi))
      THROWG(INTERNAL_ERROR, Formatter()
                                 << "NaN fluxval for Line nb: " << k2
                                 << " and GetLineProfileAtRedshift: "
                                 << GetLineProfileAtRedshift(k2, redshift, x));
  }
  return Yi;
}

Float64 CLineModelElement::GetModelDerivAmplitudeAtLambda(
    Float64 lambda, Float64 redshift, Float64 continuumFlux) const {
  if (m_OutsideLambdaRange) {
    return 0.0;
  }
  Float64 Yi = 0.0;

  Float64 x = lambda;

  for (Int32 k2 = 0; k2 < m_Lines.size(); k2++) // loop on lines
  {
    if (m_OutsideLambdaRangeList[k2]) {
      continue;
    }

    Float64 fluxval = m_SignFactors[k2] *
                      m_fittedData->m_NominalAmplitudes[k2] *
                      GetLineProfileAtRedshift(k2, redshift, x);
    Yi += m_SignFactors[k2] == -1 ? continuumFlux * fluxval : fluxval;
  }
  return Yi;
}

Float64 CLineModelElement::GetModelDerivContinuumAmpAtLambda(
    Float64 lambda, Float64 redshift, Float64 continuumFluxUnscale) const {
  if (m_OutsideLambdaRange) {
    return 0.0;
  }
  Float64 Yi = 0.0;

  Float64 x = lambda;

  for (Int32 k2 = 0; k2 < m_Lines.size(); k2++) // loop on lines
  {
    if (m_OutsideLambdaRangeList[k2]) {
      continue;
    }

    if (m_SignFactors[k2] == 1) {
      continue;
    }

    Float64 A = m_fittedData->m_FittedAmplitudes[k2];
    if (std::isnan(A))
      continue;
    // THROWG(INTERNAL_ERROR,"FittedAmplitude cannot
    // be NAN");

    Yi += m_SignFactors[k2] * continuumFluxUnscale * A *
          GetLineProfileAtRedshift(k2, redshift, x);
  }
  return Yi;
}

/* Given the value of the partial deriv of the flux of this multiline at the
 * given lamda when The continuum is not a variable of z
 */
Float64 CLineModelElement::GetModelDerivZAtLambdaNoContinuum(
    Float64 lambda, Float64 redshift, Float64 continuumFlux) const {
  if (m_OutsideLambdaRange) {
    return 0.0;
  }
  Float64 Yi = 0.0;

  Float64 x = lambda;

  for (Int32 k2 = 0; k2 < m_Lines.size(); k2++) // loop on lines
  {
    if (m_OutsideLambdaRangeList[k2]) {
      continue;
    }
    Float64 A = m_fittedData->m_FittedAmplitudes[k2];
    if (std::isnan(A))
      continue;
    // THROWG(INTERNAL_ERROR,"FittedAmplitude cannot
    // be NAN");

    Float64 dzOffset = m_Lines[k2].GetOffset() / m_speedOfLightInVacuum;
    Float64 lamdba0 = m_Lines[k2].GetPosition() * (1 + dzOffset);
    Float64 mu = NAN;
    Float64 sigma = NAN;
    getObservedPositionAndLineWidth(k2, redshift, mu, sigma,
                                    false); // do not apply Lya asym offset

    const CLineProfile &profile = getLineProfile(k2);
    Float64 lineprofile_derivz =
        profile.GetLineProfileDerivZ(x, lamdba0, redshift, sigma);

    Float64 fluxval = m_SignFactors[k2] * A * lineprofile_derivz;
    Yi += m_SignFactors[k2] == -1 ? continuumFlux * fluxval : fluxval;
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
  if (m_OutsideLambdaRange) {
    return 0.0;
  }
  Float64 Yi = 0.0;

  Float64 x = lambda;

  for (Int32 k2 = 0; k2 < m_Lines.size(); k2++) // loop on lines
  {
    if (m_OutsideLambdaRangeList[k2]) {
      continue;
    }
    Float64 A = m_fittedData->m_FittedAmplitudes[k2];
    if (std::isnan(A))
      continue;
    // THROWG(INTERNAL_ERROR,"FittedAmplitude cannot
    // be NAN");

    Float64 dzOffset = m_Lines[k2].GetOffset() / m_speedOfLightInVacuum;
    Float64 lamdba0 = m_Lines[k2].GetPosition() * (1 + dzOffset);
    Float64 mu = NAN;
    Float64 sigma = NAN;
    getObservedPositionAndLineWidth(k2, redshift, mu, sigma,
                                    false); // do not apply Lya asym offset

    const CLineProfile &profile = getLineProfile(k2);
    Float64 lineprofile_derivz = NAN;
    Float64 lineprofile = NAN;
    lineprofile = profile.GetLineProfileVal(x, mu, sigma);
    lineprofile_derivz =
        profile.GetLineProfileDerivZ(x, lamdba0, redshift, sigma);

    Float64 fluxval = m_SignFactors[k2] * A * lineprofile_derivz;
    if (m_SignFactors[k2] == 1) {
      Yi += fluxval;
    } else {
      Yi += continuumFlux * fluxval +
            m_SignFactors[k2] * A * continuumFluxDerivZ * lineprofile;
    }
  }
  return Yi;
}

/**
 * \brief For lines inside lambda range, sets the flux to the continuum flux.
 **/
void CLineModelElement::initSpectrumModel(
    CSpectrumFluxAxis &modelfluxAxis,
    const CSpectrumFluxAxis &continuumfluxAxis, Int32 lineIdx) const {

  if (m_OutsideLambdaRange)
    return;

  for (Int32 k = 0; k < m_Lines.size(); k++) { // loop on the interval
    if (m_OutsideLambdaRangeList[k])
      continue;

    if (lineIdx != undefIdx && !(m_LineIsActiveOnSupport[k][lineIdx]))
      continue;

    for (Int32 i = m_StartNoOverlap[k]; i <= m_EndNoOverlap[k]; i++)
      modelfluxAxis[i] = continuumfluxAxis[i];
  }
  return;
}

/**
 * \brief Returns the index corresponding to the first line whose GetName method
 *returns LineTagStr.
 **/
Int32 CLineModelElement::findElementIndex(const std::string &LineTagStr) const {
  Int32 idx = undefIdx;
  Int32 lines = m_Lines.size();
  for (Int32 iElts = 0; iElts < lines; iElts++) {
    std::string name = m_Lines[iElts].GetName();
    std::size_t foundstra = name.find(LineTagStr.c_str());

    if (foundstra != std::string::npos) {
      idx = iElts;
      break;
    }
  }
  return idx;
}

/**
 * \brief If the fitted amplitude of line with index subeIdx is above the limit,
 *sets it to either that limit or zero, whichever is greater.
 **/
void CLineModelElement::LimitFittedAmplitude(Int32 subeIdx, Float64 limit) {

  if (m_fittedData->m_FittedAmplitudes[subeIdx] > limit) {
    m_fittedData->m_FittedAmplitudes[subeIdx] = std::max(0.0, limit);

    // now update the amplitude of the other lines
    Float64 amplitudeRef = m_fittedData->m_FittedAmplitudes[subeIdx] /
                           m_fittedData->m_NominalAmplitudes[subeIdx];
    for (Int32 k = 0; k < m_Lines.size(); k++) {
      m_fittedData->m_FittedAmplitudes[k] =
          m_fittedData->m_NominalAmplitudes[k] * amplitudeRef;
    }
  }
  return;
}

/**
 * \brief Returns whether the line with index subeIdx is outside the lambda
 *range.
 **/
bool CLineModelElement::IsOutsideLambdaRange(
    Int32 subeIdx) const // IsOutsideLambdaRange
{
  return m_OutsideLambdaRangeList[subeIdx];
}

void CLineModelElement::debug(std::ostream &os) const {
  os << m_fittingGroupInfo << "\t";
  os << getVelocityAbsorption() << "\t" << getVelocityEmission() << "\n";
  /*
  for(int i=0;i<GetSize();i++)
    {
      os << GetLineName(i)<< "\t"<< GetNominalAmplitude(i)<< "\t" <<
  getLineProfile(i).GetName()<<"\t"<< GetSignFactor(i)<< "\n";
    }
  */
}

void CLineModelElement::dumpElement(std::ostream &os) const {
  //  debug(os); // to dump lines info
  os << "m_OutsideLambdaRange\t" << m_OutsideLambdaRange << "\n";
  os << "m_fittingGroupInfo\t" << m_fittingGroupInfo << "\n";
  os << "m_ElementType\t" << CLine::ETypeString.at(GetElementType()) << "\n";

  os << "m_VelocityEmission\t" << getVelocityEmission() << "\n";
  os << "m_VelocityAbsorption\t" << getVelocityAbsorption() << "\n";
  os << "m_OutsideLambdaRangeOverlapThreshold\t"
     << m_OutsideLambdaRangeOverlapThreshold << "\n";
  os << "m_sumCross\t" << m_sumCross << "\n";
  os << "m_sumGauss\t" << m_sumGauss << "\n";
  os << "m_dtmFree\t" << m_dtmFree << "\n";

  os << "m_absLinesLimit\t" << m_absLinesLimit << "\n";

  os << "m_LineIsActiveOnSupport \n";

  for (Int32 i = 0; i < m_LineIsActiveOnSupport.size(); i++)
    for (Int32 j = 0; j < m_LineIsActiveOnSupport.size(); j++)
      os << i << "\t" << j << "\t" << m_LineIsActiveOnSupport[i][j] << "\n";

  os << "\n";
  os << "Line \t OutsideLR\t Sign \t Amp \t AmpErrSigma \t NomAmp \t StartNO "
        "\t EndNO \t "
        "StartTheo \t EndTheo\n";

  for (Int32 i = 0; i < m_Lines.size(); i++) {
    os << i << "\t " << m_OutsideLambdaRangeList[i] << "\t " << m_SignFactors[i]
       << "\t " << m_fittedData->m_FittedAmplitudes[i] << "\t"
       << m_fittedData->m_FittedAmplitudeErrorSigmas[i] << "\t"
       << m_fittedData->m_NominalAmplitudes[i] << "\t" << m_StartNoOverlap[i]
       << "\t" << m_EndNoOverlap[i] << "\t" << m_StartTheoretical[i] << "\t"
       << m_EndTheoretical[i] << "\n";
  }
  os << "\n";
  os << "m_LineCatalogIndexes \n";
  for (Int32 i = 0; i < m_LineCatalogIndexes.size(); i++)
    os << i << "\t" << m_LineCatalogIndexes[i] << "\n";

  os << "\n";
  os << "m_asymLineIndices \n";
  for (Int32 i = 0; i < m_asymLineIndices.size(); i++)
    os << i << "\t" << m_asymLineIndices[i] << "\n";
}

bool CLineModelElement::isLineActiveOnSupport(Int32 lineIdxA,
                                              Int32 lineIdxB) const {
  return m_LineIsActiveOnSupport[lineIdxA][lineIdxB] == 1;
}
