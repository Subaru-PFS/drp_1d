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
#include "RedshiftLibrary/linemodel/spectrummodel.h"
#include "RedshiftLibrary/continuum/irregularsamplingmedian.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/linemodel/element.h"
#include "RedshiftLibrary/operator/powerlaw.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;
using namespace std;

// make a wrapper for this ?
CSpectrumModel::CSpectrumModel(
    const std::shared_ptr<CLineModelElementList> &elements,
    const std::shared_ptr<const CSpectrum> &spc, const CLineMap &restLineList,
    const std::shared_ptr<CContinuumModelSolution> &continuumModelSolution,
    const std::shared_ptr<COperatorContinuumFitting> &continuumFittingOperator,
    Int32 spcIndex)
    : m_Elements(elements), m_inputSpc(spc), m_SpectrumModel(*(spc)),
      m_RestLineList(restLineList), m_fitContinuum(continuumModelSolution),
      m_continuumFittingOperator(continuumFittingOperator),
      m_spcIndex(spcIndex) {
  const Int32 spectrumSampleCount = m_inputSpc->GetSampleCount();
  m_SpcFluxAxis.SetSize(spectrumSampleCount);
  m_spcFluxAxisNoContinuum.SetSize(spectrumSampleCount);
  m_spcFluxAxisNoContinuum.setError(m_inputSpc->GetFluxAxis().GetError());
  m_ContinuumFluxAxis.SetSize(spectrumSampleCount);
}

/**
 * \brief Returns a pointer to m_SpectrumModel.
 **/
const CSpectrum &CSpectrumModel::GetModelSpectrum() const {
  return m_SpectrumModel;
}

/**
 * \brief Returns a pointer to the (re-)estimated continuum flux.
 **/
const CSpectrumFluxAxis &CSpectrumModel::GetModelContinuum() const {
  return m_ContinuumFluxAxis;
}

void CSpectrumModel::initModelWithContinuum() {
  m_SpectrumModel.SetFluxAxis(m_ContinuumFluxAxis);
  for (Int32 i = 0; i < m_SpectrumModel.GetSampleCount(); i++) {
    m_spcFluxAxisNoContinuum[i] = m_SpcFluxAxis[i] - m_ContinuumFluxAxis[i];
  }
}

/**
 * \brief Init the argument elements from the spectrum model with continuum.
 **/
// TODO [opt] feasible without instantiating a CSpectrumFluxAxis and move it
void CSpectrumModel::reinitModelUnderElements(const TInt32List &filterEltsIdx,
                                              Int32 lineIdx) {
  CSpectrumFluxAxis modelFluxAxis = m_SpectrumModel.GetFluxAxis();
  // init spectrum model with continuum
  for (Int32 iElts : filterEltsIdx)
    (*m_Elements)[iElts]->initSpectrumModel(modelFluxAxis, m_ContinuumFluxAxis,
                                            lineIdx);
  m_SpectrumModel.SetFluxAxis(std::move(modelFluxAxis));
}

void CSpectrumModel::refreshModel(CLine::EType lineTypeFilter) {
  reinitModel();
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  CSpectrumFluxAxis modelFluxAxis = m_SpectrumModel.GetFluxAxis();

  if (m_enableAmplitudeOffsets) {
    // add amplitude offsets
    m_Elements->addToSpectrumAmplitudeOffset(m_SpectrumModel.GetSpectralAxis(),
                                             modelFluxAxis);
  }

  // create spectrum model
  Int32 nElements = m_Elements->size();
  for (Int32 iElts = 0; iElts < nElements; iElts++) {
    auto const lineType =
        (*m_Elements)[iElts]->getElementParam()->GetElementType();
    if (lineTypeFilter == CLine::EType::nType_All ||
        lineTypeFilter == lineType) {
      (*m_Elements)[iElts]->addToSpectrumModel(spectralAxis, modelFluxAxis,
                                               m_ContinuumFluxAxis, m_Redshift);
    }
  }

  m_SpectrumModel.SetFluxAxis(std::move(modelFluxAxis));
}

/**
 * \brief Adds a new model to each m_Elements entry specified on the argument.
 * Works as refreshModel.
 **/
void CSpectrumModel::refreshModelUnderElements(const TInt32List &filterEltsIdx,
                                               Int32 lineIdx) {
  reinitModelUnderElements(filterEltsIdx, lineIdx);
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  CSpectrumFluxAxis modelFluxAxis = m_SpectrumModel.GetFluxAxis();
  // create spectrum model
  for (Int32 iElts : filterEltsIdx)
    (*m_Elements)[iElts]->addToSpectrumModel(
        spectralAxis, modelFluxAxis, m_ContinuumFluxAxis, m_Redshift, lineIdx);

  m_SpectrumModel.SetFluxAxis(std::move(modelFluxAxis));
}

/**
 * \brief this function estimates the continuum after removal(interpolation)
 *of the flux samples under the lines for a given redshift
 * //todo: use lambdaRange in order to speed up the continuum estimation
 *process. (rmq: use lambdarange with some margin in order to not detetriorate
 *cont. estimation close to the borders)
 **/
void CSpectrumModel::EstimateSpectrumContinuum(Float64 opt_enhance_lines) {
  TInt32List validEltsIdx = m_Elements->GetModelValidElementsIndexes();
  // TInt32List xInds = getSupportIndexes( validEltsIdx );
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  const auto &ContinuumFluxAxis = m_ContinuumFluxAxis;

  // create new spectrum, which is corrected under the lines

  /*
  //1. interp from previous values
  for( Int32 t=0;t<spectralAxis.GetSamplesCount();t++)
  {
      Y[t] = m_SpcFluxAxis[t];
  }
  Int32 idx = 0;
  Float64 valf = m_SpcFluxAxis[0];
  Float64 corrRatio=0.8;
  for (int i = 0; i < xInds.size(); i++)
  {
      idx = xInds[i];
      if(idx>0){
          if ( std::find(xInds.begin(), xInds.end(), idx-1) == xInds.end() )
          {
              valf=m_SpcFluxAxis[idx-1];
          }
      }
      Y[idx]= corrRatio*valf + (1.0-corrRatio)*m_SpcFluxAxis[idx];
  }
  //*/

  // 2. subtract lines from model
  // model for subtraction
  CSpectrumFluxAxis spcmodel4linefittingFluxAxis =
      m_SpectrumModel.GetFluxAxis();
  // CSpectrum spcmodel4linefitting = GetModelSpectrum();
  for (Int32 i = 0; i < spcmodel4linefittingFluxAxis.GetSamplesCount(); i++) {
    spcmodel4linefittingFluxAxis[i] -= ContinuumFluxAxis[i];
  }
  // optionnaly enhance the abs model component
  if (opt_enhance_lines > 0.0) {
    bool filterOnlyAbs = false;
    for (Int32 t = 0; t < spectralAxis.GetSamplesCount(); t++)
      spcmodel4linefittingFluxAxis[t] *=
          (filterOnlyAbs && spcmodel4linefittingFluxAxis[t] >= 0.0)
              ? 1.0
              : opt_enhance_lines;
  }

  // subtract the lines component
  CSpectrumFluxAxis fluxAxisNothingUnderLines = m_SpcFluxAxis;
  for (Int32 t = 0; t < spectralAxis.GetSamplesCount(); t++) {
    fluxAxisNothingUnderLines[t] -= as_const(spcmodel4linefittingFluxAxis)[t];
  }

  // evaluate contiuum
  CSpectrum spcCorrectedUnderLines(*m_inputSpc);
  spcCorrectedUnderLines.SetFluxAxis(std::move(fluxAxisNothingUnderLines));
  m_ContinuumFluxAxis = spcCorrectedUnderLines.GetContinuumFluxAxis();
}

/**
 * \brief Returns a pointer to a spectrum containing the observed spectrum with
 *the fitted lines subtracted
 **/
CSpectrum CSpectrumModel::GetObservedSpectrumWithLinesRemoved(
    CLine::EType lineTypeFilter) {
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  const CSpectrumFluxAxis &fluxAxis = m_SpectrumModel.GetFluxAxis();

  if (lineTypeFilter == CLine::EType::nType_Emission)
    refreshModel(CLine::EType::nType_Absorption);
  else if (lineTypeFilter == CLine::EType::nType_Absorption)
    refreshModel(CLine::EType::nType_Emission);
  TFloat64List fluxAndContinuum = fluxAxis.GetSamplesVector();

  refreshModel(lineTypeFilter);

  CSpectrum spcCorrectedUnderLines(*m_inputSpc);

  TAxisSampleList Y(m_inputSpc->GetSampleCount());
  const auto &SpcFluxAxis = m_SpcFluxAxis;
  const auto &ContinuumFluxAxis = m_ContinuumFluxAxis;
  for (Int32 t = 0; t < spectralAxis.GetSamplesCount(); t++) {
    Y[t] = SpcFluxAxis[t] - fluxAxis[t] + ContinuumFluxAxis[t];
  }

  // apply smooth blending of continuun below lines
  Float64 alphaMax = 0.9; // alpha blend = 0: only lineSubtractedFlux,
                          // alpha=1: only continuum

  TInt32List nonZeroValidEltsIdx =
      m_Elements->getValidElementIndices(lineTypeFilter);

  TInt32List supportIdxes = m_Elements->getSupportIndexes(nonZeroValidEltsIdx);

  if (supportIdxes.size() > 0) {
    for (Int32 idx : supportIdxes) {
      Float64 weighting =
          GetWeightingAnyLineCenterProximity(idx, nonZeroValidEltsIdx);
      Float64 alpha = alphaMax * weighting;
      Y[idx] = (1. - alpha) * Y[idx] + alpha * fluxAndContinuum[idx];
    }
  }

  spcCorrectedUnderLines.SetFluxAxis(CSpectrumFluxAxis(std::move(Y)));

  return spcCorrectedUnderLines;
}

Float64 CSpectrumModel::GetWeightingAnyLineCenterProximity(
    Int32 sampleIndex, const TInt32List &EltsIdx) const {
  Float64 maxWeight = 0.0;
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  Float64 currentLbda = spectralAxis[sampleIndex];

  for (const Int32 iElts : EltsIdx) {
    for (const auto &range : (*m_Elements)[iElts]->getTheoreticalSupport()) {
      Float64 weight = 0.;
      if (sampleIndex > range.GetBegin() && sampleIndex < range.GetEnd()) {
        Float64 minRangeLbda = spectralAxis[range.GetBegin()];
        Float64 maxRangeLbda = spectralAxis[range.GetEnd()];

        Float64 fullInterval = maxRangeLbda - minRangeLbda;
        if (fullInterval > 0.0) {
          Float64 leftInterval = currentLbda - minRangeLbda;
          Float64 rightInterval = maxRangeLbda - currentLbda;
          Float64 coeff = std::abs(leftInterval - rightInterval);

          weight = 1. - (Float64)coeff / (Float64)fullInterval;
        }
      }

      if (maxWeight < weight)
        maxWeight = weight;
    }
  }

  return maxWeight;
}

std::pair<TInt32Range, TFloat64List>
CSpectrumModel::GetLineRangeAndProfile(Int32 eIdx, Int32 line_id,
                                       Float64 redshift) const {
  auto const &elt = (*m_Elements)[eIdx];
  auto const &spectralAxis = m_SpectrumModel.GetSpectralAxis();

  auto const &[mu, sigma] =
      elt->getObservedPositionAndLineWidth(redshift, line_id);

  // compute averaged continuum under line, weighted by the line profile
  auto const &profile = elt->getElementParam()->getLineProfile(line_id);
  // auto const &polynomCoeffs = elt->GetPolynomCoeffs();

  Float64 const winsize = profile->GetNSigmaSupport() * sigma;
  TInt32Range const indexRange = elt->EstimateIndexRange(
      spectralAxis, mu, spectralAxis.GetLambdaRange(), winsize);

  TFloat64List weights;
  weights.reserve(indexRange.GetLength() + 1);
  for (Int32 lambda_idx = indexRange.GetBegin();
       lambda_idx <= indexRange.GetEnd(); ++lambda_idx) {
    Float64 const lambda = spectralAxis[lambda_idx];
    weights.push_back(profile->GetLineProfileVal(lambda, mu, sigma));
  }

  return std::make_pair(indexRange, weights);
}

/**
 * @brief GetContinuumWeightedSumInRange
 * @param indexRange
 * Add up polynome contribution below lines
 * @return the continuum flux val at the sub element center wavelength.
 */
std::tuple<Float64, Float64, Float64>
CSpectrumModel::GetContinuumWeightedSumInRange(
    TInt32Range const &indexRange, TFloat64List const &weights,
    const TPolynomCoeffs &polynomCoeffs) const {

  auto const &spectralAxis = m_SpectrumModel.GetSpectralAxis();

  Float64 weighted_sum = 0.;
  Float64 total_weight = 0.;
  Float64 total_weight_square = 0.;
  for (size_t idx = 0; idx <= indexRange.GetLength(); ++idx) {
    Int32 const lambda_idx = indexRange.GetBegin() + idx;
    Float64 const lambda = spectralAxis[lambda_idx];
    Float64 const w = weights[idx];
    Float64 const cont =
        m_enableAmplitudeOffsets
            ? m_ContinuumFluxAxis[lambda_idx] + polynomCoeffs.getValue(lambda)
            : m_ContinuumFluxAxis[lambda_idx];
    weighted_sum += w * cont;
    total_weight += w;
    total_weight_square += w * w;
  }

  return std::make_tuple(weighted_sum, total_weight, total_weight_square);
}

/**
 * @brief GetContinuumSquaredResidualInRange
 * Estimate the squared residual (sum of squared residual) on the continuum
 * in a given wavelength index range and return the squared residual and the
 * number of terms sumed
 * 1. calculate the observed spectrum flux with lines subtracted (fitted line
 * model)
 * 2. estimate the squared summ in the given range
 * @param indexRange
 * @return NAN: not enough samples found for error estimation
 */
std::pair<Float64, Int32> CSpectrumModel::getContinuumSquaredResidualInRange(
    TInt32Range const &indexRange) {

  const CSpectrum noLinesSpectrum = GetObservedSpectrumWithLinesRemoved();
  const CSpectrumFluxAxis &noLinesFluxAxis = noLinesSpectrum.GetFluxAxis();
  const auto &ContinuumFluxAxis = m_ContinuumFluxAxis;

  // estimate sum square error between continuum and nolines-spectrum
  Float64 sum = 0.0;
  Int32 nsum = 0;
  for (Int32 t = indexRange.GetBegin(); t <= indexRange.GetEnd(); t++) {
    Float64 diff = noLinesFluxAxis[t] - ContinuumFluxAxis[t];
    sum += diff * diff;
    nsum++;
  }
  return std::make_pair(sum, nsum);
}

/**
 * \brief Returns the error of the support for subelements under the element
 *with the argument eltId as index. Accumulate "fit", the squared difference
 *between model and spectrum, divided by the square of the m_ErrorNoContinuum
 *value. Accumulate "sumErr" 1 / square of the m_ErrorNoContinuum value. return
 *the square root of fit / sumErr.
 **/
std::pair<Float64, Float64>
CSpectrumModel::getModelSquaredResidualUnderElements(TInt32List const &EltsIdx,
                                                     bool with_continuum,
                                                     bool with_weight) const {
  // before elementlistcutting this variable was
  // CElementList::m_ErrorNoContinuum, a reference initialized twice in
  // CElementList constructor, first init to m_spcFluxAxisNoContinuum.GetError()
  // and after to spectrumFluxAxis.GetError
  const CSpectrumNoiseAxis &errorNoContinuum =
      m_SpectrumModel.GetFluxAxis().GetError();
  const CSpectrumFluxAxis &fluxRef =
      with_continuum ? getSpcFluxAxis() : getSpcFluxAxisNoContinuum();

  if (EltsIdx.empty())
    return std::make_pair(NAN, NAN);

  Float64 fit = 0.0;
  const TAxisSampleList &Ymodel =
      m_SpectrumModel.GetFluxAxis().GetSamplesVector();
  const TAxisSampleList &Yspc = fluxRef.GetSamplesVector();
  Float64 diff = 0.0;
  Float64 sumErr = 0.0;

  TInt32List xInds = m_Elements->getSupportIndexes(EltsIdx);
  for (Int32 const j : xInds) {
    diff = (Yspc[j] - Ymodel[j]);
    Float64 const w =
        with_weight ? 1.0 / (errorNoContinuum[j] * errorNoContinuum[j]) : 1.0;
    fit += (diff * diff) * w;
    sumErr += w;
  }
  return std::make_pair(fit, sumErr);
}

void CSpectrumModel::setContinuumToInputSpc() {
  m_ContinuumFluxAxis.setSamplesVector(
      m_inputSpc->GetContinuumFluxAxis().GetSamplesVector());
}

void CSpectrumModel::setContinuumComponent(const std::string &component) {

  const Int32 spectrumSampleCount = m_inputSpc->GetSampleCount();

  if (component == "noContinuum") {
    // the continuum is set to zero and the observed spectrum is the spectrum
    // without continuum
    m_spcFluxAxisNoContinuum =
        m_inputSpc
            ->GetRawFluxAxis(); // m_inputSpc->GetWithoutContinuumFluxAxis();
    m_SpcFluxAxis = m_spcFluxAxisNoContinuum;
    m_SpectrumModel.SetFluxAxis(CSpectrumFluxAxis(spectrumSampleCount));
    m_ContinuumFluxAxis = CSpectrumFluxAxis(spectrumSampleCount);
  }
  if (component == "fromSpectrum") {
    // the continuum is set to the spectrum continuum and the observed
    // spectrum is the raw spectrum
    m_spcFluxAxisNoContinuum = m_inputSpc->GetWithoutContinuumFluxAxis();
    m_SpcFluxAxis = m_inputSpc->GetRawFluxAxis();
    m_ContinuumFluxAxis = m_inputSpc->GetContinuumFluxAxis();
    m_SpectrumModel.SetFluxAxis(m_ContinuumFluxAxis);
  }
  if (component == "tplFit" || component == "tplFitAuto" ||
      component == "powerLaw") {
    // the continuum is set to zero and the observed spectrum is the raw
    // spectrum
    m_SpcFluxAxis = m_inputSpc->GetRawFluxAxis();
    m_SpectrumModel.SetFluxAxis(CSpectrumFluxAxis(spectrumSampleCount));
    m_ContinuumFluxAxis = CSpectrumFluxAxis(spectrumSampleCount);
  }
}

void CSpectrumModel::setContinuumFromTplFit(Float64 alpha, Float64 tplAmp,
                                            const TFloat64List &polyCoeffs) {
  const CSpectrumSpectralAxis &spcSpectralAxis = m_inputSpc->GetSpectralAxis();
  const auto &SpcFluxAxis = m_SpcFluxAxis;
  for (Int32 k = 0; k < m_ContinuumFluxAxis.GetSamplesCount(); k++) {
    if (alpha == 1.0)
      m_ContinuumFluxAxis[k] = 0.;
    else
      m_ContinuumFluxAxis[k] =
          (1. - alpha) * m_observeGridContinuumFlux[k] * tplAmp;
    if (alpha != 0.0) {
      m_ContinuumFluxAxis[k] += alpha * m_inputSpc->GetContinuumFluxAxis()[k];
    }
    if (alpha != 1.0) {
      Float64 lbdaTerm = 1.0;
      for (Int32 kCoeff = 0; kCoeff < polyCoeffs.size(); kCoeff++) {
        m_ContinuumFluxAxis[k] += (1. - alpha) * polyCoeffs[kCoeff] * lbdaTerm;
        lbdaTerm *= spcSpectralAxis[k];
      }
    }
    m_spcFluxAxisNoContinuum[k] =
        SpcFluxAxis[k] - as_const(m_ContinuumFluxAxis)[k];
  }
}

/**
 * @brief CLineModelFitting::getFluxDirectIntegration
 * Integrates the flux (F-continuum) in a lbda range around the center observed
 * wavelength of the line. The wavelength range is defined by the instrument
 * resolution and a hardcoded nsigma factor
 * @param eIdx
 * @param subeIdx
 * @return
 */
std::pair<Float64, Float64> CSpectrumModel::getFluxDirectIntegration(
    const TInt32List &eIdx_list, const TInt32List &subeIdx_list,
    bool substract_abslinesmodel, const TFloat64Range &lambdaRange) const {

  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  Int32 nlines = eIdx_list.size();
  if (nlines != subeIdx_list.size())
    THROWG(ErrorCode::INTERNAL_ERROR, " index sizes do not match");
  TInt32RangeList indexRangeList = m_Elements->getlambdaIndexesUnderLines(
      eIdx_list, subeIdx_list, N_SIGMA_SUPPORT_DI, spectralAxis, lambdaRange,
      m_Redshift);

  if (!indexRangeList.size())
    THROWG(ErrorCode::INTERNAL_ERROR, "empty indexRanges ");

  const CSpectrumFluxAxis continuumFlux = getContinuumUnderLines(
      indexRangeList, eIdx_list, substract_abslinesmodel);
  // substarct continuum from spectrum flux
  const auto &SpcFluxAxis = getSpcFluxAxis();
  CSpectrumFluxAxis fluxMinusContinuum(SpcFluxAxis);
  for (const auto &r : indexRangeList)
    for (Int32 t = r.GetBegin(); t <= r.GetEnd(); t++)
      fluxMinusContinuum[t] -= continuumFlux[t];

  // estimate the integrated flux between obs. spectrum and continuum:
  // trapezoidal intg
  auto const &[sumFlux, sumErr] = CSpectrum::integrateFluxes_usingTrapez(
      m_inputSpc->GetSpectralAxis(), fluxMinusContinuum, indexRangeList);

  return std::make_pair(sumFlux, sumErr);
}

CSpectrumFluxAxis
CSpectrumModel::getContinuumUnderLines(const TInt32RangeList &indexRangeList,
                                       const TInt32List &eIdx_list,
                                       bool substract_abslinesmodel) const {
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  const auto &ContinuumFluxAxis = getContinuumFluxAxis();
  CSpectrumFluxAxis continuumFlux(ContinuumFluxAxis.GetSamplesCount());

  CSpectrumFluxAxis absLinesModelFlux;
  if (substract_abslinesmodel)
    absLinesModelFlux = getModel(
        eIdx_list, CLine::EType::nType_Absorption); // contains the continuum

  CSpectrumFluxAxis ampOffsetModelFlux;
  if (m_enableAmplitudeOffsets) {
    ampOffsetModelFlux = CSpectrumFluxAxis(spectralAxis.GetSamplesCount());
    m_Elements->addToSpectrumAmplitudeOffset(spectralAxis, ampOffsetModelFlux,
                                             eIdx_list);
  }

  // compute continuum
  for (const auto &r : indexRangeList)
    for (Int32 t = r.GetBegin(); t <= r.GetEnd(); t++) {
      continuumFlux[t] =
          substract_abslinesmodel ? absLinesModelFlux[t] : ContinuumFluxAxis[t];
      if (m_enableAmplitudeOffsets)
        continuumFlux[t] += ampOffsetModelFlux[t];
    }
  return continuumFlux;
}

/**
 * @brief CLineModelFitting::getLinesAboveSNR
 * Only considering a list of strong emission lines for now
 * @param snrcut
 * @return
 */
std::unordered_set<std::string>
CSpectrumModel::getLinesAboveSNR(const TFloat64Range &lambdaRange,
                                 Float64 snrcut) const {

  auto isElementInvalid = [this](Int32 eIdx, Int32 line_index) {
    return eIdx < 0 || line_index < 0 ||
           (*m_Elements)[eIdx]->IsOutsideLambdaRangeLine(line_index);
  };

  const auto lineList = {linetags::halpha_em,   linetags::oIIIa_em,
                         linetags::hbeta_em,    linetags::lya_em,
                         linetags::oII3726_em,  linetags::oII3729_em,
                         linetags::cIII1907_em, linetags::cIII1909_em};

  const TStringList lineNames = {"Ha",  "OIIIa", "Hb",   "Lya",
                                 "OII", "OII",   "CIII", "CIII"};

  std::unordered_set<std::string> str_above_cut;

  TInt32List eIdx_oii;
  TInt32List eIdx_ciii;
  TInt32List subeIdx_oii;
  TInt32List subeIdx_ciii;

  auto itLineNames = lineNames.begin();
  for (auto it = lineList.begin(), end = lineList.end(); it < end;
       ++it, ++itLineNames) {
    const auto &lineTag = *it;

    auto const it2 =
        std::find_if(m_RestLineList.cbegin(), m_RestLineList.cend(),
                     [lineTag](auto const &id_line) {
                       return id_line.second.GetName() == lineTag;
                     });
    if (it2 == m_RestLineList.cend()) // did not find line
      continue;

    auto const &[line_id, line] = *it2;

    bool isEmission = line.GetType() == CLine::EType::nType_Emission;
    if (!isEmission)
      continue;

    auto const &[eIdx, line_index] = m_Elements->findElementIndex(line_id);
    if (isElementInvalid(eIdx, line_index))
      continue;

    auto const &[mu, sigma] =
        (*m_Elements)[eIdx]->getObservedPositionAndLineWidth(m_Redshift,
                                                             line_index, false);
    Float64 fluxDI = NAN;
    Float64 snrDI = NAN;
    TInt32List eIdx_line(1, eIdx);
    TInt32List subeIdx_line(1, line_index);

    bool opt_cont_substract_abslinesmodel = isEmission;
    if (lineTag == linetags::oII3726_em || lineTag == linetags::oII3729_em) {
      opt_cont_substract_abslinesmodel = false;
      eIdx_oii.push_back(eIdx);
      subeIdx_oii.push_back(line_index);
      eIdx_line = eIdx_oii;
      subeIdx_line = subeIdx_oii;
    };
    if (lineTag == linetags::cIII1907_em || lineTag == linetags::cIII1909_em) {
      opt_cont_substract_abslinesmodel = false;
      eIdx_ciii.push_back(eIdx);
      subeIdx_ciii.push_back(line_index);
      eIdx_line = eIdx_ciii;
      subeIdx_line = subeIdx_ciii;
    }
    auto [sumFlux, sumErr] = getFluxDirectIntegration(
        eIdx_line, subeIdx_line, opt_cont_substract_abslinesmodel, lambdaRange);

    fluxDI = sumFlux;
    if (sumErr > 0.)
      snrDI = std::abs(fluxDI) / sqrt(sumErr);

    if (snrDI > snrcut)
      str_above_cut.insert(
          *itLineNames); // this ensure unicity  (for oII and cIII)
  }
  return str_above_cut;
}

// TODO should be renamed
CSpectrumFluxAxis CSpectrumModel::getModel(const TInt32List &eIdx_list,
                                           CLine::EType lineTypeFilter) const {

  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  CSpectrumFluxAxis modelfluxAxis(spectralAxis.GetSamplesCount());

  Int32 nElements = m_Elements->size();
  for (Int32 eIdx : eIdx_list) {
    const auto &elt = (*m_Elements)[eIdx];
    elt->initSpectrumModel(modelfluxAxis, getContinuumFluxAxis());

    auto const lineType = elt->getElementParam()->GetElementType();
    if (lineTypeFilter == CLine::EType::nType_All || lineTypeFilter == lineType)
      elt->addToSpectrumModel(spectralAxis, modelfluxAxis,
                              getContinuumFluxAxis(), m_Redshift);
  }

  return modelfluxAxis;
}

/**
 * Apply the template continuum by interpolating the grid as define in Init
 * Continuum
 */
Int32 CSpectrumModel::ApplyContinuumTplOnGrid(
    const std::shared_ptr<const CTemplate> &tpl, Float64 zcontinuum) {
  m_fitContinuum->tplName = tpl->GetName();
  Int32 n = tpl->GetSampleCount();

  Int32 idxDust = -1;
  if (m_fitContinuum->ebmvCoef > 0.) {
    if (tpl->CalzettiInitFailed()) {
      THROWG(ErrorCode::INTERNAL_ERROR,
             "  no calzetti calib. file in template");
    }
    idxDust =
        tpl->m_ismCorrectionCalzetti->GetEbmvIndex(m_fitContinuum->ebmvCoef);
  }
  const CSpectrumSpectralAxis &tplSpectralAxis = tpl->GetSpectralAxis();
  TFloat64Range range(tplSpectralAxis[0], tplSpectralAxis[n - 1]);

  std::string inter_opt = "spline";
  tpl->setRebinInterpMethod(inter_opt);
  Float64 overlapThreshold = 1., amplitude = 1.;
  std::shared_ptr<CModelSpectrumResult> spcmodel =
      std::make_shared<CModelSpectrumResult>();
  m_photValues =
      (std::dynamic_pointer_cast<COperatorTemplateFittingBase>(
           m_continuumFittingOperator))
          ->ComputeSpectrumModel(tpl, zcontinuum, m_fitContinuum->ebmvCoef,
                                 m_fitContinuum->meiksinIdx, amplitude,
                                 overlapThreshold, m_spcIndex, spcmodel);
  if (spcmodel == nullptr)
    THROWG(ErrorCode::INTERNAL_ERROR, "Couldnt compute spectrum model");

  // m_observeGridContinuumFlux should be a CSpectrumFluxAxis not
  // AxisSampleList
  m_observeGridContinuumFlux =
      std::move((*spcmodel).ModelFlux.at(m_inputSpc->getObsID()));

  return 0;
}

Int32 CSpectrumModel::ApplyContinuumPowerLawOnGrid(
    std::shared_ptr<CContinuumModelSolution> const &continuum) {
  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc->GetSpectralAxis();

  std::shared_ptr<CModelSpectrumResult> spcmodel =
      std::make_shared<CModelSpectrumResult>();

  (std::dynamic_pointer_cast<COperatorPowerLaw>(m_continuumFittingOperator))
      ->ComputeSpectrumModel(continuum, m_spcIndex, spcmodel);

  if (spcmodel == nullptr)
    THROWG(ErrorCode::INTERNAL_ERROR, "Couldnt compute spectrum model");

  m_observeGridContinuumFlux =
      std::move((*spcmodel).ModelFlux.at(m_inputSpc->getObsID()));
  for (Int32 k = 0; k < m_ContinuumFluxAxis.GetSamplesCount(); k++) {
    m_ContinuumFluxAxis[k] = m_observeGridContinuumFlux[k];
    m_spcFluxAxisNoContinuum[k] =
        m_SpcFluxAxis[k] - as_const(m_ContinuumFluxAxis)[k];
  }
  return 0;
}

void CSpectrumModel::initObserveGridContinuumFlux(Int32 size) {
  m_observeGridContinuumFlux.resize(size);
}
