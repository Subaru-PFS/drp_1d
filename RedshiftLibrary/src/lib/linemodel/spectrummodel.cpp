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
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;
using namespace std;

CSpectrumModel::CSpectrumModel(
    const std::shared_ptr<CLineModelElementList> &elements,
    const std::shared_ptr<const CSpectrum> &spc, const CLineMap &restLineList,
    const std::shared_ptr<CTplModelSolution> &tfv,
    const std::shared_ptr<COperatorTemplateFittingBase> &TFOperator,
    Int32 spcIndex)
    : m_Elements(elements), m_inputSpc(spc), m_SpectrumModel(*(spc)),
      m_RestLineList(restLineList), m_fitContinuum(tfv),
      m_templateFittingOperator(TFOperator), m_spcIndex(spcIndex) {
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
    auto const lineType = (*m_Elements)[iElts]->GetElementType();
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

  //
  TFloat64List fluxAndContinuum(spectralAxis.GetSamplesCount(), 0.0);
  if (lineTypeFilter == CLine::EType::nType_Emission)
    refreshModel(CLine::EType::nType_Absorption);
  else if (lineTypeFilter == CLine::EType::nType_Absorption)
    refreshModel(CLine::EType::nType_Emission);

  fluxAndContinuum = fluxAxis.GetSamplesVector();
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

/**
 * @brief GetContinuumError
 * Estimate the error on the continuum in a given window around the line center.
 * 1. calculate the observed spectrum flux with lines subtracted (fitted line
 * model)
 * 2. estimate the std in a window corresponding to the sliding median window
 * size (default/hardcoded=150A)
 * @param subeIdx
 * @param spectralAxis
 * @param spectrumfluxAxis
 * @return -1: zero samples found for error estimation, NAN: not enough samples
 * found for error estimation
 */
std::pair<Float64, Int32>
CSpectrumModel::getContinuumQuadraticError(Int32 eIdx, Int32 line_id) {

  const CSpectrum noLinesSpectrum = GetObservedSpectrumWithLinesRemoved();
  const CSpectrumFluxAxis &noLinesFluxAxis = noLinesSpectrum.GetFluxAxis();
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  const TFloat64Range lambdaRange =
      spectralAxis.GetLambdaRange(); // using the full wavelength range for this
                                     // error estimation
  const auto &ContinuumFluxAxis = m_ContinuumFluxAxis;
  Float64 winsizeAngstrom = 150.;

  Float64 mu = (*m_Elements)[eIdx]->GetObservedPosition(line_id, m_Redshift);
  TInt32Range indexRange = CLineModelElement::EstimateIndexRange(
      spectralAxis, mu, lambdaRange, winsizeAngstrom);

  // estimate sum square error between continuum and nolines-spectrum
  Float64 sum = 0.0;
  Int32 nsum = 0;
  for (Int32 t = indexRange.GetBegin(); t < indexRange.GetEnd(); t++) {
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
CSpectrumModel::getModelQuadraticErrorUnderElement(Int32 eltId,
                                                   bool with_continuum) const {
  // before elementlistcutting this variable was
  // CElementList::m_ErrorNoContinuum, a reference initialized twice in
  // CElementList constructor, first init to m_spcFluxAxisNoContinuum.GetError()
  // and after to spectrumFluxAxis.GetError
  const CSpectrumNoiseAxis &errorNoContinuum =
      m_SpectrumModel.GetFluxAxis().GetError();
  const CSpectrumFluxAxis &fluxRef =
      with_continuum ? getSpcFluxAxis() : getSpcFluxAxisNoContinuum();

  if (eltId < 0)
    return std::make_pair(NAN, NAN);

  if ((*m_Elements)[eltId]->IsOutsideLambdaRange())
    return std::make_pair(0.0, 0.0);

  Int32 numDevs = 0;
  Float64 fit = 0.0;
  const TAxisSampleList &Ymodel =
      m_SpectrumModel.GetFluxAxis().GetSamplesVector();
  const TAxisSampleList &Yspc = fluxRef.GetSamplesVector();
  Float64 diff = 0.0;
  Float64 sumErr = 0.0;
  TInt32RangeList support = (*m_Elements)[eltId]->getSupport();

  Float64 w = 0.0;
  for (const auto &s : support) {
    for (Int32 j = s.GetBegin(), e = s.GetEnd(); j <= e; j++) {
      numDevs++;
      diff = (Yspc[j] - Ymodel[j]);
      w = 1.0 / (errorNoContinuum[j] * errorNoContinuum[j]);
      fit += (diff * diff) * w;
      sumErr += w;
    }
  }
  return std::make_pair(fit, sumErr);
}

void CSpectrumModel::setContinuumToInputSpc() {
  m_ContinuumFluxAxis.setSamplesVector(
      m_inputSpc->GetContinuumFluxAxis().GetSamplesVector());
}

void CSpectrumModel::setContinuumComponent(const std::string &component) {

  const Int32 spectrumSampleCount = m_inputSpc->GetSampleCount();

  if (component == "nocontinuum") {
    // the continuum is set to zero and the observed spectrum is the spectrum
    // without continuum
    m_spcFluxAxisNoContinuum =
        m_inputSpc
            ->GetRawFluxAxis(); // m_inputSpc->GetWithoutContinuumFluxAxis();
    m_SpcFluxAxis = m_spcFluxAxisNoContinuum;
    m_SpectrumModel.SetFluxAxis(CSpectrumFluxAxis(spectrumSampleCount));
    m_ContinuumFluxAxis = CSpectrumFluxAxis(spectrumSampleCount);
  }
  if (component == "fromspectrum") {
    // the continuum is set to the spectrum continuum and the observed
    // spectrum is the raw spectrum
    m_spcFluxAxisNoContinuum = m_inputSpc->GetWithoutContinuumFluxAxis();
    m_SpcFluxAxis = m_inputSpc->GetRawFluxAxis();
    m_ContinuumFluxAxis = m_inputSpc->GetContinuumFluxAxis();
    m_SpectrumModel.SetFluxAxis(m_ContinuumFluxAxis);
  }
  if (component == "tplfit" || component == "tplfitauto") {
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
    THROWG(INTERNAL_ERROR, " index sizes do not match");
  TInt32RangeList indexRangeList = m_Elements->getlambdaIndexesUnderLines(
      eIdx_list, subeIdx_list, N_SIGMA_SUPPORT_DI, spectralAxis, lambdaRange,
      m_Redshift);

  if (!indexRangeList.size())
    THROWG(INTERNAL_ERROR, "empty indexRanges ");

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
  Float64 sumFlux = 0.0;
  Float64 sumErr = 0.0;
  integrateFluxes_usingTrapez(fluxMinusContinuum, indexRangeList, sumFlux,
                              sumErr);

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
           (*m_Elements)[eIdx]->IsOutsideLambdaRange(line_index);
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

    Float64 cont = (*m_Elements)[eIdx]->GetContinuumAtCenterProfile(
        line_index, m_inputSpc->GetSpectralAxis(), m_Redshift,
        getContinuumFluxAxis(), m_enableAmplitudeOffsets);
    Float64 mu = NAN;
    Float64 sigma = NAN;
    (*m_Elements)[eIdx]->getObservedPositionAndLineWidth(line_index, m_Redshift,
                                                         mu, sigma, false);
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

    auto const lineType = elt->GetElementType();
    if (lineTypeFilter == CLine::EType::nType_All || lineTypeFilter == lineType)
      elt->addToSpectrumModel(spectralAxis, modelfluxAxis,
                              getContinuumFluxAxis(), m_Redshift);
  }

  return modelfluxAxis;
}

/**
 * @brief Construct a new clinemodelfitting::integratefluxes usingtrapez
 * object
 *
 * @param fluxMinusContinuum
 * @param indexRangeList
 * @param sumFlux
 * @param sumErr
 */
void CSpectrumModel::integrateFluxes_usingTrapez(
    const CSpectrumFluxAxis &fluxMinusContinuum,
    const TInt32RangeList &indexRangeList, Float64 &sumFlux,
    Float64 &sumErr) const {

  sumFlux = 0.0;
  sumErr = 0.0;

  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc->GetSpectralAxis();
  if (spectralAxis.GetSamplesCount() < 2)
    THROWG(INTERNAL_ERROR, "Not enough samples in spectral axis");

  const auto &ErrorNoContinuum = m_inputSpc->GetFluxAxis().GetError();
  for (auto &r : indexRangeList) {
    for (Int32 t = r.GetBegin(), e = r.GetEnd(); t < e; t++) {
      Float64 trapweight = (spectralAxis[t + 1] - spectralAxis[t]) * 0.5;
      sumFlux +=
          trapweight * (fluxMinusContinuum[t + 1] + fluxMinusContinuum[t]);

      Float64 ea = ErrorNoContinuum[t] * ErrorNoContinuum[t];
      Float64 eb = ErrorNoContinuum[t + 1] * ErrorNoContinuum[t + 1];
      sumErr += trapweight * trapweight * (eb + ea);
    }
  }
  return;
}

/**
 * Apply the template continuum by interpolating the grid as define in Init
 * Continuum
 */
Int32 CSpectrumModel::ApplyContinuumOnGrid(
    const std::shared_ptr<const CTemplate> &tpl, Float64 zcontinuum) {
  m_fitContinuum->tplName = tpl->GetName();
  Int32 n = tpl->GetSampleCount();

  Int32 idxDust = -1;
  if (m_fitContinuum->tplEbmvCoeff > 0.) {
    if (tpl->CalzettiInitFailed()) {
      THROWG(INTERNAL_ERROR, "  no calzetti calib. file in template");
    }
    idxDust = tpl->m_ismCorrectionCalzetti->GetEbmvIndex(
        m_fitContinuum->tplEbmvCoeff);
  }
  const CSpectrumSpectralAxis &tplSpectralAxis = tpl->GetSpectralAxis();
  TFloat64Range range(tplSpectralAxis[0], tplSpectralAxis[n - 1]);

  std::string inter_opt = "spline";
  tpl->setRebinInterpMethod(inter_opt);
  Float64 overlapThreshold = 1., amplitude = 1.;
  std::shared_ptr<CModelSpectrumResult> spcmodel =
      std::make_shared<CModelSpectrumResult>();
  m_photValues = m_templateFittingOperator->ComputeSpectrumModel(
      tpl, zcontinuum, m_fitContinuum->tplEbmvCoeff,
      m_fitContinuum->tplMeiksinIdx, amplitude, overlapThreshold, m_spcIndex,
      spcmodel);
  if (spcmodel == nullptr)
    THROWG(INTERNAL_ERROR, "Couldnt compute spectrum model");

  // m_observeGridContinuumFlux should be a CSpectrumFluxAxis not
  // AxisSampleList
  m_observeGridContinuumFlux =
      std::move((*spcmodel).ModelFlux.at(m_inputSpc->getObsID()));

  return 0;
}

void CSpectrumModel::initObserveGridContinuumFlux(Int32 size) {
  m_observeGridContinuumFlux.resize(size);
}
