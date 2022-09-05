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

CSpectrumModel::CSpectrumModel(CLineModelElementList &elements,
                               std::shared_ptr<const CSpectrum> spc,
                               const CLineCatalog::TLineVector &restLineList)
    : m_Elements(elements), m_inputSpc(spc), m_spcCorrectedUnderLines(*(spc)),
      m_SpectrumModel(*(spc)), m_RestLineList(restLineList) {
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

/**
 * \brief Init the whole spectrum model with continuum.
 **/
void CSpectrumModel::reinitModel() {
  CSpectrumFluxAxis modelFluxAxis = m_SpectrumModel.GetFluxAxis();
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    m_Elements[iElts]->initSpectrumModel(modelFluxAxis, m_ContinuumFluxAxis);
  }

  m_SpectrumModel.SetFluxAxis(std::move(modelFluxAxis));
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
  Int32 iElts;
  CSpectrumFluxAxis modelFluxAxis = m_SpectrumModel.GetFluxAxis();
  // init spectrum model with continuum
  for (Int32 i = 0; i < filterEltsIdx.size(); i++) {
    iElts = filterEltsIdx[i];
    m_Elements[iElts]->initSpectrumModel(modelFluxAxis, m_ContinuumFluxAxis,
                                         lineIdx);
  }
  m_SpectrumModel.SetFluxAxis(std::move(modelFluxAxis));
}

void CSpectrumModel::refreshModel(Int32 lineTypeFilter) {

  reinitModel();
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  CSpectrumFluxAxis modelFluxAxis = m_SpectrumModel.GetFluxAxis();

  if (m_enableAmplitudeOffsets) {
    // add amplitude offsets
    m_Elements.addToSpectrumAmplitudeOffset(m_SpectrumModel.GetSpectralAxis(),
                                            modelFluxAxis);
    // addToSpectrumAmplitudeOffset(contFluxAxisWithAmpOffset);
  }

  // create spectrum model
  Int32 nElements = m_Elements.size();
  for (Int32 iElts = 0; iElts < nElements; iElts++) {
    Int32 lineType = m_Elements[iElts]->m_Lines[0].GetType();
    if (lineTypeFilter == -1 || lineTypeFilter == lineType) {
      m_Elements[iElts]->addToSpectrumModel(spectralAxis, modelFluxAxis,
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
  Int32 iElts;
  for (Int32 i = 0; i < filterEltsIdx.size(); i++) {
    iElts = filterEltsIdx[i];
    m_Elements[iElts]->addToSpectrumModel(
        spectralAxis, modelFluxAxis, m_ContinuumFluxAxis, m_Redshift, lineIdx);
  }
  m_SpectrumModel.SetFluxAxis(std::move(modelFluxAxis));
}

/**
 * \brief refreshing all the grid and Adds a new model to each m_Elements
 *entry . Calls iterate on model flux to set it equal to continuum. For each
 *entry in m_Elements, addToSpectrumModel using the reinitModel output as
 *arguments.
 **/
void CSpectrumModel::refreshModelInitAllGrid() {

  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  CSpectrumFluxAxis modelFluxAxis = m_ContinuumFluxAxis;

  // create spectrum model
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    m_Elements[iElts]->addToSpectrumModel(spectralAxis, modelFluxAxis,
                                          m_ContinuumFluxAxis, m_Redshift);
  }

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
  TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();
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
    for (Int32 t = 0; t < spectralAxis.GetSamplesCount(); t++) {
      if (filterOnlyAbs) {
        if (spcmodel4linefittingFluxAxis[t] < 0.0) {
          spcmodel4linefittingFluxAxis[t] *= opt_enhance_lines;
        }
      } else {
        spcmodel4linefittingFluxAxis[t] *= opt_enhance_lines;
      }
    }
  }

  // subtract the lines component
  CSpectrumFluxAxis fluxAxisNothingUnderLines = m_SpcFluxAxis;
  for (Int32 t = 0; t < spectralAxis.GetSamplesCount(); t++) {
    fluxAxisNothingUnderLines[t] -= as_const(spcmodel4linefittingFluxAxis)[t];
  }

  // evaluate contiuum
  CSpectrum spcCorrectedUnderLines(*(m_inputSpc));
  spcCorrectedUnderLines.SetFluxAxis(std::move(fluxAxisNothingUnderLines));
  m_ContinuumFluxAxis = spcCorrectedUnderLines.GetContinuumFluxAxis();
}

/**
 * \brief Returns a pointer to a spectrum containing the observed spectrum with
 *the fitted lines subtracted
 **/
const CSpectrum &
CSpectrumModel::GetObservedSpectrumWithLinesRemoved(Int32 lineTypeFilter) {
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  const CSpectrumFluxAxis &fluxAxis = m_SpectrumModel.GetFluxAxis();

  //
  TFloat64List fluxAndContinuum(spectralAxis.GetSamplesCount(), 0.0);
  if (lineTypeFilter == CLine::nType_Emission) {
    refreshModel(CLine::nType_Absorption);
    fluxAndContinuum = fluxAxis.GetSamplesVector();
    refreshModel(CLine::nType_Emission);
  } else if (lineTypeFilter == CLine::nType_Absorption) {
    refreshModel(CLine::nType_Emission);
    fluxAndContinuum = fluxAxis.GetSamplesVector();
    refreshModel(CLine::nType_Absorption);
  }

  CSpectrumFluxAxis fluxAxisNothingUnderLines(
      m_spcCorrectedUnderLines.GetSampleCount());
  TAxisSampleList &Y = fluxAxisNothingUnderLines.GetSamplesVector();
  const auto &SpcFluxAxis = m_SpcFluxAxis;
  const auto &ContinuumFluxAxis = m_ContinuumFluxAxis;
  for (Int32 t = 0; t < spectralAxis.GetSamplesCount(); t++) {
    Y[t] = SpcFluxAxis[t] - fluxAxis[t] + ContinuumFluxAxis[t];
  }

  bool enableSmoothBlendContinuUnderLines = true;
  if (enableSmoothBlendContinuUnderLines) {
    Float64 alphaMax = 0.9; // alpha blend = 0: only lineSubtractedFlux,
                            // alpha=1: only continuum

    TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();
    TInt32List nonZeroValidEltsIdx;
    for (Int32 j = 0; j < validEltsIdx.size(); j++) {
      Int32 lineType = m_Elements[validEltsIdx[j]]->m_Lines[0].GetType();
      if (!(lineTypeFilter == -1 || lineTypeFilter == lineType)) {
        continue;
      }

      if (m_Elements[validEltsIdx[j]]->GetElementAmplitude() > 0.0) {
        nonZeroValidEltsIdx.push_back(validEltsIdx[j]);
      }
    }
    if (nonZeroValidEltsIdx.size() > 0) {
      TInt32List supportIdxes =
          m_Elements.getSupportIndexes(nonZeroValidEltsIdx);
      if (supportIdxes.size() > 0) {
        for (Int32 i = 1; i < supportIdxes.size(); i++) {
          Float64 weighting = GetWeightingAnyLineCenterProximity(
              supportIdxes[i], nonZeroValidEltsIdx);
          Float64 alpha = alphaMax * weighting;
          Y[supportIdxes[i]] = (1. - alpha) * Y[supportIdxes[i]] +
                               alpha * (fluxAndContinuum[supportIdxes[i]]);
        }
      }
    }
  }

  m_spcCorrectedUnderLines.SetFluxAxis(std::move(fluxAxisNothingUnderLines));

  return m_spcCorrectedUnderLines;
}

Float64 CSpectrumModel::GetWeightingAnyLineCenterProximity(
    Int32 sampleIndex, const TInt32List &EltsIdx) const {
  Float64 maxWeight = 0.0;
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  Float64 currentLbda = spectralAxis[sampleIndex];

  for (Int32 i = 0; i < EltsIdx.size(); i++) {
    Int32 iElts = EltsIdx[i];

    // getTheoreticalSupport reads from the 2 class variables
    TInt32RangeList s = m_Elements[iElts]->getTheoreticalSupport();
    for (Int32 iS = 0; iS < s.size(); iS++) {
      Float64 weight = 0;
      TInt32Range range = s[iS];
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

      if (maxWeight < weight) {
        maxWeight = weight;
      }
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
Float64 CSpectrumModel::GetContinuumError(Int32 eIdx, Int32 subeIdx) {
  Int32 nMinValueForErrorEstimation = 10;

  const CSpectrum &noLinesSpectrum = GetObservedSpectrumWithLinesRemoved();
  const CSpectrumFluxAxis &noLinesFluxAxis = noLinesSpectrum.GetFluxAxis();
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  const TFloat64Range lambdaRange =
      spectralAxis.GetLambdaRange(); // using the full wavelength range for this
                                     // error estimation
  const auto &ContinuumFluxAxis = m_ContinuumFluxAxis;
  Float64 winsizeAngstrom = 150.;

  Float64 mu = m_Elements[eIdx]->GetObservedPosition(subeIdx, m_Redshift);
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

  Float64 error = NAN;
  if (nsum >= nMinValueForErrorEstimation) {
    error = sqrt(sum / nsum);
  }

  return error;
}

/**
 * \brief Returns the error of the support for subelements under the element
 *with the argument eltId as index. Accumulate "fit", the squared difference
 *between model and spectrum, divided by the square of the m_ErrorNoContinuum
 *value. Accumulate "sumErr" 1 / square of the m_ErrorNoContinuum value. return
 *the square root of fit / sumErr.
 **/
Float64 CSpectrumModel::getModelErrorUnderElement(Int32 eltId) const {
  // before elementlistcutting this variable was
  // CElementList::m_ErrorNoContinuum, a reference initialized twice in
  // CElementList constructor, first init to m_spcFluxAxisNoContinuum.GetError()
  // and after to spectrumFluxAxis.GetError
  const CSpectrumNoiseAxis &errorNoContinuum =
      m_SpectrumModel.GetFluxAxis().GetError();

  if (eltId < 0) {
    return -1.0;
  }

  Int32 numDevs = 0;
  Float64 fit = 0.0;
  const Float64 *Ymodel = m_SpectrumModel.GetFluxAxis().GetSamples();
  const Float64 *Yspc = m_inputSpc->GetFluxAxis().GetSamples();
  Float64 diff = 0.0;

  Float64 sumErr = 0.0;

  TInt32RangeList support;
  Int32 iElts = eltId;
  {
    if (m_Elements[iElts]->IsOutsideLambdaRange()) {
      return 0.0;
    }
    TInt32RangeList s = m_Elements[iElts]->getSupport();
    for (Int32 iS = 0; iS < s.size(); iS++) {
      support.push_back(s[iS]);
    }
  }

  Float64 w = 0.0;
  for (Int32 iS = 0; iS < support.size(); iS++) {
    for (Int32 j = support[iS].GetBegin(); j < support[iS].GetEnd(); j++) {
      numDevs++;
      diff = (Yspc[j] - Ymodel[j]);
      w = 1.0 / (errorNoContinuum[j] * errorNoContinuum[j]);
      fit += (diff * diff) * w;
      sumErr += w;
    }
  }
  return sqrt(fit / sumErr);
}

/**
 * \brief This function prepares the continuum for use in the fit with the
 *line elements. Rebin with PFG buffer Find and apply amplitude factor from
 *previously fitted tpl
 **/
void CSpectrumModel::PrepareContinuum() {
  const CSpectrumSpectralAxis &targetSpectralAxis =
      m_inputSpc->GetSpectralAxis();
  TAxisSampleList &Yrebin = m_ContinuumFluxAxis.GetSamplesVector();

  if (m_ContinuumComponent == "nocontinuum")
    return;

  if (!isContinuumComponentTplfitxx()) {
    Yrebin = m_inputSpc->GetContinuumFluxAxis().GetSamplesVector();
  }
}

void CSpectrumModel::SetContinuumComponent(std::string component) {
  m_ContinuumComponent = component;

  const Int32 spectrumSampleCount = m_inputSpc->GetSampleCount();

  if (m_ContinuumComponent == "nocontinuum") {
    // the continuum is set to zero and the observed spectrum is the spectrum
    // without continuum
    m_spcFluxAxisNoContinuum =
        m_inputSpc
            ->GetRawFluxAxis(); // m_inputSpc->GetWithoutContinuumFluxAxis();
    m_SpcFluxAxis = m_spcFluxAxisNoContinuum;
    m_SpectrumModel.SetFluxAxis(CSpectrumFluxAxis(spectrumSampleCount));
  }
  if (m_ContinuumComponent == "fromspectrum") {
    // the continuum is set to the spectrum continuum and the observed
    // spectrum is the raw spectrum
    m_spcFluxAxisNoContinuum = m_inputSpc->GetWithoutContinuumFluxAxis();
    m_SpcFluxAxis = m_inputSpc->GetRawFluxAxis();
    m_ContinuumFluxAxis = m_inputSpc->GetContinuumFluxAxis();
    m_SpectrumModel.SetFluxAxis(m_ContinuumFluxAxis);
  }
  if (isContinuumComponentTplfitxx()) {
    // the continuum is set to zero and the observed spectrum is the raw
    // spectrum
    m_SpcFluxAxis = m_inputSpc->GetRawFluxAxis();
    m_SpectrumModel.SetFluxAxis(CSpectrumFluxAxis(spectrumSampleCount));
  }
}

void CSpectrumModel::setContinuumFromTplFit(
    Float64 alpha, Float64 tplAmp, const TFloat64List &polyCoeffs,
    const TAxisSampleList &observeGridContinuumFlux) {
  const CSpectrumSpectralAxis &spcSpectralAxis = m_inputSpc->GetSpectralAxis();
  const auto &SpcFluxAxis = m_SpcFluxAxis;
  for (Int32 k = 0; k < m_ContinuumFluxAxis.GetSamplesCount(); k++) {
    if (alpha == 1.0)
      m_ContinuumFluxAxis[k] = 0.;
    else
      m_ContinuumFluxAxis[k] =
          (1. - alpha) * observeGridContinuumFlux[k] * tplAmp;
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
void CSpectrumModel::getFluxDirectIntegration(
    const TInt32List &eIdx_list, const TInt32List &subeIdx_list,
    bool substract_abslinesmodel, Float64 &fluxdi, Float64 &snrdi,
    const TFloat64Range &lambdaRange) const {

  fluxdi = NAN;
  snrdi = NAN;
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  Int32 nlines = eIdx_list.size();
  if (nlines != subeIdx_list.size())
    THROWG(INTERNAL_ERROR, " index sizes do not match");
  TInt32RangeList indexRangeList = m_Elements.getlambdaIndexesUnderLines(
      eIdx_list, subeIdx_list, N_SIGMA_SUPPORT_DI, spectralAxis, lambdaRange,
      m_Redshift);

  if (!indexRangeList.size())
    THROWG(INTERNAL_ERROR, "empty indexRanges ");

  const auto &ContinuumFluxAxis = getContinuumFluxAxis();
  CSpectrumFluxAxis continuumFlux(ContinuumFluxAxis.GetSamplesCount());
  CSpectrumFluxAxis absLinesModelFlux;
  if (substract_abslinesmodel)
    absLinesModelFlux =
        getModel(CLine::nType_Absorption); // contains the continuum
                                           // and abs lines model;

  // polynome are shared between overlapping CElements
  //  find polynome belonging to this CElement
  // TODO: correct this for more than one CElement
  TPolynomCoeffs polynom_coeffs{0., 0., 0.};
  if (m_enableAmplitudeOffsets)
    polynom_coeffs = m_Elements.getPolynomCoeffs(eIdx_list[0]);

  // compute continuum
  for (const auto &r : indexRangeList)
    for (Int32 t = r.GetBegin(); t <= r.GetEnd(); t++) {
      continuumFlux[t] =
          substract_abslinesmodel ? absLinesModelFlux[t] : ContinuumFluxAxis[t];
      if (m_enableAmplitudeOffsets)
        continuumFlux[t] +=
            polynom_coeffs.x0 + polynom_coeffs.x1 * spectralAxis[t] +
            polynom_coeffs.x2 * spectralAxis[t] * spectralAxis[t];
    }
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

  if (sumErr <= 0)
    return;

  fluxdi = sumFlux;
  snrdi = std::abs(fluxdi) / sqrt(sumErr);
  return;
}

/**
 * @brief CLineModelFitting::getLinesAboveSNR
 * Only considering a list of strong emission lines for now
 * @param snrcut
 * @return
 */
TStringList CSpectrumModel::getLinesAboveSNR(const TFloat64Range &lambdaRange,
                                             Float64 snrcut) const {
  TInt32List eIdx_oii;
  TInt32List subeIdx_oii;
  TInt32List eIdx_ciii;
  TInt32List subeIdx_ciii;

  Float64 snr_ha = NAN;
  Float64 snr_oii = NAN;
  Float64 snr_oiiia = NAN;
  // snr_oiiib = NAN;
  Float64 snr_hb = NAN;
  Float64 snr_ciii = NAN;
  Float64 snr_lya = NAN;

  for (Int32 iRestLine = 0; iRestLine < m_RestLineList.size(); iRestLine++) {
    Int32 subeIdx = undefIdx;
    Int32 eIdx = m_Elements.findElementIndex(iRestLine, subeIdx);
    if (eIdx < 0 || subeIdx < 0 ||
        m_Elements[eIdx]->IsOutsideLambdaRange(subeIdx))
      continue;

    TPolynomCoeffs polynom_coeffs = m_Elements.getPolynomCoeffs(eIdx);
    Float64 cont = m_Elements[eIdx]->GetContinuumAtCenterProfile(
        subeIdx, m_inputSpc->GetSpectralAxis(), m_Redshift,
        getContinuumFluxAxis(), polynom_coeffs);
    Float64 mu = NAN;
    Float64 sigma = NAN;
    m_Elements[eIdx]->getObservedPositionAndLineWidth(subeIdx, m_Redshift, mu,
                                                      sigma, false);

    Float64 flux = NAN;
    Float64 fluxError = NAN;
    Float64 fluxDI = NAN;
    Float64 snrDI = NAN;
    TInt32List eIdx_line(1, eIdx);
    TInt32List subeIdx_line(1, subeIdx);
    bool isEmission =
        m_RestLineList[iRestLine].GetType() == CLine::nType_Emission;
    bool opt_cont_substract_abslinesmodel = isEmission;
    getFluxDirectIntegration(eIdx_line, subeIdx_line,
                             opt_cont_substract_abslinesmodel, fluxDI, snrDI,
                             lambdaRange);

    if (m_RestLineList[iRestLine].GetName() == linetags::halpha_em &&
        isEmission)
      snr_ha = snrDI;

    if (m_RestLineList[iRestLine].GetName() == linetags::oIIIa_em && isEmission)
      snr_oiiia = snrDI;

    if (m_RestLineList[iRestLine].GetName() == linetags::hbeta_em && isEmission)
      snr_hb = snrDI;

    if (m_RestLineList[iRestLine].GetName() == linetags::lya_em && isEmission)
      snr_lya = snrDI;

    if ((m_RestLineList[iRestLine].GetName() == linetags::oII3726_em ||
         m_RestLineList[iRestLine].GetName() == linetags::oII3729_em) &&
        isEmission) {
      // here we only cover the fluxDI case.
      eIdx_oii.push_back(eIdx);
      subeIdx_oii.push_back(subeIdx);
      fluxDI = NAN;
      snrDI = NAN;
      opt_cont_substract_abslinesmodel = false;
      getFluxDirectIntegration(eIdx_oii, subeIdx_oii,
                               opt_cont_substract_abslinesmodel, fluxDI, snrDI,
                               lambdaRange);

      snr_oii = snrDI;
    }
    if ((m_RestLineList[iRestLine].GetName() == linetags::cIII1907_em ||
         m_RestLineList[iRestLine].GetName() == linetags::cIII1909_em) &&
        isEmission) {
      // here we only cover the fluxDI case.
      eIdx_ciii.push_back(eIdx);
      subeIdx_ciii.push_back(subeIdx);
      fluxDI = NAN;
      snrDI = NAN;
      opt_cont_substract_abslinesmodel = false;
      getFluxDirectIntegration(eIdx_ciii, subeIdx_ciii,
                               opt_cont_substract_abslinesmodel, fluxDI, snrDI,
                               lambdaRange);

      snr_ciii = snrDI;
    }
  }

  // sum snr ht cut
  Int32 nhtcut = 0;
  TStringList str_above_cut;
  if (snr_ha > snrcut) {
    nhtcut++;
    str_above_cut.push_back("Ha");
  }
  if (snr_oiiia > snrcut) {
    nhtcut++;
    str_above_cut.push_back("OIIIa");
  }
  if (snr_hb > snrcut) {
    nhtcut++;
    str_above_cut.push_back("Hb");
  }
  if (snr_oii > snrcut) {
    nhtcut++;
    str_above_cut.push_back("OII");
  }
  if (snr_lya > snrcut) {
    nhtcut++;
    str_above_cut.push_back("Lya");
  }
  if (snr_ciii > snrcut) {
    nhtcut++;
    str_above_cut.push_back("CIII");
  }

  return str_above_cut;
}

// TODO should be renamed
CSpectrumFluxAxis CSpectrumModel::getModel(Int32 lineTypeFilter) const {

  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  CSpectrumFluxAxis modelfluxAxis(spectralAxis.GetSamplesCount());

  Int32 nElements = m_Elements.size();
  for (Int32 iElts = 0; iElts < nElements; iElts++) {
    m_Elements[iElts]->initSpectrumModel(modelfluxAxis, getContinuumFluxAxis());
  }

  for (Int32 iElts = 0; iElts < nElements; iElts++) {
    Int32 lineType = m_Elements[iElts]->m_Lines[0].GetType();
    if (lineTypeFilter == -1 || lineTypeFilter == lineType) {
      m_Elements[iElts]->addToSpectrumModel(spectralAxis, modelfluxAxis,
                                            getContinuumFluxAxis(), m_Redshift);
    }
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
    for (Int32 t = r.GetBegin(); t < r.GetEnd(); t++) {
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
