
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

#include "RedshiftLibrary/linemodel/abstractfitter.h"
#include "RedshiftLibrary/linemodel/hybridfitter.h"
#ifdef LBFGSBFITTER
#include "RedshiftLibrary/linemodel/gaussianfit/lbfgsbfitter.h"
#endif
#include "RedshiftLibrary/linemodel/individualfitter.h"
#include "RedshiftLibrary/linemodel/onesfitter.h"
#include "RedshiftLibrary/linemodel/randomfitter.h"
#include "RedshiftLibrary/linemodel/svdfitter.h"
#include "RedshiftLibrary/linemodel/svdlcfitter.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;

CAbstractFitter::CAbstractFitter(
    const std::shared_ptr<CLMEltListVector> &elementsVector,
    const CCSpectrumVectorPtr &inputSpcs,
    const CTLambdaRangePtrVector &lambdaRanges,
    const CSpcModelVectorPtr &spectrumModels, const CLineMap &restLineList,
    const shared_ptr<Int32> &curObsPtr, bool enableAmplitudeOffsets,
    bool enableLambdaOffsetsFit)
    : m_ElementsVector(elementsVector), m_inputSpcs(inputSpcs),
      m_RestLineList(restLineList), m_lambdaRanges(lambdaRanges),
      m_models(spectrumModels), m_curObs(curObsPtr),
      m_enableAmplitudeOffsets(enableAmplitudeOffsets),
      m_enableLambdaOffsetsFit(enableLambdaOffsetsFit) {
  m_nbElements = m_ElementsVector->getNbElements();

  std::string method = Context.GetCurrentMethod();
  CAutoScope autoscope(Context.m_ScopeStack, "lineModel");
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();
  bool useAsymProfile = ps->GetScoped<std::string>("lya.profile") == "asym";
  if (useAsymProfile &
      (method == "lineModelSolve" || method == "lineMeasSolve")) {
    CAutoScope autoscope(Context.m_ScopeStack, "lya");
    m_opt_lya_fit_asym_min = ps->GetScoped<Float64>("asymProfile.asymFitMin");
    m_opt_lya_fit_asym_max = ps->GetScoped<Float64>("asymProfile.asymFitMax");
    m_opt_lya_fit_asym_step = ps->GetScoped<Float64>("asymProfile.asymFitStep");
    m_opt_lya_fit_width_min = ps->GetScoped<Float64>("asymProfile.widthFitMin");
    m_opt_lya_fit_width_max = ps->GetScoped<Float64>("asymProfile.widthFitMax");
    m_opt_lya_fit_width_step =
        ps->GetScoped<Float64>("asymProfile.widthFitStep");
    m_opt_lya_fit_delta_min = ps->GetScoped<Float64>("asymProfile.deltaFitMin");
    m_opt_lya_fit_delta_max = ps->GetScoped<Float64>("asymProfile.deltaFitMax");
    m_opt_lya_fit_delta_step =
        ps->GetScoped<Float64>("asymProfile.deltaStepMax");
  }
}

std::shared_ptr<CAbstractFitter> CAbstractFitter::makeFitter(
    std::string fittingMethod,
    const std::shared_ptr<CLMEltListVector> &elementsVector,
    const CCSpectrumVectorPtr &inputSpcs,
    const CTLambdaRangePtrVector &lambdaRanges,
    const CSpcModelVectorPtr &spectrumModels, const CLineMap &restLineList,
    std::shared_ptr<CContinuumManager> continuumManager,
    const std::shared_ptr<Int32> &curObsPtr, bool enableAmplitudeOffsets,
    bool enableLambdaOffsetsFit) {
  if (fittingMethod == "hybrid")
    return std::make_shared<CHybridFitter>(
        elementsVector, inputSpcs, lambdaRanges, spectrumModels, restLineList,
        curObsPtr, enableAmplitudeOffsets, enableLambdaOffsetsFit);
#ifdef LBFGSBFITTER
  else if (fittingMethod == "lbfgsb")
    return std::make_shared<CLbfgsbFitter>(
        elementsVector, inputSpcs, lambdaRanges, spectrumModels, restLineList,
        curObsPtr, enableAmplitudeOffsets, enableLambdaOffsetsFit);
#endif
  else if (fittingMethod == "svd")
    return std::make_shared<CSvdFitter>(
        elementsVector, inputSpcs, lambdaRanges, spectrumModels, restLineList,
        curObsPtr, enableAmplitudeOffsets, enableLambdaOffsetsFit);
  else if (fittingMethod == "svdlc")
    return std::make_shared<CSvdlcFitter>(
        elementsVector, inputSpcs, lambdaRanges, spectrumModels, restLineList,
        curObsPtr, continuumManager);
  else if (fittingMethod == "svdlcp2")
    return std::make_shared<CSvdlcFitter>(
        elementsVector, inputSpcs, lambdaRanges, spectrumModels, restLineList,
        curObsPtr, continuumManager, 2);

  else if (fittingMethod == "ones")
    return std::make_shared<COnesFitter>(elementsVector, inputSpcs,
                                         lambdaRanges, spectrumModels,
                                         restLineList, curObsPtr);
  else if (fittingMethod == "random")
    return std::make_shared<CRandomFitter>(elementsVector, inputSpcs,
                                           lambdaRanges, spectrumModels,
                                           restLineList, curObsPtr);
  else if (fittingMethod == "individual")
    return std::make_shared<CIndividualFitter>(elementsVector, inputSpcs,
                                               lambdaRanges, spectrumModels,
                                               restLineList, curObsPtr);
  else
    THROWG(INTERNAL_ERROR, Formatter()
                               << "Unknown fitting method " << fittingMethod);
}

void CAbstractFitter::logParameters() {
  Log.LogDetail(Formatter()
                << " m_opt_lya_fit_asym_min" << m_opt_lya_fit_asym_min);
  Log.LogDetail(Formatter()
                << " m_opt_lya_fit_asym_max" << m_opt_lya_fit_asym_max);
  Log.LogDetail(Formatter()
                << " m_opt_lya_fit_asym_step" << m_opt_lya_fit_asym_step);
  Log.LogDetail(Formatter()
                << " m_opt_lya_fit_width_min" << m_opt_lya_fit_width_min);
  Log.LogDetail(Formatter()
                << " m_opt_lya_fit_width_max" << m_opt_lya_fit_width_max);
  Log.LogDetail(Formatter()
                << " m_opt_lya_fit_width_step" << m_opt_lya_fit_width_step);
  Log.LogDetail(Formatter()
                << " m_opt_lya_fit_delta_min" << m_opt_lya_fit_delta_min);
  Log.LogDetail(Formatter()
                << " m_opt_lya_fit_delta_max" << m_opt_lya_fit_delta_max);
  Log.LogDetail(Formatter()
                << " m_opt_lya_fit_delta_step" << m_opt_lya_fit_delta_step);
}

void CAbstractFitter::fit(Float64 redshift) {
  initFit(redshift);
  doFit(redshift);
};

void CAbstractFitter::initFit(Float64 redshift) {

  resetSupport(redshift);

  m_ElementsVector->setGlobalOutsideLambdaRangeFromSpectra();
  // prepare the Lya width and asym coefficients if the asymfit profile
  // option is met

  fitLyaProfile(redshift);

  m_ElementsVector->resetElementsFittingParam(m_enableAmplitudeOffsets);
}

void CAbstractFitter::resetSupport(Float64 redshift) {

  m_ElementsVector->resetLambdaOffsets();
  m_ElementsVector->resetAsymfitParams();

  for (*m_curObs = 0; *m_curObs < m_inputSpcs->size(); (*m_curObs)++) {
    // prepare the elements support
    const CSpectrumSpectralAxis &spectralAxis = getSpectrum().GetSpectralAxis();
    for (auto const &elt_ptr : getElementList()) {

      elt_ptr->prepareSupport(spectralAxis, redshift, getLambdaRange(),
                              m_enlarge_line_supports);
    }
  }
}

void CAbstractFitter::resetElementsFittingParam() {
  for (auto param_ptr : m_ElementsVector->getElementParam())
    param_ptr->resetFittingParams();
}

void CAbstractFitter::fitLyaProfile(Float64 redshift) {
  TInt32List idxEltIGM;
  std::vector<TInt32List> idxLineIGM;
  Int32 line_idx_LyaE;
  TInt32List line_indices_LyaE_copy;

  auto const indices_Igm = m_ElementsVector->getIgmLinesIndices();

  if (indices_Igm.empty())
    return;

  // ASym Profile
  {
    auto const &[elt_idx_LyaE, line_indices_LyaE] = indices_Igm.front();
    line_idx_LyaE = line_indices_LyaE.front();

    auto const &param_LyaE = m_ElementsVector->getElementParam()[elt_idx_LyaE];
    auto const &profile = param_LyaE->getLineProfile(line_idx_LyaE);

    if (profile->isAsymFit() && !m_ElementsVector->isOutsideLambdaRangeLine(
                                    elt_idx_LyaE, line_idx_LyaE)) {
      // find the best width and asym coeff. parameters
      TAsymParams const bestfitParams =
          fitAsymParameters(redshift, elt_idx_LyaE, line_idx_LyaE);

      param_LyaE->SetAsymfitParams(bestfitParams);
    }
  }

  // deal with symIgm profiles
  {
    std::vector<std::pair<Int32, TInt32List>> line_indices_tofit;
    for (auto const &[elt_idx_igmLine, line_indices_LyaE] : indices_Igm) {
      auto const &param_EltIgm =
          m_ElementsVector->getElementParam()[elt_idx_igmLine];

      if (isOutsideLambdaRange(elt_idx_igmLine))
        continue;
      auto line_indices_filtered = line_indices_LyaE;
      auto end = std::remove_if(
          line_indices_filtered.begin(), line_indices_filtered.end(),
          [&](Int32 idx) {
            return !param_EltIgm->getLineProfile(idx)->isSymIgmFit() ||
                   m_ElementsVector->isOutsideLambdaRangeLine(elt_idx_igmLine,
                                                              idx);
          });
      line_indices_filtered.erase(end, line_indices_filtered.end());
      if (!line_indices_filtered.empty())
        line_indices_tofit.push_back({elt_idx_igmLine, line_indices_filtered});
    }

    auto bestigmidx = line_indices_tofit.empty()
                          ? undefIdx
                          : fitAsymIGMCorrection(redshift, line_indices_tofit);
    for (auto const &[elt_idx, _] : line_indices_tofit) {
      auto const &param_EltIgm = m_ElementsVector->getElementParam()[elt_idx];
      param_EltIgm->SetSymIgmParams(TSymIgmParams(bestigmidx, redshift));
    }
  }
}

void CAbstractFitter::fitAmplitude(Int32 eltIndex, Float64 redshift,
                                   Int32 lineIdx) {

  if (isOutsideLambdaRange(eltIndex)) {

    getElementParam()[eltIndex]->m_sumCross = NAN;
    getElementParam()[eltIndex]->m_sumGauss = NAN;
    getElementParam()[eltIndex]->m_dtmFree = NAN;
    return;
  }

  Int32 nLines = getElementList()[eltIndex]->GetSize();
  getElementParam()[eltIndex]->m_sumCross = 0.;
  getElementParam()[eltIndex]->m_sumGauss = 0.;
  getElementParam()[eltIndex]->m_dtmFree = 0.;

  getElementParam()[eltIndex]->m_FittedAmplitudes.assign(nLines, NAN);
  getElementParam()[eltIndex]->m_FittedAmplitudesStd.assign(nLines, NAN);

  Int32 num = 0;
  for (*m_curObs = 0; *m_curObs < m_inputSpcs->size(); (*m_curObs)++) {
    const CSpectrumSpectralAxis &spectralAxis = getSpectrum().GetSpectralAxis();
    const CSpectrumFluxAxis &noContinuumfluxAxis =
        getModel().getSpcFluxAxisNoContinuum();
    const CSpectrumFluxAxis &continuumfluxAxis =
        getModel().getContinuumFluxAxis();

    num += getElementList()[eltIndex]->computeCrossProducts(
        redshift, spectralAxis, noContinuumfluxAxis, continuumfluxAxis,
        lineIdx);
  }

  if (num == 0 || getElementParam()[eltIndex]->m_sumGauss == 0) {
    Log.LogDebug(
        Formatter()
        << "CLineModelElement::fitAmplitude: Could not fit amplitude:    "
           " num="
        << num << ", mtm=" << getElementParam()[eltIndex]->m_sumGauss
        << " at line index " << lineIdx << " of elt " << eltIndex);
    getElementParam()[eltIndex]->m_sumGauss = NAN;
    getElementParam()[eltIndex]->m_dtmFree = NAN;
    getElementParam()[eltIndex]->m_sumCross = NAN;
    return;
  }

  bool allNaN = true;
  for (*m_curObs = 0; *m_curObs < m_inputSpcs->size(); (*m_curObs)++) {
    if (!std::isnan(getElementParam()[eltIndex]->getSumGauss())) {
      allNaN = false;
      break;
    }
  }

  if (allNaN) {
    getElementParam()[eltIndex]->m_sumCross = NAN;
    return;
  }
  getElementParam()[eltIndex]->m_sumCross =
      std::max(0.0, getElementParam()[eltIndex]->m_dtmFree);
  Float64 const A = getElementParam()[eltIndex]->m_sumCross /
                    getElementParam()[eltIndex]->m_sumGauss;
  Float64 const Astd = 1.0 / sqrt(getElementParam()[eltIndex]->m_sumGauss);
  m_ElementsVector->SetElementAmplitude(eltIndex, A, Astd);

  return;
}

void CAbstractFitter::setLambdaOffset(const TInt32List &EltsIdx,
                                      Int32 offsetCount) {

  Float64 offset = m_LambdaOffsetMin + m_LambdaOffsetStep * offsetCount;
  for (Int32 iE : EltsIdx)
    getElementParam()[iE]->SetAllOffsetsEnabled(offset);

  return;
}

bool CAbstractFitter::HasLambdaOffsetFitting(TInt32List EltsIdx,
                                             bool enableOffsetFitting) const {
  bool atLeastOneOffsetToFit = false;
  if (enableOffsetFitting) {
    for (Int32 iE : EltsIdx)
      for (const auto &line : getElementParam()[iE]->GetLines())
        // check if the line is to be fitted
        if (line.IsOffsetFitEnabled()) {
          atLeastOneOffsetToFit = true;
          break;
        }
  }
  return atLeastOneOffsetToFit;
}

Int32 CAbstractFitter::GetLambdaOffsetSteps(bool atLeastOneOffsetToFit) const {
  Int32 nSteps =
      atLeastOneOffsetToFit
          ? int((m_LambdaOffsetMax - m_LambdaOffsetMin) / m_LambdaOffsetStep +
                0.5)
          : 1;

  return nSteps;
}

/**
 * @brief CAbstractFitter::fitAmplitudeAndLambdaOffset
 * Fit the amplitudes and a unique Offset for all lines by using the
 * fitAmplitude() in a loop to tabulate the offset
 * @param spectralAxis
 * @param noContinuumfluxAxis
 * @param continuumfluxAxis
 * @param redshift
 * @param lineIdx
 */
void CAbstractFitter::fitAmplitudeAndLambdaOffset(Int32 eltIndex,
                                                  Float64 redshift,
                                                  Int32 lineIdx,
                                                  bool enableOffsetFitting) {

  bool atLeastOneOffsetToFit =
      HasLambdaOffsetFitting({eltIndex}, enableOffsetFitting);
  Int32 nSteps = GetLambdaOffsetSteps(atLeastOneOffsetToFit);

  if (!atLeastOneOffsetToFit) {
    Log.LogDebug("    multiline: no offsets to fit");
    nSteps = 1;
  } else {
    Log.LogDebug(Formatter() << "    multiline: offsets to fit n=" << nSteps);
  }

  Float64 bestMerit = DBL_MAX;
  Int32 idxBestMerit = -1;
  for (Int32 iO = 0; iO < nSteps; iO++) {
    // set offset value
    if (atLeastOneOffsetToFit)
      setLambdaOffset({eltIndex}, iO);

    // fit for this offset
    fitAmplitude(eltIndex, redshift, lineIdx);

    // check fitting
    if (atLeastOneOffsetToFit) {
      Float64 fit = getLeastSquareMeritFast(eltIndex);
      if (fit < bestMerit) {
        bestMerit = fit;
        idxBestMerit = iO;
      }
    }
  }

  if (idxBestMerit == -1 || !atLeastOneOffsetToFit)
    return;

  // set offset value
  if (atLeastOneOffsetToFit)
    setLambdaOffset({eltIndex}, idxBestMerit);
  // fit again for this offset
  fitAmplitude(eltIndex, redshift, lineIdx);
}

/**
 * \brief Get the squared difference by fast method proposed by D. Vibert
 **/
Float64 CAbstractFitter::getLeastSquareMeritFast(Int32 eltIdx) const {
  Float64 fit = 0.; // TODO restore getLeastSquareContinuumMeritFast();
  Int32 istart = 0;
  Int32 iend = getElementList().size();
  if (eltIdx != undefIdx) {
    istart = eltIdx;
    iend = eltIdx + 1;
  }
  for (Int32 iElts = istart; iElts < iend; iElts++) {
    Float64 dtm = getElementParam()[iElts]->getSumCross();
    Float64 mtm = getElementParam()[iElts]->getSumGauss();
    Float64 a = getElementList().GetElementAmplitude(
        iElts); //[iElts]->GetFitAmplitude();
    Float64 term1 = a * a * mtm;
    Float64 term2 = -2. * a * dtm;
    fit += term1 + term2;
  }

  Log.LogDebug(
      Formatter() << "CLineModelFitting::getLeastSquareMerit fit fast = "
                  << fit);
  return fit;
}

TAsymParams CAbstractFitter::fitAsymParameters(Float64 redshift, Int32 idxLyaE,
                                               const Int32 &idxLineLyaE) {

  // 3. find the best width and asym coeff. parameters
  Float64 widthCoeffStep = m_opt_lya_fit_width_step;
  Float64 widthCoeffMin = m_opt_lya_fit_width_min;
  Float64 widthCoeffMax = m_opt_lya_fit_width_max;
  Int32 nWidthSteps =
      int((widthCoeffMax - widthCoeffMin) / widthCoeffStep + 1.5);
  Float64 asymCoeffStep = m_opt_lya_fit_asym_step;
  Float64 asymCoeffMin = m_opt_lya_fit_asym_min;
  Float64 asymCoeffMax = m_opt_lya_fit_asym_max;
  Int32 nAsymSteps = int((asymCoeffMax - asymCoeffMin) / asymCoeffStep + 1.5);
  Float64 deltaStep = m_opt_lya_fit_delta_step;
  Float64 deltaMin = m_opt_lya_fit_delta_min;
  Float64 deltaMax = m_opt_lya_fit_delta_max;
  Int32 nDeltaSteps = int((deltaMax - deltaMin) / deltaStep + 1.5);

  TAsymParams bestparams = ASYMF_DEFAULT_PARAMS;
  Float64 meritMin = DBL_MAX;

  TInt32List filterEltsIdxLya(1, idxLyaE);

  auto const &param_LyaE = m_ElementsVector->getElementParam()[idxLyaE];

  for (Int32 iDelta = 0; iDelta < nDeltaSteps; iDelta++) {
    Float64 delta = deltaMin + deltaStep * iDelta;
    for (Int32 iWidth = 0; iWidth < nWidthSteps; iWidth++) {
      Float64 asymWidthCoeff = widthCoeffMin + widthCoeffStep * iWidth;
      for (Int32 iAsym = 0; iAsym < nAsymSteps; iAsym++) {
        Float64 asymAlphaCoeff = asymCoeffMin + asymCoeffStep * iAsym;
        param_LyaE->SetAsymfitParams({asymWidthCoeff, asymAlphaCoeff, delta});
        fitAmplitude(idxLyaE, redshift, idxLineLyaE);
        bool fitIsvalid = param_LyaE->GetElementAmplitude() > 0.0;

        Float64 m = NAN;
        if (fitIsvalid) {
          if (1) {
            for (*m_curObs = 0; *m_curObs < m_models->size(); (*m_curObs)++)
              getModel().refreshModelUnderElements(filterEltsIdxLya,
                                                   idxLineLyaE);
            m = getModelResidualRmsUnderElements({idxLyaE}, true);

          } else {
            m = getLeastSquareMeritFast(idxLyaE);
          }
          if (m < meritMin) {
            meritMin = m;
            bestparams = param_LyaE->GetAsymfitParams(0);
          }
        }

        Log.LogDebug(Formatter()
                     << "Fitting Lya Profile: width=" << asymWidthCoeff
                     << ", asym=" << asymAlphaCoeff << ", delta=" << delta);
        Log.LogDebug(Formatter() << "Fitting Lya Profile: merit=%" << m);
        Log.LogDebug(Formatter() << "Fitting Lya Profile: idxLyaE=" << idxLyaE
                                 << ", idxLineLyaE=" << idxLineLyaE);
      }
    }
  }
  Log.LogDebug(Formatter() << "Lya Profile found: width=" << bestparams.sigma
                           << ", asym=" << bestparams.alpha
                           << ", delta=" << bestparams.delta);
  return bestparams;
}

Int32 CAbstractFitter::fitAsymIGMCorrection(
    Float64 redshift,
    std::vector<std::pair<Int32, TInt32List>> const &idxLines) {

  Float64 meritMin = DBL_MAX;
  Int32 bestIgmIdx = undefIdx;

  Int32 igmCount = m_ElementsVector->getElementParam()[idxLines.front().first]
                       ->getLineProfile(idxLines.front().second.front())
                       ->getIGMIdxCount();
  for (Int32 igmIdx = 0; igmIdx < igmCount; igmIdx++) {
    bool fitIsValid = false;
    for (auto const &[elt_idx, _] : idxLines) {
      auto const &param_EltIgm = m_ElementsVector->getElementParam()[elt_idx];
      param_EltIgm->SetSymIgmParams(TSymIgmParams(igmIdx, redshift));
      fitAmplitude(elt_idx, redshift);
      fitIsValid |= param_EltIgm->GetElementAmplitude() > 0.0;
    }
    if (fitIsValid) {
      TInt32List elt_indices;
      for (auto const [elt_idx, _] : idxLines)
        elt_indices.push_back(elt_idx);
      for (*m_curObs = 0; *m_curObs < m_inputSpcs->size(); (*m_curObs)++)
        getModel().refreshModelUnderElements(elt_indices);
      Float64 m = getModelResidualRmsUnderElements(elt_indices, true);

      if (m < meritMin) {
        meritMin = m;
        bestIgmIdx = igmIdx;
      }
    }
  }

  return bestIgmIdx;
}

Float64
CAbstractFitter::getModelResidualRmsUnderElements(TInt32List const &EltsIdx,
                                                  bool with_continuum) {
  Float64 fit_allObs = 0;
  Float64 sumErr_allObs = 0;
  Int32 nb_nan = 0;
  for (*m_curObs = 0; *m_curObs < m_inputSpcs->size(); (*m_curObs)++) {
    auto [fit, sumErr] = getModel().getModelSquaredResidualUnderElements(
        EltsIdx, with_continuum);
    if (fit == 0.0)
      continue;
    if (std::isnan(fit)) {
      nb_nan++;
      continue;
    }
    fit_allObs += fit;
    sumErr_allObs += sumErr;
  }
  if (nb_nan == m_inputSpcs->size())
    return NAN;
  return sumErr_allObs != 0.0 ? sqrt(fit_allObs / sumErr_allObs) : NAN;
}

bool CAbstractFitter::isOutsideLambdaRange(Int32 elt_index) {

  bool allOutsideLambdaRange = true;
  for (*m_curObs = 0; *m_curObs < m_inputSpcs->size(); (*m_curObs)++) {
    if (!getElementList()[elt_index]->IsOutsideLambdaRange()) {
      allOutsideLambdaRange = false;
      break;
    }
  }
  return allOutsideLambdaRange;
}
