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
#include "RedshiftLibrary/linemodel/hybridfitter.h"
#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;
using namespace std;

CHybridFitter::CHybridFitter(
    const std::shared_ptr<CLMEltListVector> &elementsVector,
    const CCSpectrumVectorPtr &inputSpcs,
    const CTLambdaRangePtrVector &lambdaRanges,
    const CSpcModelVectorPtr &spectrumModels, const CLineMap &restLineList,
    const CSpectraGlobalIndex &spcIndex, bool enableAmplitudeOffsets,
    bool enableLambdaOffsetsFit)
    : CSvdFitter(elementsVector, inputSpcs, lambdaRanges, spectrumModels,
                 restLineList, spcIndex, enableAmplitudeOffsets,
                 enableLambdaOffsetsFit)

{
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();

  if (ps->GetScoped<std::string>("lineRatioType") == "rules")
    m_opt_enable_improveBalmerFit = ps->GetScoped<bool>("improveBalmerFit");
}

void CHybridFitter::doFit(Float64 redshift) {

  m_spectraIndex.setAtBegining(); // temporary multiobs implementation
  // fit the amplitudes of each element independently, unless there is overlap
  fitAmplitudesHybrid(redshift);

  // apply a continuum iterative re-estimation with lines removed from the
  // initial spectrum
  Int32 nIt = m_cont_reestim_iterations;
  Int32 it = 0;
  while (it < nIt) {
    //    m_Regulament.Apply(m_Elements);

    //*
    // iterative continuum estimation :: RAW SLOW METHOD
    getModel().refreshModel();
    Float64 enhanceLines = 0;
    //*
    if (nIt > 2 * it && nIt > 3.0 && it <= 3) {
      enhanceLines = 2.0 - ((Float64)it * 0.33);
    }

    //*/
    /*
      if(it==0 && nIt>1){
      enhanceLines = 1.5;
      }
    */
    getModel().EstimateSpectrumContinuum(enhanceLines);
    getModel().initModelWithContinuum();

    fitAmplitudesHybrid(redshift);
    it++;
  }
}

/**
 * \brief Tries to fit subelements considering their overlap.
 * For each entry in GetElementsIndicesInsideLambdaRange:
 *   If subelement in the entry already fitted, go for the next entry.
 *   getOverlappingElements for the fitted subelements.
 *   If the overlap is smaller than 2, call fitAmplitude on the entry.
 *   If the overlap is greater than or equal to 2:
 *     Call fitAmplitudeLinSolve with the subelements as argument.
 *     Store all non-negative fits.
 *     Set to 0.0 all negative fits.
 *     If the size of non-negative fits is 1, call the entry's fitAmplitude.
 *     If the size of non-negative fits is not 1:
 *       If the size of non-negative fits is greater than 1:
 *         Call fitAmplitudesLinSolve with the indexes of the non-negative
 *subelements. If the above call return is different than 1: For each
 *non-negative subelement, if the amplitude fitted is greater than 0, call
 *fitAmplitude on its entry. Else, SetElementAmplitude to 0. Update the index of
 *already-fitted subelements.
 **/
void CHybridFitter::fitAmplitudesHybrid(Float64 redshift) {

  m_spectraIndex.setAtBegining(); // dummy implementation

  TInt32List validEltsIdx = m_ElementsVector->getValidElementIndices();
  TInt32Set indexesFitted;
  for (Int32 iElts : validEltsIdx) {

    // skip if already fitted
    if (std::find(indexesFitted.cbegin(), indexesFitted.cend(), iElts) !=
        indexesFitted.cend())
      continue;

    TInt32List overlappingInds = getElementList().getOverlappingElements(
        iElts, indexesFitted, redshift, OVERLAP_THRES_HYBRID_FIT);

    // setting the fitting group info
    for (Int32 overlapping_iElt : overlappingInds) {
      std::string fitGroupTag = boost::str(boost::format("hy%d") % iElts);
      m_ElementsVector->getElementsParams()[overlapping_iElt]
          ->SetFittingGroupInfo(fitGroupTag);
    }

    Log.LogDebug(Formatter() << "    model: hybrid fit: #" << iElts
                             << " - N overlapping=" << overlappingInds.size());
    for (Int32 ifit = 0; ifit < ssize(overlappingInds); ifit++) {
      Log.LogDebug(Formatter()
                   << "    model: hybrid fit:     overlapping #" << ifit
                   << " - eltIdx=" << overlappingInds[ifit]);
    }
    if (isIndividualFitEnabled() && overlappingInds.size() < 2) {
      m_spectraIndex.setAtBegining(); // temporary multiobs implementation
      Log.LogDebug("    model: hybrid fit:     Individual fit");
      fitAmplitudeAndLambdaOffset(iElts, redshift, undefIdx,
                                  m_enableLambdaOffsetsFit);
      m_spectraIndex.setAtBegining(); // temporary multiobs implementation

    } else {
      m_spectraIndex.setAtBegining(); // temporary multiobs implementation

      Log.LogDebug("    model: hybrid fit:     Joint fit");
      fitAmplitudesLinSolveAndLambdaOffset(overlappingInds,
                                           m_enableLambdaOffsetsFit, redshift);
      m_spectraIndex.setAtBegining(); // temporary multiobs implementation
    }

    // update the already fitted list
    for (Int32 overlapping_iElt : overlappingInds) {
      indexesFitted.insert(overlapping_iElt);
    }
  }

  if (m_opt_enable_improveBalmerFit) {
    improveBalmerFit(redshift);
  }
}

TStringList CHybridFitter::initEmissionBalmerTags() {
  return {linetags::halpha_em, linetags::hbeta_em, linetags::hgamma_em,
          linetags::hdelta_em};
}

TStringList CHybridFitter::initAbsorptionBalmerTags() {
  return {linetags::halpha_abs, linetags::hbeta_abs, linetags::hgamma_abs,
          linetags::hdelta_abs};
}

std::vector<TStringList> CHybridFitter::initAdditionalTags() {
  // Additional lines to be fitted with the Balmer lines, WARNING: only
  // EMISSION for now !!
  TStringList linetagsNII = {linetags::niia_em, linetags::niib_em};
  TStringList empty;
  return {linetagsNII, empty, empty, empty};
}

// return error: 1=can't find element index, 2=Abs_width not high enough
// compared to Em_width
void CHybridFitter::improveBalmerFit(Float64 redshift) {
  auto linetagsE = initEmissionBalmerTags();
  auto linetagsA = initAbsorptionBalmerTags();
  auto linetagsMore = initAdditionalTags();

  if (!validateTagSizes(linetagsE, linetagsA, linetagsMore)) {
    return;
  }

  for (Int32 itag = 0; itag < ssize(linetagsE); itag++) {
    processBalmerPair(itag, redshift, linetagsE, linetagsA, linetagsMore);
  }
}

void CHybridFitter::processBalmerPair(
    Int32 itag, Float64 redshift, const TStringList &linetagsE,
    const TStringList &linetagsA,
    const std::vector<TStringList> &linetagsMore) {
  auto const &[iEltE, lineE_id] = m_ElementsVector->findElementIndex(
      linetagsE[itag], CLine::EType::nType_Emission);
  auto const &[iEltA, lineA_id] = m_ElementsVector->findElementIndex(
      linetagsA[itag], CLine::EType::nType_Absorption);

  if (!isValidBalmerPair(iEltE, iEltA, lineE_id, lineA_id)) {
    return;
  }

  auto [ilinesMore, linesMoreIds] = collectAdditionalLines(itag, linetagsMore);

  if (!widthConditionSatisfied(iEltE, lineE_id, iEltA, lineA_id, redshift)) {
    return;
  }

  attemptBalmerRefit(iEltA, lineA_id, iEltE, lineE_id, ilinesMore, linesMoreIds,
                     redshift);
}

bool CHybridFitter::isValidBalmerPair(Int32 iEltE, Int32 iEltA, Int32 lineE_id,
                                      Int32 lineA_id) {
  if (iEltE == undefIdx || iEltA == undefIdx)
    return false;
  if (getElementList()[iEltE]->GetSize() > 1 ||
      getElementList()[iEltA]->GetSize() > 1)
    return false;
  if (getElementsParams()[iEltE]->isNotFittable() ||
      getElementsParams()[iEltA]->isNotFittable())
    return false;
  return true;
}

std::pair<TInt32List, TInt32List> CHybridFitter::collectAdditionalLines(
    Int32 itag, const std::vector<TStringList> &linetagsMore) {
  TInt32List ilinesMore, ids;
  for (Int32 imore = 0; imore < ssize(linetagsMore[itag]); imore++) {
    auto const &[iElt, id] = m_ElementsVector->findElementIndex(
        linetagsMore[itag][imore], CLine::EType::nType_Emission);
    if (iElt == undefIdx || getElementsParams()[iElt]->isNotFittable())
      continue;

    ilinesMore.push_back(iElt);
    ids.push_back(id);
  }
  std::sort(ilinesMore.begin(), ilinesMore.end());
  ilinesMore.erase(std::unique(ilinesMore.begin(), ilinesMore.end()),
                   ilinesMore.end());
  return {ilinesMore, ids};
}

bool CHybridFitter::widthConditionSatisfied(Int32 iEltE, Int32 lineE_id,
                                            Int32 iEltA, Int32 lineA_id,
                                            Float64 redshift) {
  const Float64 Threshold = 2.0;
  auto const &[muE, sigmaE] =
      getElementList()[iEltE]->getObservedPositionAndLineWidth(redshift,
                                                               lineE_id, false);
  auto const &[muA, sigmaA] =
      getElementList()[iEltA]->getObservedPositionAndLineWidth(redshift,
                                                               lineA_id, false);
  return sigmaA >= Threshold * sigmaE;
}

void CHybridFitter::attemptBalmerRefit(Int32 iEltA, Int32 lineA_id, Int32 iEltE,
                                       Int32 lineE_id,
                                       const TInt32List &ilinesMore,
                                       const TInt32List &idsMore,
                                       Float64 redshift) {
  Float64 modelErr_init = getModelResidualRmsUnderElements({iEltA}, true);

  // collect amps before refit
  auto [ampA, errA] = getAmplitudeAndError(iEltA, lineA_id);
  auto [ampE, errE] = getAmplitudeAndError(iEltE, lineE_id);
  auto [ampsMore, errsMore] = collectAmplitudes(ilinesMore);

  // refit
  TInt32List eltsIdx = {iEltA, iEltE};
  eltsIdx.insert(eltsIdx.end(), ilinesMore.begin(), ilinesMore.end());

  TFloat64List ampsfitted, errorsfitted;
  fitAmplitudesLinSolve(eltsIdx, ampsfitted, errorsfitted, redshift);

  // check improvement
  getModel().refreshModelUnderElements(eltsIdx);
  Float64 modelErr_withfit = getModelResidualRmsUnderElements({iEltA}, true);
  if (modelErr_withfit > modelErr_init) {
    restoreAmplitudes(iEltA, lineA_id, ampA, errA, iEltE, lineE_id, ampE, errE,
                      ilinesMore, idsMore, ampsMore, errsMore);
  }
}

std::pair<Float64, Float64> CHybridFitter::getAmplitudeAndError(Int32 iElt,
                                                                Int32 lineId) {
  Float64 amp = getElementsParams()[iElt]->GetFittedAmplitude(lineId);
  Float64 ampErr = getElementsParams()[iElt]->GetFittedAmplitudeStd(lineId);
  return {amp, ampErr};
}

std::pair<TFloat64List, TFloat64List>
CHybridFitter::collectAmplitudes(const TInt32List &ilines) {
  TFloat64List amps;
  TFloat64List ampErrors;

  for (Int32 i = 0; i < ssize(ilines); ++i) {
    Int32 idx = ilines[i];
    Float64 amp = getElementsParams()[idx]->GetFittedAmplitude(0);
    Float64 ampErr = getElementsParams()[idx]->GetFittedAmplitudeStd(0);
    amps.push_back(amp);
    ampErrors.push_back(ampErr);
  }

  return {amps, ampErrors};
}

void CHybridFitter::restoreAmplitudes(Int32 iEltA, Int32 lineA_id, Float64 ampA,
                                      Float64 errA, Int32 iEltE, Int32 lineE_id,
                                      Float64 ampE, Float64 errE,
                                      const TInt32List &ilinesMore,
                                      const TInt32List &idsMore,
                                      const TFloat64List &ampsMore,
                                      const TFloat64List &errsMore) {
  // Restore absorption and emission lines
  Float64 nominal_ampA =
      getElementsParams()[iEltA]->GetNominalAmplitude(lineA_id);
  Float64 nominal_ampE =
      getElementsParams()[iEltE]->GetNominalAmplitude(lineE_id);

  m_ElementsVector->SetElementAmplitude(iEltA, ampA / nominal_ampA,
                                        errA / nominal_ampA);
  m_ElementsVector->SetElementAmplitude(iEltE, ampE / nominal_ampE,
                                        errE / nominal_ampE);

  // Restore additional lines
  for (Int32 i = 0; i < ssize(ilinesMore); ++i) {
    Int32 idx = ilinesMore[i];
    Float64 nominal_amp =
        getElementsParams()[idx]->GetNominalAmplitude(idsMore[i]);
    m_ElementsVector->SetElementAmplitude(idx, ampsMore[i] / nominal_amp,
                                          errsMore[i] / nominal_amp);
  }
}

bool CHybridFitter::validateTagSizes(const TStringList &E, const TStringList &A,
                                     const std::vector<TStringList> &More) {
  return E.size() == A.size() && E.size() == More.size();
}