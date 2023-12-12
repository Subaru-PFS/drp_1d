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
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;
using namespace std;

CHybridFitter::CHybridFitter(
    const std::shared_ptr<CLMEltListVector> &elementsVector,
    const CCSpectrumVectorPtr &inputSpcs,
    const CTLambdaRangePtrVector &lambdaRanges,
    const CSpcModelVectorPtr &spectrumModels, const CLineMap &restLineList,
    const std::shared_ptr<Int32> &curObsPtr, bool enableAmplitudeOffsets,
    bool enableLambdaOffsetsFit)
    : CSvdFitter(elementsVector, inputSpcs, lambdaRanges, spectrumModels,
                 restLineList, curObsPtr, enableAmplitudeOffsets,
                 enableLambdaOffsetsFit)

{
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();

  if (ps->GetScoped<std::string>("lineModel.lineRatioType") == "rules")
    m_opt_enable_improveBalmerFit =
        ps->GetScoped<bool>("lineModel.improveBalmerFit");
}

void CHybridFitter::doFit(Float64 redshift) {

  *m_curObs = 0; // temporary multiobs implementation
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
 * For each entry in GetModelValidElementsIndexes:
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

  *m_curObs = 0; // dummy implementation
  if (m_enableAmplitudeOffsets)
    getElementList().resetAmplitudeOffset();

  TInt32List validEltsIdx = getElementList().GetModelValidElementsIndexes();
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
      getElementList()[overlapping_iElt]->SetFittingGroupInfo(fitGroupTag);
    }

    // Log.LogDebug( "Redshift: %f", redshift);
    Log.LogDebug("    model: hybrid fit: #%d - N overlapping=%d", iElts,
                 overlappingInds.size());
    for (Int32 ifit = 0; ifit < overlappingInds.size(); ifit++) {
      Log.LogDebug("    model: hybrid fit:     overlapping #%d - eltIdx=%d",
                   ifit, overlappingInds[ifit]);
    }
    if (isIndividualFitEnabled() && overlappingInds.size() < 2) {
      Log.LogDebug("    model: hybrid fit:     Individual fit");
      fitAmplitudeAndLambdaOffset(iElts, redshift, undefIdx,
                                  m_enableLambdaOffsetsFit);
      *m_curObs = 0; // temporary multiobs implementation

    } else {
      Log.LogDebug("    model: hybrid fit:     Joint fit");
      fitAmplitudesLinSolveAndLambdaOffset(overlappingInds,
                                           m_enableLambdaOffsetsFit, redshift);
      *m_curObs = 0; // temporary multiobs implementation
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

// return error: 1=can't find element index, 2=Abs_width not high enough
// compared to Em_width
void CHybridFitter::improveBalmerFit(Float64 redshift) {

  // Emission Balmer lines
  TStringList linetagsE;
  linetagsE.push_back(linetags::halpha_em);
  linetagsE.push_back(linetags::hbeta_em);
  linetagsE.push_back(linetags::hgamma_em);
  linetagsE.push_back(linetags::hdelta_em);
  // Absorption Balmer lines
  TStringList linetagsA;
  linetagsA.push_back(linetags::halpha_abs);
  linetagsA.push_back(linetags::hbeta_abs);
  linetagsA.push_back(linetags::hgamma_abs);
  linetagsA.push_back(linetags::hdelta_abs);
  // Additional lines to be fitted with the Balmer lines, WARNING: only
  // EMISSION for now !!
  TStringList linetagsNII;
  linetagsNII.push_back(linetags::niia_em);
  linetagsNII.push_back(linetags::niib_em);
  TStringList linetagsVoid;
  std::vector<TStringList> linetagsMore;
  linetagsMore.push_back(linetagsNII);
  linetagsMore.push_back(linetagsVoid);
  linetagsMore.push_back(linetagsVoid);
  linetagsMore.push_back(linetagsVoid);

  if (linetagsE.size() != linetagsA.size() ||
      linetagsE.size() != linetagsMore.size()) {
    return;
  }

  for (Int32 itag = 0; itag < linetagsE.size(); itag++) {
    std::string tagE = linetagsE[itag];
    std::string tagA = linetagsA[itag];

    auto const &[iElt_lineE, lineE_id] =
        m_ElementsVector->findElementIndex(tagE, CLine::EType::nType_Emission);
    auto const &[iElt_lineA, lineA_id] = m_ElementsVector->findElementIndex(
        tagA, CLine::EType::nType_Absorption);
    // Were the lines indexes found ?
    if (iElt_lineE == undefIdx || iElt_lineA == undefIdx)
      continue;

    // for now only allow this process if Em and Abs line are single lines
    if (getElementList()[iElt_lineE]->GetSize() > 1 ||
        getElementList()[iElt_lineA]->GetSize() > 1) {
      continue;
    }

    // check if line is visible:
    if (getElementList()[iElt_lineE]->IsOutsideLambdaRange())
      continue;

    // find the linesMore unique elements indexes
    TInt32List ilinesMore;
    TInt32List linesMoreIds;
    for (Int32 imore = 0; imore < linetagsMore[itag].size(); imore++) {
      std::string tagMore = linetagsMore[itag][imore];
      auto const &[iElt_lineMore, lineMore_id] =
          m_ElementsVector->findElementIndex(tagMore,
                                             CLine::EType::nType_Emission);
      if (iElt_lineMore == undefIdx)
        continue;

      ilinesMore.push_back(iElt_lineMore);
      linesMoreIds.push_back(lineMore_id);
    }
    std::sort(ilinesMore.begin(), ilinesMore.end());
    ilinesMore.erase(std::unique(ilinesMore.begin(), ilinesMore.end()),
                     ilinesMore.end());
    for (Int32 imore = 0; imore < ilinesMore.size(); imore++) {
      Log.LogDebug("    model: balmerImprove more tags = %d",
                   ilinesMore[imore]);
    }

    // try if the width is significantly different: abs > em
    Float64 AbsVSEmWidthCoeffThreshold = 2.0;
    Float64 muE = NAN;
    Float64 muA = NAN;
    Float64 sigmaE = NAN;
    Float64 sigmaA = NAN;
    getElementList()[iElt_lineE]->getObservedPositionAndLineWidth(
        lineE_id, redshift, muE, sigmaE,
        false); // do not apply Lya asym offset
    getElementList()[iElt_lineA]->getObservedPositionAndLineWidth(
        lineA_id, redshift, muA, sigmaA,
        false); // do not apply Lya asym offset
    if (sigmaA < AbsVSEmWidthCoeffThreshold * sigmaE) {
      continue;
    }

    // simulatneous fit with linsolve
    Float64 modelErr_init = getModelErrorUnderElements({iElt_lineA}, true);
    Float64 ampA = getElementList()[iElt_lineA]->GetFittedAmplitude(lineA_id);
    Float64 amp_errorA =
        getElementList()[iElt_lineA]->GetFittedAmplitudeStd(lineA_id);
    Float64 ampE = getElementList()[iElt_lineE]->GetFittedAmplitude(lineE_id);
    Float64 amp_errorE =
        getElementList()[iElt_lineE]->GetFittedAmplitudeStd(lineE_id);
    TFloat64List ampsMore;
    TFloat64List ampErrorsMore;
    for (Int32 imore = 0; imore < ilinesMore.size(); imore++) {
      Float64 amp = getElementList()[ilinesMore[imore]]->GetFittedAmplitude(0);
      Float64 ampErr =
          getElementList()[ilinesMore[imore]]->GetFittedAmplitudeStd(0);
      ampsMore.push_back(amp);
      ampErrorsMore.push_back(ampErr);
    }

    TInt32List eltsIdx;
    eltsIdx.push_back(iElt_lineA);
    eltsIdx.push_back(iElt_lineE);
    for (Int32 imore = 0; imore < ilinesMore.size(); imore++) {
      eltsIdx.push_back(ilinesMore[imore]);
    }
    TFloat64List ampsfitted;
    TFloat64List errorsfitted;
    fitAmplitudesLinSolve(eltsIdx, ampsfitted, errorsfitted, redshift);

    // decide if the fit is better than previous amps
    getModel().refreshModelUnderElements(eltsIdx);
    Float64 modelErr_withfit = getModelErrorUnderElements({iElt_lineA}, true);
    if (modelErr_withfit > modelErr_init) {
      Float64 nominal_ampA =
          getElementList()[iElt_lineA]->GetNominalAmplitude(lineA_id);
      Float64 nominal_ampE =
          getElementList()[iElt_lineE]->GetNominalAmplitude(lineE_id);
      getElementList()[iElt_lineA]->SetElementAmplitude(
          ampA / nominal_ampA, amp_errorA / nominal_ampA);
      getElementList()[iElt_lineE]->SetElementAmplitude(
          ampE / nominal_ampE, amp_errorE / nominal_ampE);
      for (Int32 imore = 0; imore < ilinesMore.size(); imore++) {
        Float64 nominal_ampMore =
            getElementList()[ilinesMore[imore]]->GetNominalAmplitude(
                linesMoreIds[imore]);
        getElementList()[ilinesMore[imore]]->SetElementAmplitude(
            ampsMore[imore] / nominal_ampMore,
            ampErrorsMore[imore] / nominal_ampMore);
      }
    }
  }
}
