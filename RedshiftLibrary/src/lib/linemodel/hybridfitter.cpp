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
    CLineModelElementList &elements,
    std::shared_ptr<const CSpectrum> inputSpectrum,
    std::shared_ptr<const TLambdaRange> lambdaRange,
    std::shared_ptr<CSpectrumModel> spectrumModel,
    const TLineVector &restLineList,
    const std::vector<TLineModelElementParam_ptr> &elementParam)
    : CSvdFitter(elements, inputSpectrum, lambdaRange, spectrumModel,
                 restLineList, elementParam)

{
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();

  if (ps->GetScoped<std::string>("linemodel.lineRatioType") == "rules")
    m_opt_enable_improveBalmerFit =
        ps->GetScoped<bool>("linemodel.improveBalmerFit");
}

void CHybridFitter::fit(Float64 redshift) {

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
    m_model->refreshModel();
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
    m_model->EstimateSpectrumContinuum(enhanceLines);
    m_model->initModelWithContinuum();

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
Int32 CHybridFitter::fitAmplitudesHybrid(Float64 redshift) {

  if (m_enableAmplitudeOffsets)
    m_Elements.resetAmplitudeOffset();

  TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();
  TInt32List indexesFitted;
  for (Int32 iValidElts = 0; iValidElts < validEltsIdx.size(); iValidElts++) {
    Int32 iElts = validEltsIdx[iValidElts];
    // skip if already fitted
    bool alreadyfitted = false;
    for (Int32 i = 0; i < indexesFitted.size(); i++) {
      if (iElts == indexesFitted[i]) {
        alreadyfitted = true;
        break;
      }
    }
    if (alreadyfitted) {
      continue;
    }

    TInt32List overlappingInds = m_Elements.getOverlappingElements(
        iElts, indexesFitted, redshift, OVERLAP_THRES_HYBRID_FIT);

    // setting the fitting group info
    for (Int32 ifit = 0; ifit < overlappingInds.size(); ifit++) {
      std::string fitGroupTag = boost::str(boost::format("hy%d") % iValidElts);
      m_Elements[overlappingInds[ifit]]->SetFittingGroupInfo(fitGroupTag);
    }

    // Log.LogDebug( "Redshift: %f", redshift);
    Log.LogDebug("    model: hybrid fit: #%d - N overlapping=%d", iValidElts,
                 overlappingInds.size());
    for (Int32 ifit = 0; ifit < overlappingInds.size(); ifit++) {
      Log.LogDebug("    model: hybrid fit:     overlapping #%d - eltIdx=%d",
                   ifit, overlappingInds[ifit]);
    }
    if (!m_enableAmplitudeOffsets && overlappingInds.size() < 2) {
      Log.LogDebug("    model: hybrid fit:     Individual fit");
      fitAmplitudeAndLambdaOffset(iElts, redshift, undefIdx,
                                  m_enableLambdaOffsetsFit);
    } else {
      Log.LogDebug("    model: hybrid fit:     Joint fit");
      TFloat64List ampsfitted;
      TFloat64List errorsfitted;
      Int32 retVal = fitAmplitudesLinSolveAndLambdaOffset(
          overlappingInds, ampsfitted, errorsfitted, m_enableLambdaOffsetsFit,
          redshift);
      // if all the amplitudes fitted don't have the same sign, do it
      // separately
      TInt32List overlappingIndsSameSign;
      if (retVal != 1 && ampsfitted.size() > 0) {
        for (Int32 ifit = 0; ifit < overlappingInds.size(); ifit++) {
          if (ampsfitted[ifit] > 0) {
            overlappingIndsSameSign.push_back(overlappingInds[ifit]);
            // m_Elements[overlappingInds[ifit]]->fitAmplitude(spectralAxis,
            // spcFluxAxisNoContinuum, redshift);
          } else {
            m_Elements.SetElementAmplitude(overlappingInds[ifit], 0.0,
                                           errorsfitted[ifit]);
          }
        }
        // fit the rest of the overlapping elements (same sign) together
        if (!m_enableAmplitudeOffsets && overlappingIndsSameSign.size() == 1) {
          fitAmplitudeAndLambdaOffset(overlappingIndsSameSign[0], redshift,
                                      undefIdx, m_enableLambdaOffsetsFit);
        } else if (overlappingIndsSameSign.size() > 0) {
          Int32 retVal2 = fitAmplitudesLinSolveAndLambdaOffset(
              overlappingIndsSameSign, ampsfitted, errorsfitted,
              m_enableLambdaOffsetsFit, redshift);

          if (retVal2 != 1) {
            for (Int32 ifit = 0; ifit < overlappingIndsSameSign.size();
                 ifit++) {
              if (ampsfitted[ifit] > 0) {
                fitAmplitudeAndLambdaOffset(overlappingIndsSameSign[ifit],
                                            redshift, undefIdx,
                                            m_enableLambdaOffsetsFit);
              } else {
                m_Elements.SetElementAmplitude(overlappingIndsSameSign[ifit],
                                               0.0, errorsfitted[ifit]);
              }
            }
          }
        }
      }
    }

    // update the already fitted list
    for (Int32 i = 0; i < overlappingInds.size(); i++) {
      indexesFitted.push_back(overlappingInds[i]);
    }
  }

  if (m_opt_enable_improveBalmerFit) {
    improveBalmerFit(redshift);
  }

  return 0;
}

// return error: 1=can't find element index, 2=Abs_width not high enough
// compared to Em_width
Int32 CHybridFitter::improveBalmerFit(Float64 redshift) {

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
    return -1;
  }

  for (Int32 itag = 0; itag < linetagsE.size(); itag++) {
    std::string tagE = linetagsE[itag];
    std::string tagA = linetagsA[itag];

    Int32 ilineE = m_Elements.findElementIndex(tagE, CLine::nType_Emission);
    Int32 ilineA = m_Elements.findElementIndex(tagA, CLine::nType_Absorption);
    // Were the lines indexes found ?
    if (ilineE < 0 || ilineA < 0) {
      continue;
    }
    // for now only allow this process if Em and Abs line are single lines
    if (m_Elements[ilineE]->GetSize() > 1 ||
        m_Elements[ilineA]->GetSize() > 1) {
      continue;
    }
    Int32 subeIdxE = 0;
    Int32 subeIdxA = 0;

    // check if line is visible:
    if (m_Elements[ilineE]->IsOutsideLambdaRange())
      continue;

    // find the linesMore unique elements indexes
    TInt32List ilinesMore;
    for (Int32 imore = 0; imore < linetagsMore[itag].size(); imore++) {
      std::string tagMore = linetagsMore[itag][imore];
      Int32 ilineMore =
          m_Elements.findElementIndex(tagMore, CLine::nType_Emission);
      if (ilineMore < 0) {
        continue;
      }
      ilinesMore.push_back(ilineMore);
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
    m_Elements[ilineE]->getObservedPositionAndLineWidth(
        subeIdxE, redshift, muE, sigmaE,
        false); // do not apply Lya asym offset
    m_Elements[ilineA]->getObservedPositionAndLineWidth(
        subeIdxA, redshift, muA, sigmaA,
        false); // do not apply Lya asym offset
    if (sigmaA < AbsVSEmWidthCoeffThreshold * sigmaE) {
      continue;
    }

    // simulatneous fit with linsolve
    Float64 modelErr_init =
        m_model->getModelErrorUnderElement(ilineA, m_model->getSpcFluxAxis());
    Float64 ampA = m_Elements[ilineA]->GetFittedAmplitude(0);
    Float64 amp_errorA = m_Elements[ilineA]->GetFittedAmplitudeErrorSigma(0);
    Float64 ampE = m_Elements[ilineE]->GetFittedAmplitude(0);
    Float64 amp_errorE = m_Elements[ilineE]->GetFittedAmplitudeErrorSigma(0);
    TFloat64List ampsMore;
    TFloat64List ampErrorsMore;
    for (Int32 imore = 0; imore < ilinesMore.size(); imore++) {
      Float64 amp = m_Elements[ilinesMore[imore]]->GetFittedAmplitude(0);
      Float64 ampErr =
          m_Elements[ilinesMore[imore]]->GetFittedAmplitudeErrorSigma(0);
      ampsMore.push_back(amp);
      ampErrorsMore.push_back(ampErr);
    }

    TInt32List eltsIdx;
    eltsIdx.push_back(ilineA);
    eltsIdx.push_back(ilineE);
    for (Int32 imore = 0; imore < ilinesMore.size(); imore++) {
      eltsIdx.push_back(ilinesMore[imore]);
    }
    TFloat64List ampsfitted;
    TFloat64List errorsfitted;
    fitAmplitudesLinSolve(eltsIdx, ampsfitted, errorsfitted, redshift);

    // decide if the fit is better than previous amps
    TInt32List elts;
    elts.push_back(ilineA);
    elts.push_back(ilineE);
    for (Int32 imore = 0; imore < ilinesMore.size(); imore++) {
      elts.push_back(ilinesMore[imore]);
    }
    m_model->refreshModelUnderElements(elts);
    Float64 modelErr_withfit =
        m_model->getModelErrorUnderElement(ilineA, m_model->getSpcFluxAxis());
    if (modelErr_withfit > modelErr_init) {
      Float64 nominal_ampA = m_Elements[ilineA]->GetNominalAmplitude(0);
      Float64 nominal_ampE = m_Elements[ilineE]->GetNominalAmplitude(0);
      m_Elements[ilineA]->SetFittedAmplitude(ampA / nominal_ampA,
                                             amp_errorA / nominal_ampA);
      m_Elements[ilineE]->SetFittedAmplitude(ampE / nominal_ampE,
                                             amp_errorE / nominal_ampE);
      for (Int32 imore = 0; imore < ilinesMore.size(); imore++) {
        Float64 nominal_ampMore =
            m_Elements[ilinesMore[imore]]->GetNominalAmplitude(0);
        m_Elements[ilinesMore[imore]]->SetFittedAmplitude(
            ampsMore[imore] / nominal_ampMore,
            ampErrorsMore[imore] / nominal_ampMore);
      }
    }
  }

  return 0;
}
