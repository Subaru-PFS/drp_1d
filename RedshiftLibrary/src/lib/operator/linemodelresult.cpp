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
#include "RedshiftLibrary/operator/linemodelresult.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/linemodel/templatesfitstore.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/statistics/deltaz.h"

#include <fstream>
#include <string>

using namespace NSEpic;

/**
 * @brief CLineModelResult::Init
 * Initializes the linemodel results, by:
 *  - setting the redshift grid
 *  - setting the size of chisquare vectors
 * @param redshifts
 * @param restLines
 * @param nTplshapes
 * @return
 */
void CLineModelResult::Init(TFloat64List redshifts,
                            CLineCatalog::TLineVector restLines,
                            Int32 nTemplates, Int32 nTplshapes,
                            TFloat64List tplshapesPriors) {
  if (tplshapesPriors.size() != nTplshapes) {
    THROWG(INTERNAL_ERROR,
           Formatter() << "LinemodelResult: tplshapeprior size ="
                       << tplshapesPriors.size() << "and tplshapes size="
                       << nTplshapes << "do not correspond");
  }

  Int32 nResults = redshifts.size();
  Status.resize(nResults);
  ChiSquare.resize(nResults);
  ScaleMargCorrection.resize(nResults);
  Redshifts = std::move(redshifts);
  restLineList = std::move(restLines);
  LineModelSolutions.resize(nResults);
  ContinuumModelSolutions.resize(nResults);

  ChiSquareContinuum.resize(nResults);
  ScaleMargCorrectionContinuum.resize(nResults);

  if (nTemplates > 0)
    ChiSquareTplContinuum.assign(nTemplates, TFloat64List(nResults, DBL_MAX));
  else
    ChiSquareTplContinuum.clear();

  // init the tplshape chisquare results
  if (nTplshapes > 0) {
    ChiSquareTplshapes.assign(nTplshapes, TFloat64List(nResults, DBL_MAX));
    ScaleMargCorrectionTplshapes.assign(nTplshapes,
                                        TFloat64List(nResults, 0.0));
    StrongELPresentTplshapes.assign(nTplshapes, TBoolList(nResults, false));
    StrongHalphaELPresentTplshapes.assign(nTplshapes,
                                          TBoolList(nResults, false));
    NLinesAboveSNRTplshapes.assign(nTplshapes, TInt32List(nResults, false));
    PriorTplshapes = std::move(tplshapesPriors);
    PriorLinesTplshapes.assign(nTplshapes, TFloat64List(nResults, 0.0));
  } else {
    ChiSquareTplshapes.clear();
    ScaleMargCorrectionTplshapes.clear();
    StrongELPresentTplshapes.clear();
    NLinesAboveSNRTplshapes.clear();
    PriorTplshapes.clear();
    PriorLinesTplshapes.clear();
  }
}

void CLineModelResult::SetChisquareTplContinuumResult(
    Int32 index_z, const shared_ptr<const CTemplatesFitStore> &tplFitStore) {
  const auto index_z_in_store =
      tplFitStore->GetRedshiftIndex(Redshifts[index_z]);
  if (index_z_in_store == -1)
    THROWG(INTERNAL_ERROR, "CLineModelResult::SetChisquareTplContinuumResult: "
                           "redshift not in fitstore");

  for (Int32 k = 0; k < tplFitStore->GetContinuumCount();
       k++) { // TODO: handle the use of more than one continuum in linemodel
    ChiSquareTplContinuum[k][index_z] =
        tplFitStore->GetFitValues(index_z_in_store, k).merit;
  }
}

void CLineModelResult::SetChisquareTplContinuumResultFromPrevious(
    Int32 index_z) {
  auto previous = index_z - 1;
  for (auto it = ChiSquareTplContinuum.begin(), e = ChiSquareTplContinuum.end();
       it != e; ++it)
    it->at(index_z) = it->at(previous);
}

void CLineModelResult::SetChisquareTplshapeResult(
    Int32 index_z, const TFloat64List &chisquareTplshape,
    const TFloat64List &scaleMargCorrTplshape,
    const TBoolList &strongEmissionLinePresentTplshape,
    const TBoolList &strongHalphaELPresentTplshapes,
    const TInt32List &nLinesAboveSNRTplshape,
    const TFloat64List &priorLinesTplshape) {
  if (chisquareTplshape.size() < 1)
    return;

  if (index_z >= Redshifts.size())
    THROWG(INTERNAL_ERROR,
           "CLineModelResult::SetChisquareTplshapeResult: invalid z index");

  if (chisquareTplshape.size() != ChiSquareTplshapes.size() ||
      chisquareTplshape.size() != scaleMargCorrTplshape.size() ||
      chisquareTplshape.size() != strongEmissionLinePresentTplshape.size() ||
      chisquareTplshape.size() != nLinesAboveSNRTplshape.size() ||
      chisquareTplshape.size() != priorLinesTplshape.size())
    THROWG(INTERNAL_ERROR, "CLineModelResult::SetChisquareTplshapeResult: "
                           "vector sizes do not match");

  for (Int32 k = 0; k < chisquareTplshape.size(); k++) {
    ChiSquareTplshapes[k][index_z] = chisquareTplshape[k];
    ScaleMargCorrectionTplshapes[k][index_z] = scaleMargCorrTplshape[k];
    StrongELPresentTplshapes[k][index_z] = strongEmissionLinePresentTplshape[k];
    StrongHalphaELPresentTplshapes[k][index_z] =
        strongHalphaELPresentTplshapes[k];
    NLinesAboveSNRTplshapes[k][index_z] = nLinesAboveSNRTplshape[k];
    PriorLinesTplshapes[k][index_z] = priorLinesTplshape[k];
  }
  return;
}

TFloat64List CLineModelResult::getChisquareTplContinuumResult(Int32 index_z) {

  TFloat64List ret;
  ret.reserve(ChiSquareTplContinuum.size());
  for (auto it = ChiSquareTplContinuum.begin(), e = ChiSquareTplContinuum.end();
       it != e; ++it)
    ret.push_back(it->at(index_z));

  return ret;
}

TFloat64List CLineModelResult::getChisquareTplshapeResult(Int32 index_z) {
  TFloat64List ret;
  ret.reserve(ChiSquareTplshapes.size());
  for (auto it = ChiSquareTplshapes.begin(), e = ChiSquareTplshapes.end();
       it != e; ++it)
    ret.push_back(it->at(index_z));

  return ret;
}

TFloat64List CLineModelResult::getScaleMargCorrTplshapeResult(Int32 index_z) {
  TFloat64List scaleMargCorrTplshape;
  if (index_z >= Redshifts.size()) {
    return scaleMargCorrTplshape;
  }
  if (ScaleMargCorrectionTplshapes.size() < 1) {
    return scaleMargCorrTplshape;
  }

  for (Int32 k = 0; k < ScaleMargCorrectionTplshapes.size(); k++) {
    scaleMargCorrTplshape.push_back(ScaleMargCorrectionTplshapes[k][index_z]);
  }

  return scaleMargCorrTplshape;
}

TBoolList CLineModelResult::getStrongELPresentTplshapeResult(Int32 index_z) {
  TBoolList _strongELPresentTplshape;
  if (index_z >= Redshifts.size() || StrongELPresentTplshapes.size() < 1) {
    return _strongELPresentTplshape;
  }

  for (Int32 k = 0; k < StrongELPresentTplshapes.size(); k++) {
    _strongELPresentTplshape.push_back(StrongELPresentTplshapes[k][index_z]);
  }

  return _strongELPresentTplshape;
}

TBoolList CLineModelResult::getHaELPresentTplshapeResult(Int32 index_z) {
  TBoolList _strongHaPresentTplshape;
  if (index_z >= Redshifts.size() ||
      StrongHalphaELPresentTplshapes.size() < 1) {
    return _strongHaPresentTplshape;
  }

  for (Int32 k = 0; k < StrongHalphaELPresentTplshapes.size(); k++) {
    _strongHaPresentTplshape.push_back(
        StrongHalphaELPresentTplshapes[k][index_z]);
  }

  return _strongHaPresentTplshape;
}

TInt32List CLineModelResult::getNLinesAboveSNRTplshapeResult(Int32 index_z) {
  TInt32List priorTplshape;
  if (index_z >= Redshifts.size()) {
    return priorTplshape;
  }
  if (NLinesAboveSNRTplshapes.size() < 1) {
    return priorTplshape;
  }

  for (Int32 k = 0; k < NLinesAboveSNRTplshapes.size(); k++) {
    priorTplshape.push_back(NLinesAboveSNRTplshapes[k][index_z]);
  }

  return priorTplshape;
}

TFloat64List CLineModelResult::getPriorLinesTplshapeResult(Int32 index_z) {
  TFloat64List priorTplshape;
  if (index_z >= Redshifts.size()) {
    return priorTplshape;
  }
  if (PriorLinesTplshapes.size() < 1) {
    return priorTplshape;
  }

  for (Int32 k = 0; k < PriorLinesTplshapes.size(); k++) {
    priorTplshape.push_back(PriorLinesTplshapes[k][index_z]);
  }

  return priorTplshape;
}

Int32 CLineModelResult::getNLinesOverCutThreshold(Int32 solutionIdx,
                                                  Float64 snrThres,
                                                  Float64 fitThres) const {
  Int32 nSol = 0;

  TInt32List indexesSols;
  for (Int32 j = 0; j < LineModelSolutions[solutionIdx].Amplitudes.size();
       j++) {
    // skip if already sol
    bool alreadysol = false;
    for (Int32 i = 0; i < indexesSols.size(); i++) {
      if (LineModelSolutions[solutionIdx].ElementId[j] == indexesSols[i]) {
        alreadysol = true;
        break;
      }
    }
    if (alreadysol) {
      continue;
    }
    if (!LineModelSolutions[solutionIdx].Lines[j].GetIsStrong()) {
      continue;
    }
    if (!LineModelSolutions[solutionIdx].Lines[j].GetIsEmission()) {
      continue;
    }

    Float64 noise = LineModelSolutions[solutionIdx].AmplitudesUncertainties[j];
    if (noise > 0) {
      Float64 snr = LineModelSolutions[solutionIdx].Amplitudes[j] / noise;
      Float64 Fittingsnr = LineModelSolutions[solutionIdx].Amplitudes[j] /
                           LineModelSolutions[solutionIdx].FittingError[j];
      if (snr >= snrThres && Fittingsnr >= fitThres) {
        nSol++;
        indexesSols.push_back(LineModelSolutions[solutionIdx].ElementId[j]);
      }
    }
  }
  return nSol;
}

/**
 * @brief CLineModelResult::getStrongLinesPresence
 * @param filterType: 1: emission only, 2 abs only, else: no filter
 * @return: a list of boolean values indicating if a strong is present (not
 * outsidelambdarange for that z) for each redshift
 */
TBoolList CLineModelResult::getStrongLinesPresence(
    Int32 filterType,
    const std::vector<CLineModelSolution> &linemodelsols) const {
  TBoolList strongIsPresent(linemodelsols.size(), false);
  for (Int32 solutionIdx = 0; solutionIdx < linemodelsols.size();
       solutionIdx++) {
    for (Int32 j = 0; j < linemodelsols[solutionIdx].Amplitudes.size(); j++) {
      if (!linemodelsols[solutionIdx].Lines[j].GetIsStrong()) {
        continue;
      }

      if (filterType == 1 &&
          !linemodelsols[solutionIdx].Lines[j].GetIsEmission())
        continue;
      else if (filterType == 2 &&
               linemodelsols[solutionIdx].Lines[j].GetIsEmission())
        continue;

      if (linemodelsols[solutionIdx].OutsideLambdaRange[j]) {
        continue;
      }

      if (linemodelsols[solutionIdx].Amplitudes[j] > 0.0) {
        strongIsPresent[solutionIdx] = true;
        break;
      }
    }
  }

  return strongIsPresent;
}

TInt32List CLineModelResult::getNLinesAboveSnrcut(
    const std::vector<CLineModelSolution> &linemodelsols) const {
  TInt32List nlinesabove(linemodelsols.size(), 0);
  for (Int32 solutionIdx = 0; solutionIdx < linemodelsols.size();
       solutionIdx++) {
    nlinesabove[solutionIdx] = linemodelsols[solutionIdx].NLinesAboveSnrCut;
  }

  return nlinesabove;
}

/**
 * WARNING: this function has not been tested at all !!! please check/debug
 * @brief CLineModelResult::getStrongestLineIsHa
 * @return: a list of boolean values indicating if the strongest line is Ha
 * (Highest amp and not outsidelambdarange for that z) for each redshift
 */
TBoolList CLineModelResult::getStrongestLineIsHa(
    const std::vector<CLineModelSolution> &linemodelsols) const {
  TBoolList isHaStrongest(linemodelsols.size(), false);
  std::string ampMaxLineTag = "";
  for (Int32 solutionIdx = 0; solutionIdx < linemodelsols.size();
       solutionIdx++) {
    Float64 ampMax = -DBL_MAX;
    ampMaxLineTag = "undefined";
    for (Int32 j = 0; j < linemodelsols[solutionIdx].Amplitudes.size(); j++) {
      if (!linemodelsols[solutionIdx].Lines[j].GetIsEmission() ||
          linemodelsols[solutionIdx].OutsideLambdaRange[j])
        continue;

      Log.LogDebug("    linemodelresult: using line for max amp search=%s",
                   linemodelsols[solutionIdx].Lines[j].GetName().c_str());
      if (linemodelsols[solutionIdx].Amplitudes[j] > ampMax) {
        ampMax = linemodelsols[solutionIdx].Amplitudes[j];
        ampMaxLineTag = linemodelsols[solutionIdx].Lines[j].GetName().c_str();
      }
    }

    isHaStrongest[solutionIdx] = (!std::isnan(ampMax) && ampMax > 0. &&
                                  ampMaxLineTag == linetags::halpha_em);
    if (isHaStrongest[solutionIdx]) {
      Log.LogDebug("CLineModelResult::GetModelHaStrongest:  z=%f found to be "
                   "true with ampMax=%e (for line=Halpha)",
                   linemodelsols[solutionIdx].Redshift, ampMax);
    }
  }

  return isHaStrongest;
}

Float64 CLineModelResult::getMinChiSquare() const {
  Float64 min = DBL_MAX;
  for (int i = 0; i < Redshifts.size(); i++) {
    if (min > ChiSquare[i]) {
      min = ChiSquare[i];
    }
  }
  return min;
}

Float64 CLineModelResult::getMaxChiSquare() const {
  Float64 max = -DBL_MAX;
  for (int i = 0; i < Redshifts.size(); i++) {
    if (max < ChiSquare[i]) {
      max = ChiSquare[i];
    }
  }
  return max;
}
