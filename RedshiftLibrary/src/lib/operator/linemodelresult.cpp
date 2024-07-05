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
#include <cfloat>
#include <fstream>
#include <string>

#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/common/vectorOperations.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/linemodel/templatesfitstore.h"
#include "RedshiftLibrary/linemodel/tplratiomanager.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/operator/linemodelresult.h"
#include "RedshiftLibrary/statistics/deltaz.h"

using namespace NSEpic;

/**
 * @brief CLineModelResult::Init
 * Initializes the linemodel results, by:
 *  - setting the redshift grid
 *  - setting the size of chisquare vectors
 * @param redshifts
 * @param restLines
 * @param nTplratios
 * @return
 */
void CLineModelResult::Init(TFloat64List redshifts, CLineMap restLines,
                            Int32 nTemplates, Int32 nTplratios,
                            TFloat64List tplratiosPriors) {
  if (tplratiosPriors.size() != nTplratios) {
    THROWG(
        ErrorCode::INTERNAL_ERROR,
        Formatter()
            << "Sizes do not match between tplratioprior and tplratios vectors:"
            << tplratiosPriors.size() << " vs " << nTplratios);
  }

  Int32 nResults = redshifts.size();
  ChiSquare.assign(nResults, NAN);
  ScaleMargCorrection.assign(nResults, NAN);
  Redshifts = std::move(redshifts);
  restLineList = std::move(restLines);
  LineModelSolutions.assign(nResults, CLineModelSolution());
  ContinuumModelSolutions.assign(nResults, CContinuumModelSolution());

  ChiSquareContinuum.assign(nResults, NAN);
  ScaleMargCorrectionContinuum.assign(nResults, NAN);

  if (nTemplates > 0)
    ChiSquareTplContinuum.assign(nTemplates, TFloat64List(nResults, DBL_MAX));
  else
    ChiSquareTplContinuum.clear();

  // init the tplratio chisquare results
  if (nTplratios > 0) {
    ChiSquareTplratios.assign(nTplratios, TFloat64List(nResults, DBL_MAX));
    ScaleMargCorrectionTplratios.assign(nTplratios,
                                        TFloat64List(nResults, 0.0));
    StrongELPresentTplratios.assign(nTplratios, TBoolList(nResults, false));
    StrongHalphaELPresentTplratios.assign(nTplratios,
                                          TBoolList(nResults, false));
    NLinesAboveSNRTplratios.assign(nTplratios, TInt32List(nResults, 0));
    PriorTplratios = std::move(tplratiosPriors);
    PriorLinesTplratios.assign(nTplratios, TFloat64List(nResults, 0.0));
  } else {
    ChiSquareTplratios.clear();
    ScaleMargCorrectionTplratios.clear();
    StrongELPresentTplratios.clear();
    NLinesAboveSNRTplratios.clear();
    PriorTplratios.clear();
    PriorLinesTplratios.clear();
  }
}

void CLineModelResult::updateVectors(Int32 idx, Int32 ndup, Int32 count) {

  insertWithDuplicates<Float64>(ChiSquare, idx, count, NAN, ndup);
  insertWithDuplicates<Float64>(ScaleMargCorrection, idx, count, NAN, ndup);
  insertWithDuplicates(LineModelSolutions, idx, count, CLineModelSolution(),
                       ndup);
  insertWithDuplicates(ContinuumModelSolutions, idx, count,
                       CContinuumModelSolution(), ndup);
  insertWithDuplicates<Float64>(ChiSquareContinuum, idx, count, NAN, ndup);
  insertWithDuplicates<Float64>(ScaleMargCorrectionContinuum, idx, count, NAN,
                                ndup);

  for (auto &xi2Cont : ChiSquareTplContinuum)
    insertWithDuplicates<Float64>(xi2Cont, idx, count, DBL_MAX, ndup);

  for (Int32 i = 0; i < ChiSquareTplratios.size(); i++) {
    insertWithDuplicates<Float64>(ChiSquareTplratios[i], idx, count, DBL_MAX,
                                  ndup);
    insertWithDuplicates<Float64>(ScaleMargCorrectionTplratios[i], idx, count,
                                  0., ndup);
    insertWithDuplicates<bool>(StrongELPresentTplratios[i], idx, count, false,
                               ndup);
    insertWithDuplicates<bool>(StrongHalphaELPresentTplratios[i], idx, count,
                               false, ndup);
    insertWithDuplicates<Int32>(NLinesAboveSNRTplratios[i], idx, count, 0,
                                ndup);
    insertWithDuplicates<Float64>(PriorLinesTplratios[i], idx, count, 0., ndup);
  }
}

void CLineModelResult::SetChisquareContinuumResult(
    Int32 index_z,
    const shared_ptr<const CContinuumFitStore> &continuumFitStore) {
  const auto index_z_in_store =
      continuumFitStore->GetRedshiftIndex(Redshifts[index_z]);
  if (index_z_in_store == -1)
    THROWG(ErrorCode::INTERNAL_ERROR, "Redshift not in fitstore");

  for (Int32 k = 0; k < continuumFitStore->GetContinuumCount();
       k++) { // TODO: handle the use of more than one continuum in linemodel
    ChiSquareTplContinuum[k][index_z] =
        continuumFitStore->GetFitValues(index_z_in_store, k).tplMerit;
  }
}

void CLineModelResult::SetChisquareContinuumResultFromPrevious(Int32 index_z) {
  auto previous = index_z - 1;
  for (auto it = ChiSquareTplContinuum.begin(), e = ChiSquareTplContinuum.end();
       it != e; ++it)
    it->at(index_z) = it->at(previous);
}

void CLineModelResult::SetChisquareTplratioResult(
    Int32 index_z, std::shared_ptr<CTplratioManager> tplratioManager) {
  if (tplratioManager->GetChisquareTplratio().size() < 1)
    return;

  if (index_z >= Redshifts.size())
    THROWG(ErrorCode::INTERNAL_ERROR, "Invalid z index");

  if (tplratioManager->GetChisquareTplratio().size() !=
          ChiSquareTplratios.size() ||
      tplratioManager->GetChisquareTplratio().size() !=
          tplratioManager->GetScaleMargTplratio().size() ||
      tplratioManager->GetChisquareTplratio().size() !=
          tplratioManager->GetStrongELPresentTplratio().size() ||
      tplratioManager->GetChisquareTplratio().size() !=
          tplratioManager->GetNLinesAboveSNRTplratio().size() ||
      tplratioManager->GetChisquareTplratio().size() !=
          tplratioManager->GetPriorLinesTplratio().size())
    THROWG(ErrorCode::INTERNAL_ERROR, "vector sizes do not match");

  for (Int32 k = 0; k < tplratioManager->GetChisquareTplratio().size(); k++) {
    ChiSquareTplratios[k][index_z] = tplratioManager->GetChisquareTplratio()[k];
    ScaleMargCorrectionTplratios[k][index_z] =
        tplratioManager->GetScaleMargTplratio()[k];
    StrongELPresentTplratios[k][index_z] =
        tplratioManager->GetStrongELPresentTplratio()[k];
    StrongHalphaELPresentTplratios[k][index_z] =
        tplratioManager->getHaELPresentTplratio()[k];
    NLinesAboveSNRTplratios[k][index_z] =
        tplratioManager->GetNLinesAboveSNRTplratio()[k];
    PriorLinesTplratios[k][index_z] =
        tplratioManager->GetPriorLinesTplratio()[k];
  }
  return;
}

void CLineModelResult::SetChisquareTplratioResultFromPrevious(Int32 index_z) {

  if (index_z >= Redshifts.size())
    THROWG(ErrorCode::INTERNAL_ERROR, "Invalid z index");

  auto previous = index_z - 1;

  for (Int32 k = 0; k < ChiSquareTplratios.size(); k++) {
    ChiSquareTplratios[k][index_z] = ChiSquareTplratios[k][previous];
    ScaleMargCorrectionTplratios[k][index_z] =
        ScaleMargCorrectionTplratios[k][previous];
    StrongELPresentTplratios[k][index_z] =
        StrongELPresentTplratios[k][previous];
    StrongHalphaELPresentTplratios[k][index_z] =
        StrongHalphaELPresentTplratios[k][previous];
    NLinesAboveSNRTplratios[k][index_z] = NLinesAboveSNRTplratios[k][previous];
    PriorLinesTplratios[k][index_z] = PriorLinesTplratios[k][previous];
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

TFloat64List CLineModelResult::getChisquareTplratioResult(Int32 index_z) {
  TFloat64List ret;
  ret.reserve(ChiSquareTplratios.size());
  for (auto it = ChiSquareTplratios.begin(), e = ChiSquareTplratios.end();
       it != e; ++it)
    ret.push_back(it->at(index_z));

  return ret;
}

TFloat64List CLineModelResult::getScaleMargCorrTplratioResult(Int32 index_z) {
  TFloat64List scaleMargCorrTplratio;
  if (index_z >= Redshifts.size()) {
    return scaleMargCorrTplratio;
  }
  if (ScaleMargCorrectionTplratios.size() < 1) {
    return scaleMargCorrTplratio;
  }

  for (Int32 k = 0; k < ScaleMargCorrectionTplratios.size(); k++) {
    scaleMargCorrTplratio.push_back(ScaleMargCorrectionTplratios[k][index_z]);
  }

  return scaleMargCorrTplratio;
}

TBoolList CLineModelResult::getStrongELPresentTplratioResult(Int32 index_z) {
  TBoolList _strongELPresentTplratio;
  if (index_z >= Redshifts.size() || StrongELPresentTplratios.size() < 1) {
    return _strongELPresentTplratio;
  }

  for (Int32 k = 0; k < StrongELPresentTplratios.size(); k++) {
    _strongELPresentTplratio.push_back(StrongELPresentTplratios[k][index_z]);
  }

  return _strongELPresentTplratio;
}

TBoolList CLineModelResult::getHaELPresentTplratioResult(Int32 index_z) {
  TBoolList _strongHaPresentTplratio;
  if (index_z >= Redshifts.size() ||
      StrongHalphaELPresentTplratios.size() < 1) {
    return _strongHaPresentTplratio;
  }

  for (Int32 k = 0; k < StrongHalphaELPresentTplratios.size(); k++) {
    _strongHaPresentTplratio.push_back(
        StrongHalphaELPresentTplratios[k][index_z]);
  }

  return _strongHaPresentTplratio;
}

TInt32List CLineModelResult::getNLinesAboveSNRTplratioResult(Int32 index_z) {
  TInt32List priorTplratio;
  if (index_z >= Redshifts.size()) {
    return priorTplratio;
  }
  if (NLinesAboveSNRTplratios.size() < 1) {
    return priorTplratio;
  }

  for (Int32 k = 0; k < NLinesAboveSNRTplratios.size(); k++) {
    priorTplratio.push_back(NLinesAboveSNRTplratios[k][index_z]);
  }

  return priorTplratio;
}

TFloat64List CLineModelResult::getPriorLinesTplratioResult(Int32 index_z) {
  TFloat64List priorTplratio;
  if (index_z >= Redshifts.size()) {
    return priorTplratio;
  }
  if (PriorLinesTplratios.size() < 1) {
    return priorTplratio;
  }

  for (Int32 k = 0; k < PriorLinesTplratios.size(); k++) {
    priorTplratio.push_back(PriorLinesTplratios[k][index_z]);
  }

  return priorTplratio;
}

Int32 CLineModelResult::getNLinesOverCutThreshold(Int32 solutionIdx,
                                                  Float64 snrThres,
                                                  Float64 fitThres) const {
  Int32 nSol = 0;
  TInt32List indexesSols;
  for (Int32 j = 0; j < LineModelSolutions[solutionIdx].Amplitudes.size();
       j++) {
    Int32 line_id = LineModelSolutions[solutionIdx].lineId[j];
    // skip if already sol
    bool alreadysol = std::find(indexesSols.begin(), indexesSols.end(),
                                LineModelSolutions[solutionIdx].ElementId[j]) !=
                      indexesSols.end();

    if (alreadysol || !restLineList.at(line_id).IsStrong() ||
        !restLineList.at(line_id).IsEmission())
      continue;

    Float64 const noise =
        LineModelSolutions[solutionIdx].AmplitudesUncertainties[j];
    if (noise <= 0)
      continue;
    Float64 const snr = LineModelSolutions[solutionIdx].Amplitudes[j] / noise;
    Float64 const Fittingsnr = LineModelSolutions[solutionIdx].Amplitudes[j] /
                               LineModelSolutions[solutionIdx].ResidualRMS[j];
    if (snr >= snrThres && Fittingsnr >= fitThres) {
      nSol++;
      indexesSols.push_back(LineModelSolutions[solutionIdx].ElementId[j]);
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

  auto lineDoesntMatchStrongFilter = [&filterType](const CLine &line) {
    return !line.IsStrong() ||
           // i have a doubt on these below conditions. didier ?
           (filterType == 1 && !line.IsEmission()) ||
           (filterType == 2 && line.IsEmission());
  };
  TBoolList strongIsPresent(linemodelsols.size(), false);
  for (Int32 solutionIdx = 0; solutionIdx < linemodelsols.size();
       solutionIdx++) {
    // search for the first strong line, with valid amplitude
    for (Int32 j = 0; j < linemodelsols[solutionIdx].Amplitudes.size(); j++) {
      Int32 line_id = linemodelsols[solutionIdx].lineId[j];
      if (lineDoesntMatchStrongFilter(restLineList.at(line_id)))
        continue;

      if (linemodelsols[solutionIdx].isLineValid(j)) {
        strongIsPresent[solutionIdx] = true;
        break;
      }
    }
  }

  return strongIsPresent;
}

TInt32List CLineModelResult::getNLinesAboveSnrcut(
    const std::vector<CLineModelSolution> &linemodelsols) const {
  TInt32List nlinesabove(linemodelsols.size());
  std::transform(
      linemodelsols.cbegin(), linemodelsols.cend(), nlinesabove.begin(),
      [](const CLineModelSolution &sol) { return sol.NLinesAboveSnrCut; });
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
  std::string ampMaxLineTag = undefStr;
  for (Int32 solutionIdx = 0; solutionIdx < linemodelsols.size();
       solutionIdx++) {
    Float64 ampMax = -DBL_MAX;
    ampMaxLineTag = undefStr;
    for (Int32 j = 0; j < linemodelsols[solutionIdx].Amplitudes.size(); j++) {
      Int32 line_id = linemodelsols[solutionIdx].lineId[j];
      if (!restLineList.at(line_id).IsEmission() ||
          linemodelsols[solutionIdx].OutsideLambdaRange[j])
        continue;

      Log.LogDebug(Formatter()
                   << "    linemodelresult: using line for max amp search="
                   << restLineList.at(line_id).GetName());
      if (linemodelsols[solutionIdx].Amplitudes[j] > ampMax) {
        ampMax = linemodelsols[solutionIdx].Amplitudes[j];
        ampMaxLineTag = restLineList.at(line_id).GetName().c_str();
      }
    }

    isHaStrongest[solutionIdx] = (!std::isnan(ampMax) && ampMax > 0. &&
                                  ampMaxLineTag == linetags::halpha_em);
    if (isHaStrongest[solutionIdx]) {
      Log.LogDebug(Formatter() << "CLineModelResult::GetModelHaStrongest:  z="
                               << linemodelsols[solutionIdx].Redshift
                               << " found to be true with ampMax=" << ampMax
                               << " (for line=Halpha)");
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
