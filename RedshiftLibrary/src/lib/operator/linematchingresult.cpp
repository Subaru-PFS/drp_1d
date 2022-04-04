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
#include "RedshiftLibrary/operator/linematchingresult.h"
#include "RedshiftLibrary/line/rules.h"
#include "RedshiftLibrary/log/log.h"
#include <algorithm> // std::sort
#include <math.h>
#include <stdio.h>

using namespace NSEpic;

/**
 * Attribution constructor. Will only set the debugging bypass option.
 */
CLineMatchingResult::CLineMatchingResult() { m_bypassDebug = true; }

/**
 * Empty destructor.
 */
CLineMatchingResult::~CLineMatchingResult() {}

/**
 * @brief CLineMatchingResult::GetBestRedshift
 *
 * select the best redshift given a set of rules:
 * - the solution with the most strong lines
 * - (TODO) if no solution with strong lines, the missing strong lines' absence
 * should be excused by:
 *      * high noise in the theoretical position of the strong lines
 *      * lambda range not able to cover the strong lines
 * @return
 */
bool CLineMatchingResult::GetBestRedshift(Float64 &Redshift,
                                          Int32 &MatchingNumber) const {
  Int32 thresMatchingNumber = 2; // minimum matching number for a solution
  TSolutionSetList selectedResults =
      GetSolutionsListOverNumber(thresMatchingNumber - 1);

  if (selectedResults.size() > 0) {
    int iStrongMax = -1;
    int nStrongMax = -1;
    for (Int32 iSol = 0; iSol < selectedResults.size(); iSol++) {
      int currentNStrongRestLines = getNStrongRestLines(selectedResults[iSol]);
      if (currentNStrongRestLines > nStrongMax) {
        iStrongMax = iSol;
        nStrongMax = currentNStrongRestLines;
      } else {
        if (currentNStrongRestLines == nStrongMax &&
            selectedResults[iSol].size() > selectedResults[iStrongMax].size()) {
          iStrongMax = iSol;
          nStrongMax = currentNStrongRestLines;
        }
      }
    }

    if (iStrongMax > -1) {
      Redshift = GetMeanRedshiftSolution(selectedResults[iStrongMax]);
      MatchingNumber = selectedResults[iStrongMax].size();
    }
  } else {
    return false;
  }
  return true;
}

/**
 * @brief CLineMatchingResult::GetBestMatchNumRedshift
 *
 * get the best redshift in the sense of the highest matching number
 *
 * @return
 */
bool CLineMatchingResult::GetBestMatchNumRedshift(Float64 &Redshift,
                                                  Int32 &MatchingNumber) const {
  MatchingNumber = GetMaxMatchingNumber();
  TSolutionSetList selectedResults =
      GetSolutionsListOverNumber(MatchingNumber - 1);

  if (selectedResults.size() > 0) {
    Redshift = GetMeanRedshiftSolution(selectedResults[0]);
  } else {
    return false;
  }
  return true;
}

/**
 * Returns the set of results that have more than "number" entries.
 */
CLineMatchingResult::TSolutionSetList
CLineMatchingResult::GetSolutionsListOverNumber(Int32 number) const {
  // select results by matching number
  TSolutionSetList selectedResults;
  for (Int32 iSol = 0; iSol < SolutionSetList.size(); iSol++) {
    TSolutionSet currentSet = SolutionSetList[iSol];
    if (currentSet.size() > number) {
      selectedResults.push_back(currentSet);
    }
  }
  return selectedResults;
}

/**
 * Returns a list of the mean redshifts for solutions that have at least
 * "number" lines.
 */
TFloat64List
CLineMatchingResult::GetAverageRedshiftListOverNumber(Int32 number) const {
  TFloat64List selectedRedshift;
  CLineMatchingResult::TSolutionSetList selectedResults =
      GetSolutionsListOverNumber(number);
  for (Int32 j = 0; j < selectedResults.size(); j++) {
    Float64 z = GetMeanRedshiftSolution(selectedResults[j]);
    selectedRedshift.push_back(z);
  }
  return selectedRedshift;
}

/**
 * Returns a custom-computed round value for the mean redshift for each solution
 * that has over "number" lines.
 */
TFloat64List CLineMatchingResult::GetRoundedRedshiftCandidatesOverNumber(
    Int32 number, Float64 step) const {
  TFloat64List selectedRedshift = GetAverageRedshiftListOverNumber(number);
  TFloat64List roundedRedshift;
  for (Int32 j = 0; j < selectedRedshift.size(); j++) {
    Float64 zround = Float64(int(selectedRedshift[j] / step + 0.5f) * step);
    roundedRedshift.push_back(zround);
  }
  return roundedRedshift;
}

/**
 * Retuns a list of redshifts that are "step"pedly within the "rangeWidth"
 * interval around each mean redshift for solutions with over "number" lines.
 */
TFloat64List CLineMatchingResult::GetExtendedRedshiftCandidatesOverNumber(
    Int32 number, Float64 step, Float64 rangeWidth) const {
  TFloat64List roundedRedshift =
      GetRoundedRedshiftCandidatesOverNumber(number, step);
  TFloat64List extendedRedshifts;
  Int32 halfk = rangeWidth / step / 2.0;
  for (Int32 j = 0; j < roundedRedshift.size(); j++) {
    for (Int32 k = -halfk; k < halfk; k++) {
      Float64 z = roundedRedshift[j] + k * step;
      extendedRedshifts.push_back(z);
    }
  }
  std::sort(extendedRedshifts.begin(), extendedRedshifts.end());
  extendedRedshifts.erase(
      std::unique(extendedRedshifts.begin(), extendedRedshifts.end()),
      extendedRedshifts.end());
  return extendedRedshifts;
}

/**
 * If "index" is valid, return the mean of the redshifts in the solution set at
 * "index".
 */
Float64 CLineMatchingResult::GetMeanRedshiftSolutionByIndex(Int32 index) const {
  if (index > SolutionSetList.size() - 1) {
    return -1.0;
  }

  Float64 redshiftMean = 0.0;
  const TSolutionSet &currentSet = SolutionSetList[index];
  for (Int32 i = 0; i < currentSet.size(); i++) {
    redshiftMean += currentSet[i].Redshift;
  }
  redshiftMean /= currentSet.size();
  return redshiftMean;
}

/**
 * Returns the mean redshift for the solution set "s".
 */
Float64
CLineMatchingResult::GetMeanRedshiftSolution(const TSolutionSet &s) const {
  TSolutionSet currentSet = s;
  Float64 redshiftMean = 0.0;
  for (Int32 i = 0; i < currentSet.size(); i++) {
    redshiftMean += currentSet[i].Redshift;
  }
  redshiftMean /= (float)currentSet.size();

  return redshiftMean;
}

/**
 * Returns the number of lines that are near (tolerance = 0.11 Angstroms) to
 * strong emission lines.
 */
Int32 CLineMatchingResult::getNStrongRestLines(const TSolutionSet &s) const {
  CLineCatalog::TLineVector strongRestLineList = m_RestCatalog.GetFilteredList(
      CLine::nType_Emission, CLine::nForce_Strong);
  Int32 ncatalog = strongRestLineList.size();

  TSolutionSet currentSet = s;
  Int32 nStrong = 0;
  Float64 tol = 0.11;
  for (Int32 i = 0; i < currentSet.size(); i++) {
    Int32 found = 0;
    for (Int32 c = 0; c < ncatalog; c++) {
      if (fabs(currentSet[i].RestLine.GetPosition() -
               strongRestLineList[c].GetPosition()) < tol) {
        found = 1;
        break;
      }
    }
    if (found == 1) {
      nStrong++;
    }
  }

  return nStrong;
}

/**
 * Returns the largest number of lines in a single solution for the current set
 * of solutions.
 */
Int32 CLineMatchingResult::GetMaxMatchingNumber() const {
  if (SolutionSetList.size() < 1) {
    return -1.0;
  }

  Int32 maxNumber = 0;
  for (Int32 i = 0; i < SolutionSetList.size(); i++) {
    TSolutionSet currentSet = SolutionSetList[i];
    if (maxNumber < currentSet.size()) {
      maxNumber = currentSet.size();
    }
  }
  return maxNumber;
}

/**
 * Attempts to apply all rules to each solution. If a rule is applicable, that
 * solution is removed from the current set of solutions.
 */
void CLineMatchingResult::FilterWithRules(CSpectrum spc,
                                          TFloat64Range lambdaRange,
                                          Float64 winsize) {
  Log.LogDebug("void CLineMatchingResult::FilterWithRules( CSpectrum spc, "
               "TFloat64Range lambdaRange, Float64 winsize )");
  if (SolutionSetList.size() < 1) {
    Log.LogDebug("SolutionSetList.size()<1, returning");
    return;
  }

  CRules rules(spc, m_DetectedCatalog, m_RestCatalog, lambdaRange, winsize);
  TSolutionSetList _solutionSetList;
  Log.LogDebug("There are %d solutions to test.", SolutionSetList.size());
  for (Int32 i = 0; i < SolutionSetList.size(); i++) {
    TSolutionSet currentSet = SolutionSetList[i];
    Float64 z = GetMeanRedshiftSolution(currentSet);
    Int32 ruleId = rules.check(z, currentSet);

    if (ruleId <= 0) {
      if (!m_bypassDebug)
        Log.LogDebug("Solution %d passed!", i);
      _solutionSetList.push_back(currentSet);
    } else {
      if (!m_bypassDebug)
        Log.LogDebug("Solution %d failed rule %d", i, ruleId);
      FilteredSolutionSetList.push_back(currentSet);
      FilterTypeList.push_back(ruleId);
    }
  }
  SolutionSetList = _solutionSetList;
}
