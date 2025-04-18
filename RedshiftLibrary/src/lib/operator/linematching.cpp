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
#include <numeric>

#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/operator/linematching.h"

using namespace NSEpic;

/**
 * Executes the algorithm for matching the input sets of lines against the
 * catalogue of atomic transitions. Parameters:
 * - detectedLineCatalog, the set of lines detected in the spectrum.
 * - restLineCatalog, the set of reference lines.
 * - redshiftRange.
 * - nThreshold.
 * - tol.
 * - typeFilter.
 * - detectedForceFilter.
 * - restForceFilter.
 * Returns a smart pointer of a CLineMatchingResult, containing the set of
 * matches.
 */
std::shared_ptr<CLineMatchingResult> CLineMatching::Compute(
    const CLineDetectedCatalog &detectedLineCatalog,
    const CLineCatalog &restLineCatalog, const TFloat64Range &redshiftRange,
    Int32 nThreshold, Float64 tol, CLine::EType typeFilter,
    CLine::EForce detectedForceFilter, CLine::EForce restForceFilter) const {
  auto const detectedLineList =
      detectedLineCatalog.GetFilteredList(typeFilter, detectedForceFilter);
  auto const restLineList =
      restLineCatalog.GetFilteredList(typeFilter, restForceFilter);
  auto const redshiftGetter = getRedshift();

  CLineMatchingResult::TSolutionSetList solutions;
  Int32 const N = restLineList.size();
  Log.LogDebug(Formatter() << "CLineMatching: Lines detected. n=" << N);

  for (auto const &[iDetectedLine, DetectedLine] : detectedLineList) {
    for (auto const &[iRestLine, RestLine] : restLineList) {
      // for each detected line / rest line couple, enumerate how many other
      // couples fit within the tolerance
      CLineMatchingResult::TSolutionSet solution;
      Float64 redShift = redshiftGetter(DetectedLine, RestLine);
      if (redShift < 0) // we don't care about blueshifts.
        continue;
      solution.push_back(CLineMatchingResult::SSolution(
          DetectedLine, {iRestLine, RestLine}, redShift));

      if (N == 1) {
        // if there is only one detected line, then there are N=#restlines
        // possible redshifts
        solutions.push_back(solution);
        break;
      }
      updateSolution(iDetectedLine, redShift, tol, solution, // to update
                     detectedLineList, restLineList);
      if (solution.size() > 0)
        solutions.push_back(solution);
    }
  }

  Log.LogDebug(Formatter() << "CLineMatching: non unique solutions found n="
                           << solutions.size());
  const CLineMatchingResult::TSolutionSetList newSolutions =
      refineSolutions(solutions, redshiftRange, nThreshold);

  if (!newSolutions.size())
    return nullptr;

  auto result = std::make_shared<CLineMatchingResult>();
  result->SolutionSetList = newSolutions;
  result->m_RestCatalog = restLineCatalog;
  result->m_DetectedCatalog = detectedLineCatalog;
  return result;
}

std::function<Float64(CLineDetected const &, CLine const &)>
CLineMatching::getRedshift() const {
  return [](CLineDetected const &detectedLine, CLine const &restLine) {
    return (detectedLine.GetPosition() - restLine.GetPosition()) /
           restLine.GetPosition();
  };
}

void CLineMatching::updateSolution(
    Int32 iDetectedLine, Float64 redShift, Float64 tol,
    CLineMatchingResult::TSolutionSet &solution, // to update
    CLineDetectedMap const &detectedLineList,
    CLineMap const &restLineList) const {
  auto redshiftGetter = getRedshift();

  for (auto const &[iDetectedLine2, DetectedLine2] : detectedLineList) {
    if (iDetectedLine == iDetectedLine2)
      continue;

    for (auto const &[iRestLine2, RestLine2] : restLineList) {
      Float64 redShift2 = redshiftGetter(DetectedLine2, RestLine2);
      if (redShift2 < 0) // we don't care about blueshifts.
        continue;
      Float64 redshiftTolerance = tol * (1 + (redShift + redShift2) * 0.5);
      if (fabs((redShift - redShift2)) > redshiftTolerance)
        continue;

      // avoid repeated solution sets
      if (!isLineAlreadyPresent(DetectedLine2, solution)) {
        solution.push_back(CLineMatchingResult::SSolution(
            DetectedLine2, {iRestLine2, RestLine2}, redShift2));
      }
    }
  }
  return;
}

/**
 * @brief Check if line is already present in the solution
 *
 * @param line
 * @param solution
 * @return true
 * @return false
 */
bool CLineMatching::isLineAlreadyPresent(
    const CLineDetected &line,
    const CLineMatchingResult::TSolutionSet &solution) const {
  for (const auto &currentSet : solution)
    if (currentSet.DetectedLine.GetPosition() == line.GetPosition())
      return true;

  return false;
}

// delete duplicated solutions
const CLineMatchingResult::TSolutionSetList CLineMatching::refineSolutions(
    const CLineMatchingResult::TSolutionSetList &solutions,
    const TFloat64Range &redshiftRange, Int32 nThreshold) const {

  CLineMatchingResult::TSolutionSetList newSolutions;
  for (auto currentSet : solutions) {
    if (ssize(currentSet) < nThreshold)
      continue;
    std::sort(currentSet.begin(), currentSet.end());
    bool found = false;
    for (auto newsol : newSolutions)
      if (AreSolutionSetsEqual(newsol, currentSet)) {
        found = true;
        break;
      }
    if (found) // found duplicate
      continue;
    // if solution is new, let's check it is within the range
    Float64 redshiftMean =
        accumulate(currentSet.begin(), currentSet.end(), 0.,
                   [](Float64 v, const CLineMatchingResult::SSolution &s2) {
                     return v + s2.Redshift;
                   }) /
        currentSet.size();

    if (redshiftMean > redshiftRange.GetBegin() &&
        redshiftMean < redshiftRange.GetEnd())
      newSolutions.push_back(currentSet);
  }
  Log.LogDebug(Formatter() << "CLineMatching: unique solutions found n="
                           << newSolutions.size());
  return newSolutions;
}
/**
 * Given 2 TSolutionSets, returns true if they are equivalent (have the same
 * line values), false otherwise. The algorithm only works reliably if the
 * inputs are sorted.
 */
bool CLineMatching::AreSolutionSetsEqual(
    const CLineMatchingResult::TSolutionSet &s1,
    const CLineMatchingResult::TSolutionSet &s2) const {
  if (s1.size() != s2.size())
    return false;

  bool diffFound = false;
  for (Int32 iSet = 0; iSet < ssize(s1); iSet++) {
    if (s1[iSet].DetectedLine != s2[iSet].DetectedLine ||
        s1[iSet].RestLineId != s2[iSet].RestLineId ||
        s1[iSet].Redshift != s2[iSet].Redshift) {
      diffFound = true;
      break;
    }
  }

  return !diffFound;
}
