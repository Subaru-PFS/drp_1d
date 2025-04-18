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
#include <climits>
#include <numeric>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/log/log.h"

using namespace NSEpic;
using namespace std;

/**
 * Member attribution constructor.
 */
CExtremum::CExtremum(Int32 maxPeakCount, Float64 peakSeparation,
                     Float64 meritcut, bool invertForMinSearch,
                     bool allow_extrema_at_border, const TFloat64Range &xRange)
    : m_MaxPeakCount(maxPeakCount), m_extrema_separation(peakSeparation),
      m_meritCut(meritcut), m_SignSearch(invertForMinSearch ? -1.0 : 1.0),
      m_allow_extrema_at_border(allow_extrema_at_border),
      m_PeakSeparationActive(peakSeparation > 0), m_XRange(xRange) {}

/**
 * create an index vector sorting maxY elements
 */
void CExtremum::SortIndexes(TFloat64List const &maxY) const {
  m_sortedIndexes.resize(maxY.size());
  iota(m_sortedIndexes.begin(), m_sortedIndexes.end(), 0);
  sort(m_sortedIndexes.begin(), m_sortedIndexes.end(),
       [&maxY](Int32 i1, Int32 i2) { return maxY[i1] > maxY[i2]; });
}

/**
 * Wrapper around InternalFind, this function validates and sets the search
 * interval limits.
 */
TPointList CExtremum::Find(const TFloat64List &xAxis, const TFloat64List &yAxis,
                           const bool isFirstPass) const {
  Int32 n = xAxis.size();

  if (n == 0)
    THROWG(ErrorCode::INTERNAL_ERROR, "input X vector is empty");

  if (n != ssize(yAxis))
    THROWG(ErrorCode::INTERNAL_ERROR,
           "input X and Y vector do not have the same size");

  TPointList maxPoint;

  Int32 BeginIndex = 0;
  Int32 EndIndex = n - 1;

  // Find index in xAxis that correspond to the boundary specified by m_XRange
  if (!m_XRange.GetIsEmpty())
    m_XRange.getClosedIntervalIndices(xAxis, BeginIndex, EndIndex);

  auto [maxX, maxY] = FindAllPeaks(xAxis, yAxis, BeginIndex, EndIndex);

  if (maxX.size() == 0) {
    if (isFirstPass) {
      THROWG(ErrorCode::PDF_PEAK_NOT_FOUND,
             Formatter() << "CExtremum::" << __func__
                         << " Failed to identify pdf candidates ");
    } else {
      // we should not raise an exception here
      // (we can accept a missing candidate in second pass window)
      Flag.warning(WarningCode::PDF_PEAK_NOT_FOUND,
                   Formatter() << "          CExtremum::" << __func__
                               << ": FindAllPeaks returned empty MaxX");
      return maxPoint;
    }
  }

  // TODO: add a boolean referring to the metric to use for the sort
  // by default, using the pdf value
  SortIndexes(maxY);

  Int32 keepMinN = 1;
  if (m_meritCut > 0.0 && ssize(maxX) > keepMinN) {
    Cut_Threshold(maxX, maxY, keepMinN);
  }

  // refine: eliminate very close candidates when possible.
  if (m_PeakSeparationActive) {
    maxPoint = FilterOutNeighboringPeaksAndTruncate(maxX, maxY, keepMinN);

    // verify that peaks are well separated by at least secondpassradius
    assertPeakSeparation(maxPoint);
  } else {
    maxPoint = Truncate(maxX, maxY);
  }

  return maxPoint;
}

void CExtremum::assertPeakSeparation(TFloat64List &maxX) const {

  std::sort(maxX.begin(), maxX.end());
  for (Int32 i = 0; i < ssize(maxX) - 1; i++) {
    TFloat64Range windowh(
        maxX[i + 1] - m_extrema_separation / 2 * (1 + maxX[i + 1]),
        maxX[i + 1] + m_extrema_separation / 2 * (1 + maxX[i + 1]));
    windowh.IntersectWith(TFloat64Range(maxX));
    TFloat64Range windowl(maxX[i] - m_extrema_separation / 2 * (1 + maxX[i]),
                          maxX[i] + m_extrema_separation / 2 * (1 + maxX[i]));
    windowl.IntersectWith(TFloat64Range(maxX));
    Float64 overlap;
    overlap = windowh.GetBegin() - windowl.GetEnd();
    if (overlap < 0)
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << " Peaks " << maxX[i] << " and " << maxX[i + 1]
                         << " are not enough separated.");
  }
}

void CExtremum::assertPeakSeparation(TPointList &maxPoint) const {
  TFloat64List maxX(maxPoint.size());
  for (Int32 i = 0; i < ssize(maxPoint); i++)
    maxX[i] = maxPoint[i].X;
  assertPeakSeparation(maxX);
}

/**
 * \Brief: removes extrema based on meritCut.
 * First: Sort based on Y values
 * Second: Remove non-relevant candidates
 * Third: Sort based on X values to return candidates following their initial
 * order It returns the new list of extrema
 */
void CExtremum::Cut_Threshold(TFloat64List &maxX, TFloat64List &maxY,
                              Int32 keepMinN) const {
  if (maxX.size() == 0)
    THROWG(ErrorCode::INTERNAL_ERROR, " empty MaxX arg");

  if (ssize(maxX) <= keepMinN)
    return;

  TFloat64List newX;
  TFloat64List newY;
  for (Int32 i = 0; i < ssize(maxX); i++) {
    Float64 meritDiff = maxY[m_sortedIndexes[0]] - maxY[i];
    if (meritDiff <= m_meritCut) {
      newX.push_back(maxX[i]);
      newY.push_back(maxY[i]);
    }
  }
  if (ssize(newX) < keepMinN)
    for (Int32 isort = ssize(newX); isort < keepMinN; isort++) {
      Int32 i = m_sortedIndexes[isort];
      auto it = std::lower_bound(newX.begin(), newX.end(), maxX[i]);
      Int32 index = it - newX.begin();
      newX.insert(it, maxX[i]);
      newY.insert(newY.begin() + index, maxY[i]);
    }

  maxX = std::move(newX);
  maxY = std::move(newY);
}

TPointList CExtremum::FilterOutNeighboringPeaksAndTruncate(
    TFloat64List const &maxX, TFloat64List const &maxY, Int32 keepmin) const {

  TPointList maxPoint;

  if (ssize(maxX) <= keepmin) {
    for (Int32 i = 0; i < ssize(maxX); i++) {
      maxPoint.push_back(SPoint(maxX[m_sortedIndexes[i]],
                                m_SignSearch * maxY[m_sortedIndexes[i]]));
    }
    return maxPoint;
  }

  TBoolList peakTracking(maxX.size(), true);
  Float64 wind_high, wind_low;
  Int32 nkeep = 0;

  for (Int32 i : m_sortedIndexes) {
    if (nkeep == m_MaxPeakCount) {
      break;
    }

    if (!peakTracking[i]) {
      continue;
    }
    nkeep++;
    maxPoint.push_back(SPoint(maxX[i], m_SignSearch * maxY[i]));
    wind_high = maxX[i] + (maxX[i] + 1) * m_extrema_separation /
                              (1 - m_extrema_separation / 2);
    wind_low = maxX[i] - (maxX[i] + 1) * m_extrema_separation /
                             (1 + m_extrema_separation / 2);
    TFloat64Range window(wind_low, wind_high);
    Int32 i_min = i, i_max = i;
    try {
      window.getClosedIntervalIndices(maxX, i_min, i_max);
    } catch (const AmzException &exception) {
    }
    for (Int32 j = i_min; j <= i_max; j++) {
      if (j == i)
        continue;
      peakTracking[j] = false;
    }
  }
  return maxPoint;
}

// find first and last non Nan element
void CExtremum::getFirstandLastnonNANElementIndices(const TFloat64List &yAxis,
                                                    Int32 &BeginIndex,
                                                    Int32 &EndIndex) const {
  // find first and last non Nan element
  for (; BeginIndex != EndIndex; ++BeginIndex)
    if (!std::isnan(yAxis[BeginIndex]))
      break;

  for (; BeginIndex != EndIndex; --EndIndex)
    if (!std::isnan(yAxis[EndIndex]))
      break;

  // check at least 3 points left (to get an extrema)
  if ((EndIndex - BeginIndex) < 2)
    THROWG(ErrorCode::INTERNAL_ERROR, "less than 3 contiguous non nan values");
}

TFloat64List CExtremum::applySign(const TFloat64List &yAxis,
                                  bool invertSearch) const {
  const Float64 SignSearch = invertSearch ? -1.0 * m_SignSearch : m_SignSearch;
  TFloat64List tmpY(yAxis.size());
  std::transform(yAxis.begin(), yAxis.end(), tmpY.begin(),
                 [SignSearch](const Float64 &val) { return val * SignSearch; });
  return tmpY;
}

std::pair<TFloat64List, TFloat64List>
CExtremum::FindAllPeaks(const TFloat64List &xAxis, const TFloat64List &yAxis,
                        Int32 BeginIndex, Int32 EndIndex,
                        bool invertSearch) const {
  TFloat64List tmpY = applySign(yAxis, invertSearch);

  auto goup = [&tmpY](Int32 i) { return tmpY[i] < tmpY[i + 1]; };
  auto godown = [&tmpY](Int32 i) { return tmpY[i] > tmpY[i + 1]; };
  auto goflat = [&tmpY](Int32 i) { return tmpY[i] == tmpY[i + 1]; };

  auto isPeak = [&tmpY, goup, godown](Int32 i) {
    return goup(i - 1) *
           godown(i); // using multiplication rather than && to reduce cc
  };
  auto isStartHighPlank = [&tmpY, goup, goflat](Int32 i) {
    return goup(i - 1) * goflat(i);
  };
  auto isPeakAfterPlank = [&tmpY, godown, goflat](Int32 i, bool plank) {
    return plank * goflat(i - 1) * godown(i);
  };
  auto isEndPlank = [&tmpY, goup, goflat](Int32 i) {
    return goflat(i - 1) * goup(i);
  };

  getFirstandLastnonNANElementIndices(yAxis, BeginIndex, EndIndex);

  Int32 maxCount = 0;
  bool plank = false;
  Int32 cnt_plk = 0;

  auto allow_extrema_border =
      m_allow_extrema_at_border; // to not have to pass this to lambda
  auto isLowBorderPeak = [&tmpY, allow_extrema_border,
                          godown](Int32 BeginIndex) {
    return allow_extrema_border * godown(BeginIndex);
  };
  auto isLowBorderPlank = [&tmpY, allow_extrema_border,
                           goflat](Int32 BeginIndex) {
    return allow_extrema_border * goflat(BeginIndex);
  };

  auto result = std::make_pair(TFloat64List(), TFloat64List());
  auto &[maxX, maxY] = result;

  // First element
  if (isLowBorderPeak(BeginIndex)) {
    maxX.push_back(xAxis[BeginIndex]);
    maxY.push_back(tmpY[BeginIndex]);
    maxCount++;
  } else if (isLowBorderPlank(BeginIndex)) {
    plank = true;
    cnt_plk++;
  }

  for (Int32 i = BeginIndex + 1; i < EndIndex; i++) {
    if (isPeak(i)) {
      maxX.push_back(xAxis[i]);
      maxY.push_back(tmpY[i]);
      maxCount++;
      continue;
    }
    if (isStartHighPlank(i)) {
      // plank start point
      plank = true;
      cnt_plk++;
      continue;
    }
    if (!goflat(i - 1))
      continue;

    // high plank: signal is decreasing after the plank. Peak identified
    if (isPeakAfterPlank(i, plank)) {
      // check if we already identified a high plank
      cnt_plk++;
      Int32 idx_plk = i - round(cnt_plk / 2);
      // plank end point
      maxX.push_back(xAxis[idx_plk]);
      maxY.push_back(tmpY[idx_plk]);
      maxCount++;
      plank = false;
      cnt_plk = 0;
      continue;
    }
    // low plank: pdf is increasing after the plank. No peak here!
    if (isEndPlank(i)) {
      plank = false; // end
      cnt_plk = 0;
      continue;
    }
    // Plank is getting larger
    cnt_plk++;
  }

  if (!m_allow_extrema_at_border)
    return result;
  // last element: check if plank is extended
  if (godown(EndIndex - 1))
    return result;

  cnt_plk += Int32(plank); // if no plank, cnt_plk+=0
  Int32 idx = EndIndex - round(cnt_plk / 2);

  maxX.push_back(xAxis[idx]);
  maxY.push_back(tmpY[idx]);
  maxCount++;

  return result;
}

/**
 * Brief: Reduce number of peaks based on maxCount passed from param.json if
 * present
 */
TPointList CExtremum::Truncate(TFloat64List const &maxX,
                               TFloat64List const &maxY) const {
  TPointList maxPoint;
  Int32 n = maxX.size();
  n = std::min(m_MaxPeakCount, n);
  for (Int32 j = 0; j < n; j++) {
    maxPoint.push_back(SPoint(maxX[m_sortedIndexes[j]],
                              m_SignSearch * maxY[m_sortedIndexes[j]]));
  }
  return maxPoint;
}
