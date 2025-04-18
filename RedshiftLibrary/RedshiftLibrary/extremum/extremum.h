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
#ifndef _REDSHIFT_EXTREMUM_EXTREMUM_
#define _REDSHIFT_EXTREMUM_EXTREMUM_

#include <iostream>
#include <vector>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"

using namespace std;

namespace Extremum_test {
class Extremum_cut_isolated;
class Extremum_cut_onePeak;
class Extremum_FilterOutNeighboringPeaks;
class Extremum_FilterOutNeighboringPeaks_consecutifClosePeaks;
class Extremum_FilterOutNeighboringPeaks_consecutifClosePeaks2;
class Extremum_FilterOutNeighboringPeaks_increasingPDF;
class Extremum_FilterOutNeighboringPeaks2;
class Extremum_FilterOutNeighboringPeaks_negativeX;
class Extremum_Truncate;
} // namespace Extremum_test
namespace NSEpic {

/**
 * \ingroup Redshift
 * Analyse a given input 2D array and find the n stronger extremum.
 */
class CExtremum {

public:
  CExtremum(Int32 maxPeakCount = 10, Float64 peakSeparation = 0.005 * 2,
            Float64 meritcut = -1, // <0 = no cut thresholding
            bool invertForMinSearch = false,
            bool allow_extrema_at_border = true,
            const TFloat64Range &xRange = TFloat64Range());

  void SetMaxPeakCount(Int32 n) { m_MaxPeakCount = n; };
  void SetXRange(const TFloat64Range &r) { m_XRange = r; };
  void SetMeritCut(Float64 n) { m_meritCut = n; };
  TPointList Find(const TFloat64List &xAxis, const TFloat64List &yAxis,
                  const bool isFirstpass) const;

private:
  // for unit testing private methods
  friend class Extremum_test::Extremum_cut_isolated;
  friend class Extremum_test::Extremum_cut_onePeak;
  friend class Extremum_test::Extremum_FilterOutNeighboringPeaks;
  friend class Extremum_test::
      Extremum_FilterOutNeighboringPeaks_consecutifClosePeaks;
  friend class Extremum_test::
      Extremum_FilterOutNeighboringPeaks_consecutifClosePeaks2;
  friend class Extremum_test::Extremum_FilterOutNeighboringPeaks_increasingPDF;
  friend class Extremum_test::Extremum_FilterOutNeighboringPeaks2;
  friend class Extremum_test::Extremum_FilterOutNeighboringPeaks_negativeX;
  friend class Extremum_test::Extremum_Truncate;

  std::pair<TFloat64List, TFloat64List>
  FindAllPeaks(const TFloat64List &xAxis, const TFloat64List &yAxis,
               Int32 BeginIndex, Int32 EndIndex,
               bool invertSearch = false) const;

  TPointList FilterOutNeighboringPeaksAndTruncate(TFloat64List const &maxX,
                                                  TFloat64List const &maxY,
                                                  Int32 keepmin) const;

  TPointList Truncate(TFloat64List const &xAxis,
                      TFloat64List const &yAxis) const;

  void Cut_Threshold(TFloat64List &maxX, TFloat64List &maxY,
                     Int32 keepMinN) const;

  void SortIndexes(TFloat64List const &maxY) const;

  void assertPeakSeparation(TFloat64List &maxX) const;
  void assertPeakSeparation(TPointList &maxPoint) const;
  void getFirstandLastnonNANElementIndices(const TFloat64List &yAxis,
                                           Int32 &BeginIndex,
                                           Int32 &EndIndex) const;
  TFloat64List applySign(const TFloat64List &yAxis, bool invertSearch) const;

  Int32 m_MaxPeakCount;
  Float64 m_extrema_separation;
  Float64 m_meritCut;
  const Float64 m_SignSearch;
  bool m_allow_extrema_at_border = true;
  bool m_PeakSeparationActive = true;
  TFloat64Range m_XRange;
  mutable TInt32List m_sortedIndexes;
};

} // namespace NSEpic

#endif
