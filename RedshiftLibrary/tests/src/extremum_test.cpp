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
#include <cmath>
#include <iostream>
#include <numeric>

#include <boost/test/unit_test.hpp>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/extremum/extremum.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Extremum_test)

void print_point(TPointList points) {
  for (SPoint p : points) {
    BOOST_TEST_MESSAGE("(" << p.X << ", " << p.Y << "), ");
  }
  BOOST_TEST_MESSAGE("=======");
}

void check_points(TPointList points, TPointList test_points) {
  std::cout << " returned peaks: " << points.size() << "\n";
  for (Int32 i = 0; i < ssize(points); i++) {
    std::cout << " returned points: " << i << " " << points[i].X << " "
              << points[i].Y << " " << test_points[i].X << " "
              << test_points[i].Y << "\n";
  }

  BOOST_CHECK(points.size() == test_points.size());
  for (int i = 0; i < ssize(points); ++i) {
    BOOST_CHECK(
        (points[i].X == test_points[i].X && points[i].Y == test_points[i].Y));
  }
}

BOOST_AUTO_TEST_CASE(Extremum1) {
  CExtremum peaks1;

  TPointList maxPoint;
  TFloat64List x = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  TFloat64List y = {1.0, 5.0, 0.0, 3.0, 9.0, 0.1, 0.2, 0.5, 8.0, 1.0, 4.0};
  TFloat64List y_2 = {5.0, 0.0, 0.0, 3.0, 9.0, 0.1, 0.2, 0.5, 8.0, 1.0, 4.0};
  TFloat64List x_empty = {};

  peaks1.SetMaxPeakCount(2);
  peaks1.SetXRange(TFloat64Range(-10.0, 10.0));

  maxPoint = peaks1.Find(x, y, false);
  check_points(maxPoint, TPointList({{0.4, 9}, {0.8, 8}}));

  maxPoint = peaks1.Find(x, y_2, false);
  check_points(maxPoint, TPointList({{0.4, 9}, {0.8, 8}}));

  BOOST_CHECK_THROW(maxPoint = peaks1.Find(x_empty, y, true), AmzException);
  print_point(maxPoint);
}

BOOST_AUTO_TEST_CASE(Extremum_PlankDetection) {
  CExtremum peaks1;
  TPointList maxPoint;
  TFloat64List x = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  TFloat64List y_plank = {5.0, 0.0, 9.0, 9.0, 9.0, 0.1,
                          0.2, 0.5, 8.0, 1.0, 4.0};
  TFloat64List y_plank_chair = {5.0, 9.5, 9.0, 9.0, 9.0, 0.1,
                                0.2, 0.5, 8.0, 1.0, 4.0};
  TFloat64List y_plank_stair = {5.0, 0.0, 9.0, 9.0, 9.0, 9.5,
                                0.2, 0.5, 8.0, 1.0, 4.0};
  TFloat64List y_plank_low = {5.0, 0.0, 0.0, 0.0, 0.0, 9.5,
                              0.2, 0.5, 8.0, 1.0, 4.0};

  Float64 radius = 0.005;
  peaks1.SetMaxPeakCount(2);
  peaks1.SetXRange(TFloat64Range(-10.0, 10.0));
  bool isFirstPass = false;

  maxPoint = peaks1.Find(x, y_plank, isFirstPass);
  check_points(maxPoint, TPointList({{0.3, 9}, {0.8, 8}}));

  maxPoint = peaks1.Find(x, y_plank_chair, isFirstPass);
  check_points(maxPoint, TPointList({{0.1, 9.5}, {0.8, 8}}));

  maxPoint = peaks1.Find(x, y_plank_stair, isFirstPass);
  check_points(maxPoint, TPointList({{0.5, 9.5}, {0.8, 8}}));

  maxPoint = peaks1.Find(x, y_plank_low, isFirstPass);
  check_points(maxPoint, TPointList({{0.5, 9.5}, {0.8, 8}}));
}

BOOST_AUTO_TEST_CASE(Extremum_cut_isolated) {
  CExtremum peaks1;
  TPointList maxPoint;
  TFloat64List x = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  TFloat64List y = {1.0, 5.0, 0.0, 3.0, 1.0, 0.1, 6.0, 0.5, 8.0, 1.0, 4.0};

  Float64 radius = 0.005;
  peaks1.SetMaxPeakCount(5);
  peaks1.SetXRange(TFloat64Range(-10.0, 10.0));
  peaks1.SetMeritCut(3);

  peaks1.SortIndexes(y);

  peaks1.Cut_Threshold(x, y, 2);
  for (Int32 i = 0; i < ssize(x); i++) {
    maxPoint.push_back(SPoint(x[i], y[i]));
  }
  check_points(maxPoint, TPointList({{0.1, 5}, {0.6, 6}, {0.8, 8}}));
}

BOOST_AUTO_TEST_CASE(Extremum_cut_onePeak) {
  CExtremum peaks1;
  TPointList maxPoint;
  TFloat64List x = {0.0, 0.1, 0.2};
  TFloat64List y = {1.0, 6.0, 1.2};

  Float64 radius = 0.005;
  peaks1.SetMaxPeakCount(5);
  peaks1.SetXRange(TFloat64Range(-10.0, 10.0));
  peaks1.SetMeritCut(3);

  peaks1.SortIndexes(y);

  peaks1.Cut_Threshold(x, y, 2);
  for (Int32 i = 0; i < ssize(x); i++) {
    maxPoint.push_back(SPoint(x[i], y[i]));
  }
  check_points(maxPoint, TPointList({{0.1, 6}, {0.2, 1.2}}));
}
BOOST_AUTO_TEST_CASE(Extremum_FilterOutNeighboringPeaks) {

  TFloat64List x = {0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  TFloat64List y = {1.0, 5.0, 0.0, 3.0, 1.0, 0.1, 6.0, 0.5, 8.0, 1.0, 4.0};

  Float64 radius = 0.2;
  CExtremum peaks1(5, radius, -1, false, true, TFloat64Range(-10.0, 10.0));
  peaks1.SortIndexes(y);
  TPointList maxPoint = peaks1.FilterOutNeighboringPeaksAndTruncate(x, y, 2);
  check_points(maxPoint, TPointList({{0.8, 8}, {0.1, 5.0}, {0.4, 1.0}}));
}

BOOST_AUTO_TEST_CASE(Extremum_FilterOutNeighboringPeaks_consecutifClosePeaks) {

  TFloat64List x = {0.01, 0.02, 0.025, 0.03, 0.04, 0.05,
                    0.06, 0.07, 0.08,  0.09, 0.1};
  TFloat64List y = {1.0, 5.0, 0.0, 3.0, 1.0, 0.1, 6.0, 0.5, 8.0, 1.0, 4.0};

  Float64 radius = 0.05;
  CExtremum peaks1(5, radius, -1, false, true, TFloat64Range(-10.0, 10.0));
  peaks1.SortIndexes(y);
  TPointList maxPoint = peaks1.FilterOutNeighboringPeaksAndTruncate(x, y, 2);
  check_points(maxPoint, TPointList({{0.08, 8}, {0.02, 5}}));
}

BOOST_AUTO_TEST_CASE(Extremum_FilterOutNeighboringPeaks_consecutifClosePeaks2) {

  TFloat64List x = {0.01, 0.02, 0.025, 0.03, 0.04, 0.05,
                    0.06, 0.07, 0.08,  0.09, 0.1};
  TFloat64List y = {1.0, 3.0, 0.0, 4.0, 1.0, 0.1, 6.0, 0.5, 8.0, 1.0, 4.0};

  Float64 radius = 0.005;
  CExtremum peaks1(5, radius, -1, false, true, TFloat64Range(-10.0, 10.0));
  peaks1.SortIndexes(y);
  TPointList maxPoint = peaks1.FilterOutNeighboringPeaksAndTruncate(x, y, 2);
  check_points(
      maxPoint,
      TPointList({{0.08, 8},
                  {0.06, 6.0},
                  {0.03, 4.0},
                  {0.1, 4.0},
                  {0.02, 3}})); // here we lose a good peak (0.6, 6) since it is
                                // close to 0.8, within sliding windows
}

BOOST_AUTO_TEST_CASE(Extremum_FilterOutNeighboringPeaks_increasingPDF) {

  TFloat64List x = {0.01, 0.02, 0.025, 0.03, 0.04, 0.05,
                    0.06, 0.07, 0.08,  0.09, 0.1};
  TFloat64List y = {1.0, 3.0, 4.0, 5.0, 6.0, 7.1, 8.0, 9.5, 10.0, 11.0, 14.0};

  Float64 radius = 0.005;
  CExtremum peaks1(5, radius, -1, false, true, TFloat64Range(-10.0, 10.0));
  peaks1.SortIndexes(y);
  TPointList maxPoint = peaks1.FilterOutNeighboringPeaksAndTruncate(x, y, 2);

  check_points(
      maxPoint,
      TPointList({{0.1, 14.0},
                  {0.09, 11.0},
                  {0.08, 10.0},
                  {0.07, 9.5},
                  {0.06, 8.0}})); // here we lose a good peak (0.6, 6) since it
                                  // is close to 0.8, within sliding windows

  // bigger radius
  radius = 0.1;
  maxPoint = peaks1.FilterOutNeighboringPeaksAndTruncate(x, y, 2);
  check_points(
      maxPoint,
      TPointList({{0.1, 14.0},
                  {0.09, 11.0},
                  {0.08, 10.0},
                  {0.07, 9.5},
                  {0.06, 8.0}})); // here we lose a good peak (0.6, 6) since it
                                  // is close to 0.8, within sliding windows
}
// testing case where nb of peaks = keepmin
BOOST_AUTO_TEST_CASE(Extremum_FilterOutNeighboringPeaks2) {

  TFloat64List x = {0.01, 0.1};
  TFloat64List y = {1.0, 5.0};

  Float64 radius = 0.2;
  CExtremum peaks1(5, radius, -1, false, true, TFloat64Range(-10.0, 10.0));
  peaks1.SortIndexes(y);
  // testing only the sliding window algo
  Int32 keepmin = 2;
  TPointList maxPoint =
      peaks1.FilterOutNeighboringPeaksAndTruncate(x, y, keepmin);
  check_points(
      maxPoint,
      TPointList({{0.1, 5.0},
                  {0.01, 1.0}})); // here we lose a good peak (0.6, 6) since it
                                  // is close to 0.8, within sliding windows
}

// testing case where nb of peaks = keepmin
BOOST_AUTO_TEST_CASE(Extremum_FilterOutNeighboringPeaks_negativeX) {

  TFloat64List x = {-0.000750, -0.000650};
  TFloat64List y = {3.0, 4.0};

  Float64 peakseparationDist = 0.005 * 2;
  CExtremum peaks1(1, peakseparationDist, -1, false, true,
                   TFloat64Range(-10.0, 10.0));
  peaks1.SortIndexes(y);
  // testing only the sliding window algo
  Int32 keepmin = 1;
  TPointList maxPoint =
      peaks1.FilterOutNeighboringPeaksAndTruncate(x, y, keepmin);
  check_points(
      maxPoint,
      TPointList(
          {{-0.000650, 4.0}})); // here we lose a good peak (0.6, 6) since it is
                                // close to 0.8, within sliding windows
}

BOOST_AUTO_TEST_CASE(Extremum_Truncate) {

  TFloat64List x = {0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  TFloat64List y = {1.0, 5.0, 0.0, 3.0, 1.0, 0.1, 6.0, 0.5, 8.0, 1.0, 4.0};

  Int32 maxCount = 2;
  CExtremum peaks1(maxCount, 0.005 * 2, -1, false, true,
                   TFloat64Range(-10.0, 10.0));
  peaks1.SortIndexes(y);
  TPointList maxPoint = peaks1.Truncate(x, y);

  check_points(maxPoint, TPointList({{0.8, 8}, {0.6, 6}}));
}
BOOST_AUTO_TEST_SUITE_END()
