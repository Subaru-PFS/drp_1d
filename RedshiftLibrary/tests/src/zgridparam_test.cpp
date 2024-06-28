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
#include <boost/test/unit_test.hpp>

#include "RedshiftLibrary/common/zgridparam.h"

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(CZGridParam_test)

BOOST_AUTO_TEST_CASE(CZGridListParams_test) {
  // test default ctor
  CZGridParam p;

  // test ctor with empty args
  CZGridListParams lp({});

  BOOST_CHECK(lp.size() == 0);

  // test getters
  CZGridParam p0(TFloat64Range(0, 10), 3, NAN);
  CZGridParam p1(TFloat64Range(2, 5), 1, 3);
  lp = CZGridListParams({p0, p1});

  BOOST_CHECK(lp.size() == 2);

  p = lp[0];
  BOOST_CHECK(p.zmin == p0.zmin);
  BOOST_CHECK(p.zmax == p0.zmax);
  BOOST_CHECK(p.zstep == p0.zstep);
  BOOST_CHECK(std::isnan(p.zcenter));

  p = lp[1];
  BOOST_CHECK(p.zmin == p1.zmin);
  BOOST_CHECK(p.zmax == p1.zmax);
  BOOST_CHECK(p.zstep == p1.zstep);
  BOOST_CHECK(p.zcenter == p1.zcenter);

  // test iterators
  for (auto it = lp.cbegin(); it != lp.cend(); ++it) {
    if (it == lp.cbegin()) {
      BOOST_CHECK(it->zmin == p0.zmin);
      BOOST_CHECK(it->zmax == p0.zmax);
      BOOST_CHECK(it->zstep == p0.zstep);
      BOOST_CHECK(std::isnan(it->zcenter));
    } else {
      BOOST_CHECK(it->zmin == p1.zmin);
      BOOST_CHECK(it->zmax == p1.zmax);
      BOOST_CHECK(it->zstep == p1.zstep);
      BOOST_CHECK(it->zcenter == p1.zcenter);
    }
  }

  for (auto it = lp.begin(); it != lp.end(); ++it) {
    if (it == lp.cbegin()) {
      BOOST_CHECK(it->zmin == p0.zmin);
      BOOST_CHECK(it->zmax == p0.zmax);
      BOOST_CHECK(it->zstep == p0.zstep);
      BOOST_CHECK(std::isnan(it->zcenter));
      BOOST_CHECK_NO_THROW(it->zmin = 0.);
      BOOST_CHECK(it->zmin == 0.);
    } else {
      BOOST_CHECK(it->zmin == p1.zmin);
      BOOST_CHECK(it->zmax == p1.zmax);
      BOOST_CHECK(it->zstep == p1.zstep);
      BOOST_CHECK(it->zcenter == p1.zcenter);
      BOOST_CHECK_NO_THROW(it->zmin = 0.);
      BOOST_CHECK(it->zmin == 0.);
    }
  }
}

BOOST_AUTO_TEST_CASE(regular_grid_test) {

  CZGridListParams zparams({});
  TFloat64List grid = zparams.getZGrid(false);
  BOOST_CHECK(grid.size() == 0);

  zparams = CZGridListParams(
      {CZGridParam(TFloat64Range(0, 10), 3, NAN)}); // 0, 3, 6, 9
  grid = CZGridListParams(zparams).getZGrid(false);

  TFloat64List resultingVect{0, 3, 6, 9};
  BOOST_CHECK(grid.size() == resultingVect.size());
  BOOST_CHECK(grid == resultingVect);
}

BOOST_AUTO_TEST_CASE(mixedGrid_test) {

  TZGridListParams zparams = {CZGridParam(TFloat64Range(0, 10), 3, NAN),
                              CZGridParam(TFloat64Range(2, 5), 1, 3)};
  const TFloat64List mixedGrid = CZGridListParams(zparams).getZGrid(false);

  TFloat64List resultingVect{0, 2, 3, 4, 5, 6, 9};
  BOOST_CHECK(mixedGrid.size() == resultingVect.size());
  BOOST_CHECK(mixedGrid == resultingVect);
}

BOOST_AUTO_TEST_CASE(mixedGrid_zcenterNotincluded_test) {

  TZGridListParams zparams = {
      CZGridParam(TFloat64Range(0, 10), 3, NAN), // 0, 3, 6, 9
      CZGridParam(TFloat64Range(2, 5), 1, 4)};
  const TFloat64List mixedGrid = CZGridListParams(zparams).getZGrid(false);

  TFloat64List resultingVect{0, 2, 3, 4, 5, 6, 9};
  BOOST_CHECK(mixedGrid.size() == resultingVect.size());
  BOOST_CHECK(mixedGrid == resultingVect);
}

BOOST_AUTO_TEST_CASE(mixedGrid_zcenterNotincluded2_test) {

  TZGridListParams zparams = {
      CZGridParam(TFloat64Range(1, 10), 3, NAN), // 1, 4, 7, 10
      CZGridParam(TFloat64Range(2, 5), 0.8, 4)};
  const TFloat64List mixedGrid = CZGridListParams(zparams).getZGrid(false);

  TFloat64List resultingVect{1, 2.4, 3.2, 4, 4.8, 7, 10};
  BOOST_CHECK(mixedGrid.size() == resultingVect.size());
  BOOST_CHECK(mixedGrid == resultingVect);
}

BOOST_AUTO_TEST_CASE(mixedGrid_check_error_on_single_element_test) {
  // If no intersection in the given grids, error
  TZGridListParams zparams = {
      CZGridParam(TFloat64Range(1, 10), 3, NAN), // 1, 4, 7, 10
      CZGridParam(TFloat64Range(4., 5.), 0.5, 4.5)};
  BOOST_CHECK_THROW(CZGridListParams(zparams).getZGrid(true), AmzException);
}

BOOST_AUTO_TEST_CASE(mixedGrid_truncated_subgrids_test) {
  CZGridParam p0(TFloat64Range(1.01, 10.01), 3, NAN);
  CZGridParam p1(TFloat64Range(0.5, 3.), 0.5, 2.5);
  TZGridListParams zparams = {p0, p1};

  TFloat64List mixedGrid = CZGridListParams(zparams).getZGrid(false);
  BOOST_CHECK(mixedGrid.front() == p0.zmin);

  CZGridParam p1bis(TFloat64Range(6.5, 10.5), 0.5, 8.5);
  zparams = {p0, p1bis};
  mixedGrid = CZGridListParams(zparams).getZGrid(false);
  BOOST_CHECK(mixedGrid.back() == p0.zmax);
}

BOOST_AUTO_TEST_CASE(grid_idempotence_test) {

  // with no center
  CZGridParam zparam(TFloat64Range(0, 10), 3, NAN); // 0, 3, 6, 9
  TFloat64List grid = zparam.getZGrid(false);

  CZGridParam zparam2(TFloat64Range(grid), 3, NAN);
  TFloat64List grid2 = zparam2.getZGrid(false);

  BOOST_CHECK(grid == grid2);

  // with center
  zparam = CZGridParam(TFloat64Range(0, 11), 3, 4);
  grid = zparam.getZGrid(false); // 1 4 7 10

  zparam2 = CZGridParam(TFloat64Range(grid), 3, 4);
  grid2 = zparam2.getZGrid(false);

  BOOST_CHECK(grid == grid2);

  // in log without center
  zparam = CZGridParam(TFloat64Range(0, 11), log(1 + 0.1), NAN);
  grid = zparam.getZGrid(true);
  BOOST_CHECK(TFloat64Range(grid) != TFloat64Range(0, 11));
  zparam2 = CZGridParam(TFloat64Range(grid), log(1 + 0.1), NAN);
  grid2 = zparam2.getZGrid(true);
  BOOST_CHECK(grid == grid2);

  // in log with center
  zparam = CZGridParam(TFloat64Range(0, 11), log(1 + 0.1), 5);
  grid = zparam.getZGrid(true);
  BOOST_CHECK(TFloat64Range(grid) != TFloat64Range(0, 5));
  zparam2 = CZGridParam(TFloat64Range(grid), log(1 + 0.1), 5);
  grid2 = zparam2.getZGrid(true);
  BOOST_CHECK(grid == grid2);
}

BOOST_AUTO_TEST_SUITE_END()
