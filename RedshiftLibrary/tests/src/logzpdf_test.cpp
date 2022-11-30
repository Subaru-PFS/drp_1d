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
#include "RedshiftLibrary/operator/logZPdfResult.h"
#include <boost/test/unit_test.hpp>

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(logzpdfTest)

BOOST_AUTO_TEST_CASE(regular_grid) {
  TZGridListParams zparams(1);
  zparams[0] = CZGridParam(TFloat64Range(0, 10), 3, NAN); // 0, 3, 6, 9
  //  zparams[1] = CZGridParam(TFloat64Range(2, 5), 1, 4);
  const TFloat64List grid = CZGridListParams(zparams).getZGrid(false);

  TFloat64List resultingVect{0, 3, 6, 9};
  BOOST_CHECK(grid.size() == resultingVect.size());
  BOOST_CHECK(grid == resultingVect);
}
BOOST_AUTO_TEST_CASE(mixedGrid_test) {

  TZGridListParams zparams(2);
  zparams[0] = CZGridParam(TFloat64Range(0, 10), 3, NAN);
  zparams[1] = CZGridParam(TFloat64Range(2, 5), 1, 3);
  const TFloat64List mixedGrid = CZGridListParams(zparams).getZGrid(false);

  TFloat64List resultingVect{0, 2, 3, 4, 5, 6, 9};
  BOOST_CHECK(mixedGrid.size() == resultingVect.size());
  BOOST_CHECK(mixedGrid == resultingVect);
}

BOOST_AUTO_TEST_CASE(mixedGrid_test_zcenterNotincluded) {

  TZGridListParams zparams(2);
  zparams[0] = CZGridParam(TFloat64Range(0, 10), 3, NAN); // 0, 3, 6, 9
  zparams[1] = CZGridParam(TFloat64Range(2, 5), 1, 4);
  const TFloat64List mixedGrid = CZGridListParams(zparams).getZGrid(false);

  TFloat64List resultingVect{0, 2, 3, 4, 5, 6, 9};
  BOOST_CHECK(mixedGrid.size() == resultingVect.size());
  BOOST_CHECK(mixedGrid == resultingVect);
}

BOOST_AUTO_TEST_CASE(mixedGrid_test_zcenterNotincluded2) {

  TZGridListParams zparams(2);
  zparams[0] = CZGridParam(TFloat64Range(1, 10), 3, NAN); // 1, 4, 7, 10
  zparams[1] = CZGridParam(TFloat64Range(2, 5), 0.8, 4);
  const TFloat64List mixedGrid = CZGridListParams(zparams).getZGrid(false);

  TFloat64List resultingVect{1, 2.4, 3.2, 4, 4.8, 7, 10};
  BOOST_CHECK(mixedGrid.size() == resultingVect.size());
  BOOST_CHECK(mixedGrid == resultingVect);
}

BOOST_AUTO_TEST_CASE(check_error_on_single_element) {
  TZGridListParams zparams(3);
  zparams[0] = CZGridParam(TFloat64Range(1, 10), 3, NAN); // 1, 4, 7, 10
  zparams[1] = CZGridParam(TFloat64Range(2., 3.), 0.5, 2.5);
  zparams[2] = CZGridParam(TFloat64Range(4., 5.), 0.5, 4.5);

  BOOST_CHECK_THROW(
      CZGridListParams(zparams).getZGrid(true),
      GlobalException); // log grid with single element in main grid
}
BOOST_AUTO_TEST_SUITE_END()
