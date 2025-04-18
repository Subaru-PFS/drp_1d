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
#include <fstream>
#include <iostream>
#include <string>

#include <boost/test/execution_monitor.hpp>
#include <boost/test/unit_test.hpp>

#include "RedshiftLibrary/common/vectorOperations.h"
#include "RedshiftLibrary/operator/linemodelresult.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(vectorOperations)

BOOST_AUTO_TEST_CASE(insertWithDuplicates_float) {

  TFloat64List vect{10, 20, 30};
  Int32 idx = 1;
  Int32 ndup = 1;
  insertWithDuplicates<Float64>(vect, idx, 3, 0., ndup);

  TFloat64List resultingVect{10, 0, 0, 0, 30};
  BOOST_CHECK(vect.size() == resultingVect.size());
  BOOST_CHECK(vect == resultingVect);
}

BOOST_AUTO_TEST_CASE(insertWithDuplicates_float_begin) {

  TFloat64List vect{10, 20, 30};
  Int32 idx = 0;
  Int32 ndup = 1;
  insertWithDuplicates<Float64>(vect, idx, 3, 0., ndup);

  TFloat64List resultingVect{0., 0., 0., 20, 30};

  BOOST_CHECK(vect.size() == resultingVect.size());
  BOOST_CHECK(vect == resultingVect);
}

BOOST_AUTO_TEST_CASE(insertWithDuplicates_float_end) {

  TFloat64List vect{10};
  Int32 idx = 0;
  Int32 ndup = 1;
  insertWithDuplicates<Float64>(vect, idx, 3, 0., ndup);

  TFloat64List resultingVect{0., 0., 0.};

  BOOST_CHECK(vect.size() == resultingVect.size());
  BOOST_CHECK(vect == resultingVect);
}

BOOST_AUTO_TEST_CASE(insertWithDuplicates_mid) {

  TFloat64List vect{10, 20, 30};
  const TFloat64List extendedVect{18, 20, 22};
  Int32 idx = 1;
  Int32 ndup = 1;
  insertWithDuplicates(vect, idx, extendedVect, ndup);

  TFloat64List resultingVect{10, 18, 20, 22, 30};

  BOOST_CHECK(vect.size() == resultingVect.size());
  BOOST_CHECK(vect == resultingVect);
}

BOOST_AUTO_TEST_CASE(insertWithDuplicates_begin) {

  TFloat64List vect{10, 20, 30};
  const TFloat64List extendedVect{8, 10, 12};
  Int32 idx = 0;
  Int32 ndup = 1;
  insertWithDuplicates(vect, idx, extendedVect, ndup);

  TFloat64List resultingVect{8, 10, 12, 20, 30};

  BOOST_CHECK(vect.size() == resultingVect.size());
  BOOST_CHECK(vect == resultingVect);
}

BOOST_AUTO_TEST_CASE(insertWithDuplicates_end) {

  TFloat64List vect{10, 20, 30};
  const TFloat64List extendedVect{28, 30, 32};
  Int32 idx = 2;
  Int32 ndup = 1;
  insertWithDuplicates(vect, idx, extendedVect, ndup);

  TFloat64List resultingVect{10, 20, 28, 30, 32};

  BOOST_CHECK(vect.size() == resultingVect.size());
  BOOST_CHECK(vect == resultingVect);
}
BOOST_AUTO_TEST_CASE(insertWithDuplicates_commononborder1) {
  TFloat64List vect{10, 15, 20};
  const TFloat64List extendedVect{10, 12, 15, 16};
  Int32 idx = 0;
  Int32 ndup = 2;
  insertWithDuplicates(vect, idx, extendedVect, ndup);

  TFloat64List resultingVect{10, 12, 15, 16, 20};
  BOOST_CHECK(vect.size() == resultingVect.size());
  BOOST_CHECK(vect == resultingVect);
}

BOOST_AUTO_TEST_CASE(insertWithDuplicates_commononborder) {
  TFloat64List vect{10, 15, 20, 25, 30};
  const TFloat64List extendedVect{15, 16, 18, 20, 22, 25, 26, 30};
  Int32 idx = 1;
  Int32 ndup = 4;
  insertWithDuplicates(vect, idx, extendedVect, ndup);

  TFloat64List resultingVect{10, 15, 16, 18, 20, 22, 25, 26, 30};
  BOOST_CHECK(vect.size() == resultingVect.size());
  BOOST_CHECK(vect == resultingVect);
}

BOOST_AUTO_TEST_CASE(linearInterp_float) {
  auto r = linearInterp<Float64>(1, 2, 0);
  BOOST_CHECK_CLOSE(r, 1, 1e-7);

  r = linearInterp<Float64>(1, 2, 1);
  BOOST_CHECK_CLOSE(r, 2, 1e-7);

  r = linearInterp<Float64>(1, 2, 0.5);
  BOOST_CHECK_CLOSE(r, 1.5, 1e-7);
}

BOOST_AUTO_TEST_CASE(linearInterp_int) {
  auto r = linearInterp<Int32>(1, 2, 0);
  BOOST_CHECK(r == 1);

  r = linearInterp<Int32>(1, 2, 1);
  BOOST_CHECK(r == 2);

  r = linearInterp<Int32>(1, 2, 0.5001);
  BOOST_CHECK(r == 2);

  r = linearInterp<Int32>(1, 2, 0.4999);
  BOOST_CHECK(r == 1);
}

BOOST_AUTO_TEST_CASE(createLinearInterpVector_float) {
  auto r = createLinearInterpVector<Float64>(1, 2, 5);
  TFloat64List ref = {1, 1.25, 1.5, 1.75, 2};
  BOOST_CHECK(r == ref);

  r = createLinearInterpVector<Float64>(1, 2, 1);
  ref = {1};
  BOOST_CHECK(r == ref);

  r = createLinearInterpVector<Float64>(1, 2, 0);
  ref = {};
  BOOST_CHECK(r == ref);
}

BOOST_AUTO_TEST_CASE(interpolateBetweenDuplicates_float) {
  TFloat64List source = {-1, 0, 1, 2, 3};
  Int32 insertionIdx = 1;
  TInt32List overwrittenSourceIndices = {0, 4, 8};
  Int32 nInterp = 9;
  Int32 largeStepFactor = 4;
  auto r = interpolateBetweenDuplicates<Float64>(
      source, insertionIdx, overwrittenSourceIndices, nInterp, largeStepFactor);
  TFloat64List ref = {0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2};
  BOOST_CHECK(r == ref);

  // incomplete first and last segments
  insertionIdx = 2;
  overwrittenSourceIndices = {2};
  nInterp = 5;
  r = interpolateBetweenDuplicates<Float64>(
      source, insertionIdx, overwrittenSourceIndices, nInterp, largeStepFactor);
  ref = {0.5, 0.75, 1, 1.25, 1.5};
  BOOST_CHECK(r == ref);

  insertionIdx = 0;
  BOOST_CHECK_THROW(interpolateBetweenDuplicates<Float64>(
                        source, insertionIdx, overwrittenSourceIndices, nInterp,
                        largeStepFactor),
                    AmzException);

  insertionIdx = 5;
  BOOST_CHECK_THROW(interpolateBetweenDuplicates<Float64>(
                        source, insertionIdx, overwrittenSourceIndices, nInterp,
                        largeStepFactor),
                    AmzException);
}
BOOST_AUTO_TEST_SUITE_END()