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
#include <iostream>

#include <boost/test/unit_test.hpp>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/mask.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Mask)

BOOST_AUTO_TEST_CASE(Mask1) {
  CMask mask;
  const NSEpic::Mask *data;

  mask.SetSize(4);

  BOOST_CHECK(mask.GetMasksCount() == 4);

  mask[0] = 0;
  mask[1] = 1;
  mask[2] = 2;
  mask[3] = 3;

  for (Int32 n = 0; n < (Int32)mask.GetMasksCount(); n++) {
    BOOST_CHECK(mask[n] == n);
  }

  data = mask.GetMasks();
  BOOST_CHECK(data[1] == 1 && data[3] == 3);
}

BOOST_AUTO_TEST_CASE(Mask2) {
  CMask mask1(2);
  CMask mask2(2);

  BOOST_CHECK(mask1.GetMasksCount() == 2);

  mask1[0] = 0xf0;
  mask1[1] = 0x80;
  mask2[0] = 0x0f;
  mask2[1] = 0xf0; // mask1=f080 mask2=0ff0

  mask2 &= mask1; // mask1=f080 mask2=0080

  BOOST_CHECK(mask2[0] == 0 && mask2[1] == 0x80);

  BOOST_CHECK(mask1.IntersectWith(mask2) == true); // mask1=0080 mask2=0080
  BOOST_CHECK(mask1[0] == 0 && mask1[1] == 0x80);

  mask1[1] = 0x7f;                                 // mask1=007f mask2=0080
  BOOST_CHECK(mask1.IntersectWith(mask2) == true); // mask1=0000 mask2=0080
  BOOST_CHECK(mask1[0] == 0 && mask1[1] == 0);
  BOOST_CHECK(mask2[0] == 0 && mask2[1] == 0x80);

  mask1[0] = 0x0f;
  mask1[1] = 0xf0;
  mask2[0] = 0xf0;
  mask2[1] = 0x0f; // mask1=0ff0 mask2=f00f
  BOOST_CHECK_CLOSE(mask1.CompouteOverlapFraction(mask2), 1.0, 1e-6);

  mask1[0] = 0;
  mask1[1] = 0; // mask1=0000 mask2=f00f
  BOOST_CHECK_CLOSE(mask1.CompouteOverlapFraction(mask2), 0, 1e-6);
  // TODO : should better raise an exception

  mask1.SetSize(3);
  mask1[2] = 0x01; // mask1=000001 mask2=f00f
  BOOST_CHECK_CLOSE(mask1.CompouteOverlapFraction(mask2), -1.0, 1e-6);
  // TODO : should better raise an exception

  mask1 &= mask2;
  BOOST_CHECK(mask1.GetMasksCount() == 3);
  BOOST_CHECK(mask1[0] == 0 && mask1[1] == 0 && mask1[2] == 1);

  BOOST_CHECK(mask1.IntersectWith(mask2) == false);

  BOOST_CHECK_CLOSE(mask1.IntersectAndComputeOverlapFraction(mask2), -1.0,
                    1e-6);

  mask1.SetSize(2);
  mask1[0] = 0x0f;
  mask1[1] = 0xf0; // mask1=0ff0 mask2=f00f
  BOOST_CHECK_CLOSE(mask1.IntersectAndComputeOverlapFraction(mask2), 0.0, 1e-6);

  mask1[0] = 0x12;
  mask1[1] = 0xef; // mask1=12ef mask2=f00f
  BOOST_CHECK_CLOSE(mask1.IntersectAndComputeOverlapFraction(mask2), 0.1206,
                    0.19);

  mask1[0] = 0;
  mask1[1] = 0; // mask1=0000 mask2=f00f
  BOOST_CHECK_CLOSE(mask1.IntersectAndComputeOverlapFraction(mask2), 0, 1e-6);
  // TODO : should better raise an exception

  BOOST_CHECK(mask2.GetUnMaskedSampleCount() == 0xff);

  // TODO: GetMaskedSampleCount is dead code
  mask1[0] = 0;
  mask1[0] = 1;
  BOOST_CHECK(mask1.GetMaskedSampleCount() == 1);
}

BOOST_AUTO_TEST_SUITE_END()
