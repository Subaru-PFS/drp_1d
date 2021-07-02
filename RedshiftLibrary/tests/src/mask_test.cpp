#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/common/datatypes.h"

#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Mask)

BOOST_AUTO_TEST_CASE(Mask1)
{
  CMask mask;
  const NSEpic::Mask *data;

  mask.SetSize(4);

  BOOST_CHECK( mask.GetMasksCount() == 4 );

  mask[0] = 0;
  mask[1] = 1;
  mask[2] = 2;
  mask[3] = 3;

  for (UInt8 n = 0; n < mask.GetMasksCount(); n++) {
    BOOST_CHECK( mask[n] == n );
  }

  data = mask.GetMasks();
  BOOST_CHECK( data[1] == 1 && data[3] == 3 );
}


BOOST_AUTO_TEST_CASE(Mask2)
{
  CMask mask1(2);
  CMask mask2(2);

  BOOST_CHECK( mask1.GetMasksCount() == 2 );

  mask1[0] = 0xf0;
  mask1[1] = 0x80;
  mask2[0] = 0x0f;
  mask2[1] = 0xf0; // mask1=f080 mask2=0ff0

  mask2 &= mask1; // mask1=f080 mask2=0080

  BOOST_CHECK( mask2[0] == 0 && mask2[1] == 0x80 );

  BOOST_CHECK( mask1.IntersectWith(mask2) == true );  // mask1=0080 mask2=0080
  BOOST_CHECK( mask1[0] == 0 && mask1[1] == 0x80 );

  mask1[1] = 0x7f;                                   // mask1=007f mask2=0080
  BOOST_CHECK( mask1.IntersectWith(mask2) == true ); // mask1=0000 mask2=0080
  BOOST_CHECK( mask1[0] == 0 && mask1[1] == 0 );
  BOOST_CHECK( mask2[0] == 0 && mask2[1] == 0x80 );

  mask1[0] = 0x0f;
  mask1[1] = 0xf0;
  mask2[0] = 0xf0;
  mask2[1] = 0x0f; // mask1=0ff0 mask2=f00f
  BOOST_CHECK_CLOSE( mask1.CompouteOverlapRate(mask2), 1.0, 1e-6 );

  mask1[0] = 0;
  mask1[1] = 0;    // mask1=0000 mask2=f00f
  BOOST_CHECK_CLOSE( mask1.CompouteOverlapRate(mask2), 0, 1e-6 );
  // TODO : should better raise an exception

  mask1.SetSize(3);
  mask1[2] = 0x01; // mask1=000001 mask2=f00f
  BOOST_CHECK_CLOSE( mask1.CompouteOverlapRate(mask2), -1.0, 1e-6  );
  // TODO : should better raise an exception

  mask1 &= mask2;
  BOOST_CHECK( mask1.GetMasksCount() == 3);
  BOOST_CHECK( mask1[0] == 0 && mask1[1] == 0 && mask1[2] == 1);


  BOOST_CHECK( mask1.IntersectWith(mask2) == false );

  BOOST_CHECK_CLOSE( mask1.IntersectAndComputeOverlapRate(mask2), -1.0, 1e-6);

  mask1.SetSize(2);
  mask1[0] = 0x0f;
  mask1[1] = 0xf0; // mask1=0ff0 mask2=f00f
  BOOST_CHECK_CLOSE( mask1.IntersectAndComputeOverlapRate(mask2), 0.0, 1e-6);

  mask1[0] = 0x12;
  mask1[1] = 0xef; // mask1=12ef mask2=f00f
  BOOST_CHECK_CLOSE( mask1.IntersectAndComputeOverlapRate(mask2), 0.1206, 0.19);

  mask1[0] = 0;
  mask1[1] = 0;    // mask1=0000 mask2=f00f
  BOOST_CHECK_CLOSE( mask1.IntersectAndComputeOverlapRate(mask2), 0, 1e-6 );
  // TODO : should better raise an exception

  BOOST_CHECK( mask2.GetUnMaskedSampleCount() == 0xff );

  // TODO: GetMaskedSampleCount is dead code
  mask1[0] = 0;
  mask1[0] = 1;
  BOOST_CHECK( mask1.GetMaskedSampleCount() == 1 );
}

BOOST_AUTO_TEST_SUITE_END()
