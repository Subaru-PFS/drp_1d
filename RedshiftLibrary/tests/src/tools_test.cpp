#include <RedshiftLibrary/spectrum/tools.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/common/mask.h>

#include <boost/test/unit_test.hpp>

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(test_tools)

BOOST_AUTO_TEST_CASE(Interpolate)
{
  CSpectrumTools tools;
  Float64 xOrg[] = {1,2,3,4};
  Float64 yOrg[] = {10,20,30,40};
  Float64 xInt[] = {1.5,2.5,3.5,4.5};
  Float64 xIntNoInt[] = {1,2.5,3,4};

  CSpectrumAxis axisXorg( xOrg, 4);
  CSpectrumAxis axisYorg( yOrg, 4);
  CSpectrumAxis axisXint( xInt, 4);
  CSpectrumAxis axisXintNoInt( xIntNoInt, 4);
  CSpectrumAxis axisYint(4);
  Int32 offset;
  Int32 nOrg = 4;
  CMask mask(4);

  // first test
  offset = 0;
  nOrg = 4;
  tools.Interpolate(axisXorg, axisYorg, offset, nOrg, axisXint, axisYint, mask);
  BOOST_TEST_MESSAGE("out: " << axisYint[0] << " " << axisYint[1] << " " << axisYint[2] << " " << axisYint[3]);
  BOOST_TEST_MESSAGE("mask: " << (mask[0] ? "1" : "0")  << " " << (mask[1] ? "1" : "0") << " " << (mask[2] ? "1" : "0") << " " << (mask[3] ? "1" : "0"));
  BOOST_CHECK_CLOSE(axisYint[0], 15., 1e-18);
  BOOST_CHECK_CLOSE(axisYint[1], 25., 1e-18);
  BOOST_CHECK_CLOSE(axisYint[2], 35., 1e-18);
  BOOST_CHECK_CLOSE(axisYint[3], 0., 1e-18);
  BOOST_CHECK(mask[0] == 1 && mask[1] == 1 && mask[2] == 1 && mask[3] == 0);

  // second test
  axisYint[0] = 0;
  axisYint[1] = 0;
  axisYint[2] = 0;
  axisYint[3] = 0;
  offset = 1;
  nOrg = 2;
  tools.Interpolate(axisXorg, axisYorg, offset, nOrg, axisXint, axisYint, mask);
  BOOST_TEST_MESSAGE("out: " << axisYint[0] << " " << axisYint[1] << " " << axisYint[2] << " " << axisYint[3]);
  BOOST_TEST_MESSAGE("mask: " << (mask[0] ? "1" : "0")  << " " << (mask[1] ? "1" : "0") << " " << (mask[2] ? "1" : "0") << " " << (mask[3] ? "1" : "0"));
  BOOST_CHECK_CLOSE(axisYint[0], 0., 1e-18);
  BOOST_CHECK_CLOSE(axisYint[1], 25., 1e-18);
  BOOST_CHECK_CLOSE(axisYint[2], 0., 1e-18);
  BOOST_CHECK_CLOSE(axisYint[3], 0., 1e-18);
  BOOST_CHECK(mask[0] == 0 && mask[1] == 1 && mask[2] == 0 && mask[3] == 0);

  // third test
  axisYint[0] = 0;
  axisYint[1] = 0;
  axisYint[2] = 0;
  axisYint[3] = 0;
  offset = 0;
  nOrg = 4;
  tools.Interpolate(axisXorg, axisYorg, offset, nOrg, axisXintNoInt, axisYint, mask);
  BOOST_TEST_MESSAGE("out: " << axisYint[0] << " " << axisYint[1] << " " << axisYint[2] << " " << axisYint[3]);
  BOOST_TEST_MESSAGE("mask: " << (mask[0] ? "1" : "0")  << " " << (mask[1] ? "1" : "0") << " " << (mask[2] ? "1" : "0") << " " << (mask[3] ? "1" : "0"));
  BOOST_CHECK_CLOSE(axisYint[0], 10., 1e-18);
  BOOST_CHECK_CLOSE(axisYint[1], 25., 1e-18);
  BOOST_CHECK_CLOSE(axisYint[2], 30., 1e-18);
  BOOST_CHECK_CLOSE(axisYint[3], 40., 1e-18);
  BOOST_CHECK(mask[0] == 1 && mask[1] == 1 && mask[2] == 1 && mask[3] == 1);

}

BOOST_AUTO_TEST_SUITE_END()
