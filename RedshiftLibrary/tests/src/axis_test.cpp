#include <boost/test/unit_test.hpp>
#include "RedshiftLibrary/spectrum/axis.h"

using namespace NSEpic;

BOOST_AUTO_TEST_CASE(SampleCount)
{
  // No element
  CSpectrumAxis n0Axis = CSpectrumAxis();
  BOOST_CHECK(n0Axis.GetSamplesCount() == 0);
  // One element
  CSpectrumAxis n10Axis = CSpectrumAxis(1);
  BOOST_CHECK(n10Axis.GetSamplesCount() == 1);
  Float64 n1Array[] = {0.5};
  CSpectrumAxis n11Axis = CSpectrumAxis(n1Array, 1);
  BOOST_CHECK(n11Axis.GetSamplesCount() == 1);
  // Two elements
  CSpectrumAxis n20Axis = CSpectrumAxis(2);
  BOOST_CHECK(n20Axis.GetSamplesCount() == 2);
  Float64 n2Array[] = {0.,1.};
  CSpectrumAxis n21Axis = CSpectrumAxis(n2Array, 2);
  BOOST_CHECK(n21Axis.GetSamplesCount() == 2);
}

BOOST_AUTO_TEST_CASE(GetSample)
{
  // No element
  CSpectrumAxis n0Axis = CSpectrumAxis();
  BOOST_CHECK(n0Axis.GetSamples() == NULL);
  // One element
  Float64 n1Array[] = {0.5};
  CSpectrumAxis n1Axis = CSpectrumAxis(n1Array, 1);
  BOOST_CHECK(n1Axis.GetSamples()[0] == n1Array[0]);
  // Two elements
  Float64 n2Array[] = {0.,1.};
  CSpectrumAxis n2Axis = CSpectrumAxis(n2Array, 2);
  BOOST_CHECK(n2Axis.GetSamples()[0] == n2Array[0]);
  BOOST_CHECK(n2Axis.GetSamples()[1] == n2Array[1]);
}

BOOST_AUTO_TEST_CASE(Operator)
{
  CSpectrumAxis n0Axis = CSpectrumAxis();
  BOOST_CHECK(n0Axis.GetSamples() == NULL);
  // One element
  Float64 n1Array[] = {0.5};
  CSpectrumAxis n10Axis = CSpectrumAxis(n1Array, 1);
  BOOST_CHECK(n10Axis[0] == n1Array[0]);
  const CSpectrumAxis n11Axis = CSpectrumAxis(n1Array, 1);
  BOOST_CHECK(n11Axis[0] == n1Array[0]);
  // Two elements
  Float64 n2Array[] = {0.,1.};
  CSpectrumAxis n2Axis = CSpectrumAxis(n2Array, 2);
  BOOST_CHECK(n2Axis[0] == n2Array[0]);
  BOOST_CHECK(n2Axis[1] == n2Array[1]);
}
