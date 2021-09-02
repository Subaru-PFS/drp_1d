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
