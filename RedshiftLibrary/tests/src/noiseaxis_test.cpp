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
#include <cmath>

#include "RedshiftLibrary/spectrum/noiseaxis.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(noiseaxis_test)

BOOST_AUTO_TEST_CASE(constructor_test) {

  CSpectrumNoiseAxis noiseAxis;
  BOOST_CHECK(noiseAxis.GetSamplesCount() == 0);

  CSpectrumNoiseAxis noiseAxis2(1);
  BOOST_CHECK(noiseAxis2.GetSamplesCount() == 1);
  BOOST_CHECK(noiseAxis2[0] == 0.0);

  Float64 n1Array[] = {0.5};
  CSpectrumNoiseAxis noiseAxis3 = CSpectrumNoiseAxis(n1Array, 1);
  BOOST_CHECK(noiseAxis3.GetSamplesCount() == 1);
  BOOST_CHECK(noiseAxis3[0] == 0.5);

  CSpectrumNoiseAxis noiseAxis4 = CSpectrumNoiseAxis(2, 0.5);
  BOOST_CHECK(noiseAxis4.GetSamplesCount() == 2);
  BOOST_CHECK(noiseAxis4[0] == 0.5);
  BOOST_CHECK(noiseAxis4[1] == 0.5);

  TFloat64List n1Vector = {0.5, 1.5};
  CSpectrumNoiseAxis noiseAxis5 = CSpectrumNoiseAxis(n1Vector);
  BOOST_CHECK(noiseAxis5.GetSamplesCount() == 2);
  BOOST_CHECK(noiseAxis5[0] == 0.5);
  BOOST_CHECK(noiseAxis5[1] == 1.5);

  CSpectrumNoiseAxis noiseAxis6 = CSpectrumNoiseAxis(TFloat64List{0.5, 1.5});
  BOOST_CHECK(noiseAxis6.GetSamplesCount() == 2);
  BOOST_CHECK(noiseAxis6[0] == 0.5);
  BOOST_CHECK(noiseAxis6[1] == 1.5);

  // copy
  CSpectrumNoiseAxis noiseAxis2b(noiseAxis2);
  BOOST_CHECK(noiseAxis2b.GetSamplesCount() == 1);
  BOOST_CHECK(noiseAxis2b[0] == 0.0);

  // move
  CSpectrumNoiseAxis noiseAxis2c = std::move(noiseAxis2);
  BOOST_CHECK(noiseAxis2c.GetSamplesCount() == 1);
  BOOST_CHECK(noiseAxis2c[0] == 0.0);

  // construct from base CSpectrumAxis
  CSpectrumAxis axis;
  CSpectrumNoiseAxis noiseAxis2d(axis);
  BOOST_CHECK(noiseAxis2d.GetSamplesCount() == 0);

  CSpectrumNoiseAxis noiseAxis2e(std::move(axis));
  BOOST_CHECK(noiseAxis2e.GetSamplesCount() == 0);

  // SetSize
  noiseAxis.resize(2);
  BOOST_CHECK(noiseAxis.GetSamplesCount() == 2);

  noiseAxis.resize(3, 2);
  BOOST_CHECK(noiseAxis.GetSamplesCount() == 3);
  BOOST_CHECK(noiseAxis[0] == 0.0);
  BOOST_CHECK(noiseAxis[1] == 0.0);
  BOOST_CHECK(noiseAxis[2] == 2.0);

  // Invert
  noiseAxis.Invert();
  BOOST_CHECK(noiseAxis.GetSamplesCount() == 3);
  BOOST_CHECK(noiseAxis[0] == INFINITY);
  BOOST_CHECK(noiseAxis[1] == INFINITY);
  BOOST_CHECK(noiseAxis[2] == 0.5);

  // sum
  BOOST_CHECK_THROW(noiseAxis2 + noiseAxis4, AmzException);
  CSpectrumNoiseAxis noiseAxis7(TFloat64List{
      0.,
      3,
  });
  CSpectrumNoiseAxis sumAxis2(TFloat64List{2., 4.});
  sumAxis2 += noiseAxis7;
  BOOST_CHECK(sumAxis2[0] == 2.);
  BOOST_CHECK(sumAxis2[1] == 5.);

  sumAxis2 = CSpectrumNoiseAxis(TFloat64List{2., 4.});
  sumAxis2 = sumAxis2 + noiseAxis7;
  BOOST_CHECK(sumAxis2[0] == 2.);
  BOOST_CHECK(sumAxis2[1] == 5.);

  //////////// DEBUGING //////////////////
  CSpectrumAxis axis7(noiseAxis7);
  CSpectrumNoiseAxis noiseAxis8(TFloat64List{2., 4.});
  sumAxis2 = axis7 + noiseAxis8;
  BOOST_CHECK(sumAxis2[0] == 2.);
  BOOST_CHECK(sumAxis2[1] == 5.);

  CSpectrumAxis &axis7_ref(noiseAxis7);
  sumAxis2 = axis7_ref + noiseAxis8;
  BOOST_CHECK(sumAxis2[0] == 2.);
  BOOST_CHECK(sumAxis2[1] == 5.);

  // diff
  BOOST_CHECK_THROW(noiseAxis2 - noiseAxis4, AmzException);
  sumAxis2 = CSpectrumNoiseAxis(TFloat64List{2., 4.});
  sumAxis2 -= noiseAxis7;
  BOOST_CHECK(sumAxis2[0] == 2.);
  BOOST_CHECK(sumAxis2[1] == 5.);

  sumAxis2 = CSpectrumNoiseAxis(TFloat64List{2., 4.});
  sumAxis2 = sumAxis2 - noiseAxis7;
  BOOST_CHECK(sumAxis2[0] == 2.);
  BOOST_CHECK(sumAxis2[1] == 5.);

  // extract
  CSpectrumNoiseAxis noiseAxis2f = noiseAxis.extract(1, 2);
  BOOST_CHECK(noiseAxis2f.GetSamplesCount() == 2);
  BOOST_CHECK(noiseAxis2f[0] == INFINITY);
  BOOST_CHECK(noiseAxis2f[1] == 0.5);

  // checkNoise
  TBoolList isValid{false, false, true};
  BOOST_CHECK(noiseAxis.checkNoise() == isValid);
  isValid = {true, true, true};
  noiseAxis[0] = 0.5;
  noiseAxis[1] = 0.5;
  BOOST_CHECK(noiseAxis.checkNoise() == isValid);
}

BOOST_AUTO_TEST_SUITE_END()