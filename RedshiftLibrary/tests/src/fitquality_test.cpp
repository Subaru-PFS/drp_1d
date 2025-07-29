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
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/statistics/fitquality.h"
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(FitQuality_test)
BOOST_AUTO_TEST_CASE(reducedChi2_test) {
  BOOST_CHECK_NO_THROW(NSFitQuality::reducedChi2(10., 0));
  BOOST_CHECK_CLOSE(NSFitQuality::reducedChi2(100., 10), 10., 1e-4);
}

BOOST_AUTO_TEST_CASE(pValue_test) {
  BOOST_CHECK_NO_THROW(NSFitQuality::pValue(10., 1));
  BOOST_CHECK_CLOSE(NSFitQuality::pValue(5., 4), 0.1717971442967335, 1e-4);
}

BOOST_AUTO_TEST_CASE(mean_test) {
  BOOST_CHECK(std::isnan(NSFitQuality::mean({})));
  BOOST_CHECK_CLOSE(NSFitQuality::mean({1., 2., 3.}), 2., 1e-4);
}

BOOST_AUTO_TEST_CASE(stdev_test) {
  BOOST_CHECK(std::isnan(NSFitQuality::stdev({1}, 2.)));
  BOOST_CHECK_CLOSE(NSFitQuality::stdev({1., 2., 3.}, 2.), 1., 1e-4);
}

BOOST_AUTO_TEST_CASE(skewness_test) {
  const auto data = {0., 1., 3., 3., 5., 7., 10.};
  const auto mean = NSFitQuality::mean(data);
  const auto stdev = NSFitQuality::stdev(data, mean);

  BOOST_CHECK(std::isnan(NSFitQuality::skewness({}, mean, stdev)));
  BOOST_CHECK_CLOSE(NSFitQuality::skewness(data, mean, stdev),
                    0.40431016151273425, 1e-4);
}

std::vector<double> linspace(double start, double end, std::size_t num) {
  std::vector<double> result;
  result.reserve(num);

  if (num == 0)
    return result;
  if (num == 1) {
    result.push_back(start);
    return result;
  }

  double step = (end - start) / (num - 1);
  for (std::size_t i = 0; i < num; ++i) {
    result.push_back(start + step * i);
  }

  return result;
}

BOOST_AUTO_TEST_CASE(kurtosis_test) {
  const auto data = {0., 1., 3., 3., 5., 7., 10.};
  const auto mean = NSFitQuality::mean(data);
  const auto stdev = NSFitQuality::stdev(data, mean);

  BOOST_CHECK(std::isnan(NSFitQuality::kurtosisGsl({}, mean, stdev)));
  BOOST_CHECK_CLOSE(NSFitQuality::kurtosisGsl(data, mean, stdev),
                    -1.41141727279147, 1e-4);
  const auto linspaceArray = linspace(0, 1, 500);
  const auto linspaceArrayMean = NSFitQuality::mean(linspaceArray);
  const auto linspaceArrayStd =
      NSFitQuality::stdev(linspaceArray, linspaceArrayMean);
  BOOST_CHECK_CLOSE(NSFitQuality::kurtosisGsl(linspaceArray, linspaceArrayMean,
                                              linspaceArrayStd),
                    -1.2072023616766467, 1e-4);
}

BOOST_AUTO_TEST_CASE(andersonDarlingTest_test) {
  const auto data = {0., 1., 3., 3., 5., 7., 10.};
  const auto mean = NSFitQuality::mean(data);
  const auto stdev = NSFitQuality::stdev(data, mean);
  BOOST_CHECK_CLOSE(NSFitQuality::andersonDarlingTest(data, mean, stdev),
                    0.22083660833332985, 1e-4);
}

BOOST_AUTO_TEST_CASE(ksTest_test) {
  const auto data = {-0.2, -0.1, 0., 0., 0.1, 0.2};
  const auto mean = NSFitQuality::mean(data);
  const auto stdev = NSFitQuality::stdev(data, mean);
  BOOST_CHECK_CLOSE(NSFitQuality::ksTest(data, 0, 1, true), 0.42074029056089701,
                    1e-4);
  BOOST_CHECK_CLOSE(NSFitQuality::ksTest(data, mean, stdev, true),
                    0.16666666666666663, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()