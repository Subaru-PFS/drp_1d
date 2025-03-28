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

#include "RedshiftLibrary/statistics/fitquality.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/exception.h"
#include <boost/accumulators/statistics/skewness.hpp>
#include <boost/math/distributions/empirical_cumulative_distribution_function.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/statistics/anderson_darling.hpp>
#include <gsl/gsl_statistics_double.h>

using boost::math::cdf;
using boost::math::complement;
using boost::math::empirical_cumulative_distribution_function;
using boost::math::statistics::anderson_darling_normality_statistic;

namespace NSEpic::NSFitQuality {
Float64 reducedChi2(const Float64 chi2, const Int32 nPixels) {
  if (nPixels <= 0)
    THROWG(ErrorCode::INTERNAL_ERROR, "nPixels must be > 0 ");
  return chi2 / nPixels;
}

Float64 pValue(const Float64 chi2, const Int32 nPixels) {
  if (nPixels <= 1)
    THROWG(ErrorCode::INTERNAL_ERROR, "nPixels must be > 2 ");
  if (chi2 > DBL_MAX)
    return 0;
  boost::math::chi_squared chi2Dist(nPixels - 1);
  Float64 p = cdf(complement(chi2Dist, chi2));
  return p;
}

Float64 mean(const TFloat64List &data) {
  const Int32 n = data.size();
  if (n == 0)
    return NAN;
  return gsl_stats_mean(data.data(), 1, n);
}

Float64 var(const TFloat64List &data, const Float64 mean) {
  const Int32 n = data.size();
  if (n < 2)
    return NAN;
  return gsl_stats_variance_m(data.data(), 1, n, mean);
}

Float64 stdev(const TFloat64List &data, const Float64 mean) {
  return std::sqrt(var(data, mean));
}

Float64 skewness(const TFloat64List &data, const Float64 mean,
                 const Float64 stdev) {
  Int32 n = data.size();
  if (n < 1)
    return NAN;
  return gsl_stats_skew_m_sd(data.data(), 1, n, mean, stdev);
}

Float64 kurtosis(const TFloat64List &data, const Float64 mean) {
  // Fisher-Pearson kurtosis measurement
  const Int32 n = data.size();
  if (n < 1)
    return NAN;
  Float64 invN = 1.0 / n;
  // Compute second and fourth central moments
  Float64 sum2 = 0.0;
  Float64 sum4 = 0.0;
  for (Float64 r : data) {
    Float64 diff = r - mean;
    Float64 diff2 = diff * diff;
    sum2 += diff2;
    sum4 += diff2 * diff2;
  }

  // Compute kurtosis excess
  Float64 invNSum2 = invN * sum2;
  return (invN * sum4) / (invNSum2 * invNSum2) - 3.0;
}

Float64 kurtosisGsl(const TFloat64List &data, const Float64 mean,
                    const Float64 stdev) {
  // Fisher-Pearson kurtosis measurement
  const Int32 n = data.size();
  if (n < 1)
    return NAN;
  return gsl_stats_kurtosis_m_sd(data.data(), 1, n, mean, stdev);
}

Float64 andersonDarlingTest(const TFloat64List data) {
  // NB this method can also take mean and std if needed
  return anderson_darling_normality_statistic(data);
}

Float64 ksTest(const TFloat64List &data, const Float64 mean,
               const Float64 stdev, const bool sorted) {
  auto dataBis = data;
  auto empiricalCdf =
      empirical_cumulative_distribution_function(std::move(dataBis), sorted);
  boost::math::normal dist(mean, stdev);

  Float64 maxDiff = 0.0;
  for (Int32 i = 0; i < data.size(); i++) {
    const Float64 diff =
        std::abs(empiricalCdf(data[i]) - boost::math::cdf(dist, data[i]));
    if (diff > maxDiff)
      maxDiff = diff;
  }

  return maxDiff;
}

Float64 computeResidual(const Float64 expData, const Float64 refData,
                        const Float64 expDataError) {
  if (expDataError == 0.0)
    return NAN;
  return (expData - refData) / expDataError;
};

} // namespace NSEpic::NSFitQuality
