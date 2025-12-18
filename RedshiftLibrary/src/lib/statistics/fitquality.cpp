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
#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/operator/continuumfitting.h"
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

Float64 chi2(const TFloat64List &residuals) {
  return std::accumulate(
      residuals.begin(), residuals.end(), 0.0,
      [](Float64 sum, Float64 val) { return sum + val * val; });
}

Float64 reducedChi2(const Float64 chi2, const Int32 nPixels) {
  if (nPixels <= 0) {
    Flag.warning(WarningCode::TOO_LITTLE_PIXELS,
                 Formatter()
                     << "          NSFitQuality::" << __func__
                     << ": nPixels must be > 0 for reduced chi2 computation");
    return NAN;
  }
  return chi2 / nPixels;
}

Float64 pValue(const Float64 chi2, const Int32 nPixels) {
  if (nPixels <= 1) {
    Flag.warning(WarningCode::TOO_LITTLE_PIXELS,
                 Formatter() << "          NSFitQuality::" << __func__
                             << ": nPixels must be > 2 for pValue computation");
    return NAN;
  }
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

Float64 kurtosisGsl(const TFloat64List &data, const Float64 mean,
                    const Float64 stdev) {
  // Fisher-Pearson kurtosis measurement
  const Int32 n = data.size();
  if (n < 1)
    return NAN;
  return gsl_stats_kurtosis_m_sd(data.data(), 1, n, mean, stdev);
}

Float64 andersonDarlingTest(const TFloat64List &data, Float64 mean,
                            Float64 stdev) {
  // NB this method can also take mean and std if needed
  const Int32 n = data.size();
  if (n < 2)
    return NAN;
  return anderson_darling_normality_statistic(data, mean, stdev);
}

Float64 ksTest(const TFloat64List &data, const Float64 mean,
               const Float64 stdev, const bool sorted) {
  if (data.size() < 2)
    return NAN;
  auto dataBis = data;
  auto empiricalCdf =
      empirical_cumulative_distribution_function(std::move(dataBis), sorted);
  boost::math::normal dist(mean, stdev);

  Float64 maxDiff = 0.0;
  for (Int32 i = 0; i < ssize(data); i++) {
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

TFitQuality computeFitQuality(TFloat64List spcFlux, TFloat64List modelFlux,
                              TFloat64List spcFluxError, Float64 chi2,
                              Int32 nPixelsOfResiduals,
                              Int32 nPixelsUsedForFit) {
  // kEnd set to -1 means take the full spectrum
  if (ssize(spcFlux) != ssize(modelFlux) ||
      ssize(spcFlux) != ssize(spcFluxError)) {
    THROWG(ErrorCode::INTERNAL_ERROR, "spcFlux, modelFlux and "
                                      "spcFluxError must be of the same size");
  }

  TList<TFloat64List> spcFluxVect;
  TList<TFloat64List> modelFluxVect;
  TList<TFloat64List> spcFluxErrorVect;
  spcFluxVect.push_back(std::move(spcFlux));
  modelFluxVect.push_back((std::move(modelFlux)));
  spcFluxErrorVect.push_back(std::move(spcFluxError));

  return computeFitQuality(spcFluxVect, modelFluxVect, spcFluxErrorVect, chi2,
                           nPixelsOfResiduals, nPixelsUsedForFit);
}

TFitQuality computeFitQuality(const std::vector<TFloat64List> &spcFlux,
                              const std::vector<TFloat64List> &modelFlux,
                              const std::vector<TFloat64List> &spcFluxError,
                              Float64 chi2, Int32 nPixelsOfResiduals,
                              Int32 nPixelsUsedForFit,
                              const std::vector<CMask> &mask) {
  // It is expected that the input vectors are of the same size

  const bool useMask = mask.empty() ? false : true;

  if (ssize(spcFlux) != ssize(modelFlux) ||
      ssize(spcFlux) != ssize(spcFluxError) ||
      (useMask && ssize(spcFlux) != ssize(mask))) {
    THROWG(ErrorCode::INTERNAL_ERROR,
           "m_spectra, spcFlux, modelFlux,  "
           "spcFluxError and mask must be of the same size");
    for (Int32 spcIdx = 0; spcIdx < ssize(spcFlux); spcIdx++) {
      if (ssize(spcFlux[spcIdx]) != ssize(modelFlux[spcIdx]) ||
          ssize(spcFlux[spcIdx]) != ssize(spcFluxError[spcIdx]) ||
          (useMask && ssize(spcFlux[spcIdx]) != mask[spcIdx].GetMasksCount())) {
        THROWG(ErrorCode::INTERNAL_ERROR,
               Formatter()
                   << "spcFlux, modelFlux, spcFluxError and mask must be of "
                      "the same size at spectrum index "
                   << spcIdx);
      }
    }
  }

  const Int32 nSpectra = ssize(spcFlux);
  // Compute the maximum number of pixels used to compute residuals in order to
  // reserve enough space in vector
  Int32 nTotPixels = 0;
  for (auto const &vect : spcFlux)
    nTotPixels += vect.size();

  const auto isValid = [useMask](const std::vector<CMask> &mask, Int32 spcIdx,
                                 Int32 pixelIdx) {
    return !useMask || mask[spcIdx][pixelIdx];
  };
  TFloat64List residuals;
  residuals.reserve(nTotPixels);
  Int32 sumNPixels = 0;
  for (Int32 spcIdx = 0; spcIdx != nSpectra; spcIdx++) {
    for (Int32 pixelIdx = 0; pixelIdx != ssize(spcFlux[spcIdx]); pixelIdx++) {
      if (isValid(mask, spcIdx, pixelIdx)) {
        const Float64 residual = NSFitQuality::computeResidual(
            spcFlux[spcIdx][pixelIdx], modelFlux[spcIdx][pixelIdx],
            spcFluxError[spcIdx][pixelIdx]);
        if (std::isnan(residual))
          continue;
        residuals.push_back(residual);
        sumNPixels += 1;
      }
    }
  }
  std::sort(residuals.begin(), residuals.end());

  if (std::isnan(chi2))
    chi2 = NSFitQuality::chi2(residuals);
  else {
    // if chi2 is given as input then nPixelsOfResiduals should be given as
    // well.
    if (nPixelsOfResiduals == undefIdx)
      THROWG(ErrorCode::INTERNAL_ERROR,
             "undefined nPixelsOfResiduals argument, it should "
             "be given since chi2 was given");
    sumNPixels = nPixelsOfResiduals;
  }

  TFitQuality fitQuality;
  fitQuality.reducedChiSquare = NSFitQuality::reducedChi2(chi2, sumNPixels);
  fitQuality.pValue = NSFitQuality::pValue(chi2, sumNPixels);

  fitQuality.meanResiduals = NSFitQuality::mean(residuals);
  fitQuality.stdResiduals =
      NSFitQuality::stdev(residuals, fitQuality.meanResiduals);
  fitQuality.skewnessResiduals = NSFitQuality::skewness(
      residuals, fitQuality.meanResiduals, fitQuality.stdResiduals);
  fitQuality.kurtosisResiduals = NSFitQuality::kurtosisGsl(
      residuals, fitQuality.meanResiduals, fitQuality.stdResiduals);
  fitQuality.ksResiduals = NSFitQuality::ksTest(residuals, 0, 1, true);
  fitQuality.ksStdResiduals = NAN;
  fitQuality.ksStdMeanResiduals = NAN;
  if (fitQuality.stdResiduals > DBL_MIN) {
    fitQuality.ksStdResiduals =
        NSFitQuality::ksTest(residuals, 0, fitQuality.stdResiduals, true);
    fitQuality.ksStdMeanResiduals = NSFitQuality::ksTest(
        residuals, fitQuality.meanResiduals, fitQuality.stdResiduals, true);
  }
  fitQuality.andersonResiduals = NSFitQuality::andersonDarlingTest(
      residuals, fitQuality.meanResiduals, fitQuality.stdResiduals);
  fitQuality.nPixelsOfResiduals = sumNPixels;
  fitQuality.nPixelsUsedForFit = nPixelsUsedForFit;
  return fitQuality;
}

} // namespace NSEpic::NSFitQuality
