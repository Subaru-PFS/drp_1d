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
#include <boost/range/combine.hpp>

#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/operator/continuumfitting.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/statistics/fitquality.h"

using namespace NSEpic;
using namespace std;

COperatorContinuumFitting::COperatorContinuumFitting()
    : m_kStart(Context.getSpectra().size()),
      m_kEnd(Context.getSpectra().size()),
      m_maskBuilder(std::make_shared<CMaskBuilder>()),
      m_spectra(Context.getSpectra()),
      m_lambdaRanges(Context.getClampedLambdaRanges()){};

/**
 * \brief this function estimates the likelihood_cstLog term withing the
 * wavelength range
 **/
Float64 COperatorContinuumFitting::EstimateLikelihoodCstLog() const {
  Float64 cstLog = 0.0;
  for (auto const &[spectrum_ptr, lambdaRange_ptr] :
       boost::combine(m_spectra, m_lambdaRanges)) {
    const CSpectrumSpectralAxis &spcSpectralAxis =
        spectrum_ptr->GetSpectralAxis();
    const TFloat64List &error =
        spectrum_ptr->GetFluxAxis().GetError().GetSamplesVector();

    Int32 numDevs = 0;

    Float64 sumLogNoise = 0.0;

    Int32 imin;
    Int32 imax;
    lambdaRange_ptr->getClosedIntervalIndices(
        spcSpectralAxis.GetSamplesVector(), imin, imax);
    for (Int32 j = imin; j <= imax; j++) {
      numDevs++;
      sumLogNoise += log(error[j]);
    }
    cstLog += -numDevs * 0.5 * log(2 * M_PI) - sumLogNoise;
  }
  return cstLog;
}

Float64
COperatorContinuumFitting::computeNPixels(const Int32 spcIdx,
                                          const TInt32List &kStart,
                                          const TInt32List &kEnd) const {
  return kEnd[spcIdx] - kStart[spcIdx] + 1;
}

void COperatorContinuumFitting::addQualityFitResidualsToResult(
    TContinuumResult &result, const TFloat64List &spcFlux,
    const TFloat64List &modelFlux, const TFloat64List &spcFluxError,
    const Int32 kStart, Int32 kEnd) const {
  // kEnd set to -1 means take the full spectrum
  if (ssize(spcFlux) != ssize(modelFlux) ||
      ssize(spcFlux) != ssize(spcFluxError)) {
    THROWG(ErrorCode::INTERNAL_ERROR, "m_spectra, spcFlux, modelFlux and "
                                      "spcFluxError must be of the same size");
  }
  const Int32 nPixels = ssize(spcFlux);
  if (kEnd == -1)
    kEnd = nPixels - 1;

  addQualityFitResidualsToResult(result, std::vector<TFloat64List>(1, spcFlux),
                                 std::vector<TFloat64List>(1, modelFlux),
                                 std::vector<TFloat64List>(1, spcFluxError),
                                 TInt32List(1, kStart), TInt32List(1, kEnd));
}

void COperatorContinuumFitting::addQualityFitResidualsToResult(
    TContinuumResult &result, const std::vector<TFloat64List> &spcFlux,
    const std::vector<TFloat64List> &modelFlux,
    const std::vector<TFloat64List> &spcFluxError, const TInt32List &kStartArg,
    const TInt32List &kEndArg) const {
  // It is expected that the input vectors are of the same size

  const auto &kStart = kStartArg.empty() ? m_kStart : kStartArg;
  const auto &kEnd = kEndArg.empty() ? m_kEnd : kEndArg;

  if (ssize(spcFlux) != ssize(modelFlux) || ssize(spcFlux) != ssize(kStart) ||
      ssize(spcFlux) != ssize(kEnd) || ssize(spcFlux) != ssize(spcFluxError)) {
    THROWG(ErrorCode::INTERNAL_ERROR, "m_spectra, spcFlux, modelFlux and "
                                      "spcFluxError must be of the same size");
    for (Int32 spcIdx = 0; spcIdx < ssize(spcFlux); spcIdx++) {
      if (ssize(spcFlux[spcIdx]) != ssize(modelFlux[spcIdx]) ||
          ssize(spcFlux[spcIdx]) != ssize(spcFluxError[spcIdx])) {
        THROWG(ErrorCode::INTERNAL_ERROR,
               Formatter() << "spcFlux, modelFlux and spcFluxError must be of "
                              "the same size at spectrum index "
                           << spcIdx);
      }
    }
  }

  const Int32 nSpectra = ssize(spcFlux);
  // Compute the maximum number of pixels used to compute residuals in order to
  // reserve enough space in vector
  // TODO see if can do this better
  Int32 nTotPixels = 0;
  for (Int32 spcIdx = 0; spcIdx < nSpectra; spcIdx++) {
    nTotPixels += computeNPixels(spcIdx, kStart, kEnd);
    ;
  }

  TFloat64List residuals;
  residuals.reserve(nTotPixels);
  Int32 sumNPixels = 0;
  for (Int32 spcIdx = 0; spcIdx < nSpectra; spcIdx++) {
    for (Int32 pixelIdx = kStart[spcIdx]; pixelIdx <= kEnd[spcIdx];
         pixelIdx++) {
      const Float64 residual = NSFitQuality::computeResidual(
          spcFlux[spcIdx][pixelIdx], modelFlux[spcIdx][pixelIdx],
          spcFluxError[spcIdx][pixelIdx]);
      if (std::isnan(residual))
        continue;
      residuals.push_back(residual);
      sumNPixels += 1;
    }
  }

  std::sort(residuals.begin(), residuals.end());
  result.fitQuality.meanResiduals = NSFitQuality::mean(residuals);
  result.fitQuality.stdResiduals =
      NSFitQuality::stdev(residuals, result.fitQuality.meanResiduals);
  result.fitQuality.skewnessResiduals =
      NSFitQuality::skewness(residuals, result.fitQuality.meanResiduals,
                             result.fitQuality.stdResiduals);
  result.fitQuality.kurtosisResiduals =
      NSFitQuality::kurtosis(residuals, result.fitQuality.meanResiduals);
  result.fitQuality.ksResiduals = NSFitQuality::ksTest(residuals, 0, 1, true);
  result.fitQuality.ksStdResiduals = NAN;
  result.fitQuality.ksStdMeanResiduals = NAN;
  if (result.fitQuality.stdResiduals > DBL_MIN) {
    result.fitQuality.ksStdResiduals = NSFitQuality::ksTest(
        residuals, 0, result.fitQuality.stdResiduals, true);
    result.fitQuality.ksStdMeanResiduals =
        NSFitQuality::ksTest(residuals, result.fitQuality.meanResiduals,
                             result.fitQuality.stdResiduals, true);
  }
  result.fitQuality.andersonResiduals =
      NSFitQuality::andersonDarlingTest(residuals);
  result.fitQuality.nPixels = sumNPixels;
}
