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

#include "RedshiftLibrary/operator/curve.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"

using namespace NSEpic;

TCurve::TCurve() {}

TCurve::TCurve(TList<Float64> lambda, TList<Float64> flux,
               TList<Float64> fluxError) {
  setLambda(lambda);
  setFlux(flux);
  setFluxError(fluxError);
}

void TCurve::checkIdx(Int32 idx) const {
  if (idx > size())
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "TCurve::" << __func__ << ": Trying to get index "
                       << idx << "outside of curve with size " << size());
}

TCurveElement TCurve::get_at_index(Int32 idx) const {
  checkIdx(idx);
  return {lambda[idx], flux[idx], fluxError[idx]};
}

Float64 TCurve::getLambdaAt(Int32 pixelIdx) const {
  checkIdx(pixelIdx);
  return lambda[pixelIdx];
};
Float64 TCurve::getFluxAt(Int32 pixelIdx) const {
  checkIdx(pixelIdx);
  return flux[pixelIdx];
};
Float64 TCurve::getFluxErrorAt(Int32 pixelIdx) const {
  checkIdx(pixelIdx);
  return fluxError[pixelIdx];
};

void TCurve::push_back(TCurveElement const &elem) {
  lambda.push_back(elem.lambda);
  flux.push_back(elem.flux);
  fluxError.push_back(elem.fluxError);
}

Int32 TCurve::size() const { return lambda.size(); }

void TCurve::setLambda(TFloat64List inputLambda) {
  lambda = std::move(inputLambda);
}

void TCurve::setFlux(TList<Float64> inputFlux) { flux = std::move(inputFlux); }

void TCurve::setFluxError(TList<Float64> inputFluxError) {
  fluxError = std::move(inputFluxError);
}

void TCurve::reserve(Int32 size) {
  lambda.reserve(size);
  flux.reserve(size);
  fluxError.reserve(size);
}

void TCurve::sort() {
  // Reorder lambda, and move flux and fluxError elements accordingly

  if (CSpectrumSpectralAxis(lambda).isSorted())
    return;

  TList<Int32> sortingIndices(lambda.size());
  for (size_t pixelIdx = 0; pixelIdx < sortingIndices.size(); pixelIdx++) {
    sortingIndices[pixelIdx] = pixelIdx;
  }

  std::sort(sortingIndices.begin(), sortingIndices.end(),
            [this](size_t i1, size_t i2) { return lambda[i1] < lambda[i2]; });

  // Use the sorted indices to reorder lambda, flux, and fluxError
  std::vector<Float64> lambdaSorted(lambda.size()), FluxSorted(flux.size()),
      fluxErrorSorted(fluxError.size());
  for (std::size_t i = 0; i < sortingIndices.size(); ++i) {
    lambdaSorted[i] = lambda[sortingIndices[i]];
    FluxSorted[i] = flux[sortingIndices[i]];
    fluxErrorSorted[i] = fluxError[sortingIndices[i]];
  }

  // Assign the sorted vectors back to the original vectors
  lambda = std::move(lambdaSorted);
  flux = std::move(FluxSorted);
  fluxError = std::move(fluxErrorSorted);
}