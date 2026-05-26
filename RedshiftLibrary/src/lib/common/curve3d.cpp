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
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================

#include "RedshiftLibrary/common/curve3d.h"
#include "RedshiftLibrary/common/curve.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/size.h"

using namespace NSEpic;

T3DCurve::T3DCurve(Int32 nIgm, Int32 nIsm) : nIgm(nIgm), nIsm(nIsm){};

T3DCurve::T3DCurve(TCurve &&curve, Int32 nIgm_, Int32 nIsm_)
    : TCurve(std::move(curve)), nIgm(nIgm_), nIsm(nIsm_),
      flux(T3DList<Float64>(1, T2DList<Float64>(1))),
      fluxError(T3DList<Float64>(1, T2DList<Float64>(1))),
      isExtincted(T3DList<bool>(1, T2DList<bool>(1))) {

  flux[0][0] = std::move(TCurve::flux);
  fluxError[0][0] = std::move(TCurve::fluxError);
  isExtincted[0][0] = std::move(TCurve::isExtincted);

  flux[0].resize(nIsm, flux[0][0]);
  flux.resize(nIgm, flux[0]);

  fluxError[0].resize(nIsm, fluxError[0][0]);
  fluxError.resize(nIgm, fluxError[0]);

  isExtincted[0].resize(nIsm, isExtincted[0][0]);
  isExtincted.resize(nIgm, isExtincted[0]);
};

template <typename T>
void T3DCurve::checkCurveElement(T3DList<T> const &inputElement,
                                 std::string elementName) {

  Int32 dim1 = inputElement.size();
  Int32 dim2;

  for (Int32 i = 0; i < dim1; i++) {
    dim2 = inputElement[i].size();
    if (dim1 != nIgm || dim2 != nIsm) {
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "Wrong curve " << elementName << " size, should be "
                         << nIgm << ", " << nIsm << " but is " << dim1 << ", "
                         << dim2);
    }
  }
}

void T3DCurve::setFlux(T3DList<Float64> inputFlux) {
  checkCurveElement(inputFlux, "flux");
  flux = std::move(inputFlux);
}

void T3DCurve::setFluxError(T3DList<Float64> inputFluxError) {
  checkCurveElement(inputFluxError, "fluxError");
  fluxError = std::move(inputFluxError);
}

void T3DCurve::setIsExtincted(T3DList<bool> inputIsExtincted) {
  checkCurveElement(inputIsExtincted, "isExtincted");
  isExtincted = std::move(inputIsExtincted);
}

TCurve &&T3DCurve::toCurve(Int16 igmIdx, Int16 ismIdx) && {
  if (ismIdx > nIsm)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "ismIdx = " << ismIdx << " < " << nIsm);
  if (igmIdx > nIgm)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "ismIdx = " << igmIdx << " < " << nIgm);
  TCurve &curve = *this;
  curve.setFlux(std::move(flux[igmIdx][ismIdx]));
  curve.setFluxError(std::move(fluxError[igmIdx][ismIdx]));
  curve.setIsExtincted(std::move(isExtincted[igmIdx][ismIdx]));
  return std::move(curve);
}

TCurve T3DCurve::toCoefCurve(Int16 igmIdx, Int16 ismIdx) const {
  if (ismIdx > nIsm)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "ismIdx = " << ismIdx << " < " << nIsm);
  checkIgmIdx(igmIdx);
  checkIsmIdx(ismIdx);
  TCurve curve;
  curve.reserve(size());
  for (Int32 pixelIdx = 0; pixelIdx < size(); pixelIdx++) {
    if (pixelIsCoefValid(igmIdx, ismIdx, pixelIdx))
      curve.push_back({lambda[pixelIdx], flux[igmIdx][ismIdx][pixelIdx],
                       fluxError[igmIdx][ismIdx][pixelIdx]});
  }
  return curve;
}

bool T3DCurve::pixelIsCoefValid(Int16 igmIdx, Int16 ismIdx,
                                Int32 pixelIdx) const {
  return pixelIdx < size() && isSnrCompliant[pixelIdx] && mask[pixelIdx] &&
         !isExtincted[igmIdx][ismIdx][pixelIdx];
}

void T3DCurve::checkIgmIdx(Int16 igmIdx) const {
  if (igmIdx > nIgm)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "igmIdx = " << igmIdx << " > " << nIgm);
}

void T3DCurve::checkIsmIdx(Int16 ismIdx) const {
  if (ismIdx > nIsm)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "ismIdx = " << ismIdx << " > " << nIsm);
}

void T3DCurve::checkIdxs(Int16 igmIdx, Int16 ismIdx, Int32 pixelIdx) const {
  checkIgmIdx(igmIdx);
  checkIsmIdx(ismIdx);
  checkIdx(pixelIdx);
}

void T3DCurve::setFluxAt(Int16 igmIdx, Int16 ismIdx, Int32 pixelIdx,
                         Float64 value) {
  checkIdxs(igmIdx, ismIdx, pixelIdx);
  flux[igmIdx][ismIdx][pixelIdx] = value;
};

void T3DCurve::setFluxErrorAt(Int16 igmIdx, Int16 ismIdx, Int32 pixelIdx,
                              Float64 value) {
  checkIdxs(igmIdx, ismIdx, pixelIdx);
  fluxError[igmIdx][ismIdx][pixelIdx] = value;
};

Float64 T3DCurve::getFluxAt(Int16 igmIdx, Int16 ismIdx, Int32 pixelIdx) const {
  checkIdxs(igmIdx, ismIdx, pixelIdx);
  return flux[igmIdx][ismIdx][pixelIdx];
};

Float64 T3DCurve::getFluxErrorAt(Int16 igmIdx, Int16 ismIdx,
                                 Int32 pixelIdx) const {
  checkIdxs(igmIdx, ismIdx, pixelIdx);
  return fluxError[igmIdx][ismIdx][pixelIdx];
};

bool T3DCurve::getIsExtinctedAt(Int16 igmIdx, Int16 ismIdx,
                                Int32 pixelIdx) const {
  checkIdxs(igmIdx, ismIdx, pixelIdx);
  return isExtincted[igmIdx][ismIdx][pixelIdx];
};