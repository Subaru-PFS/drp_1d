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
#ifndef _REDSHIFT_CURVE_3D_
#define _REDSHIFT_CURVE_3D_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/operator/powerlawbase.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"

// TODO clean here

namespace NSEpic {

struct TCurveElement {
  Float64 lambda;
  Float64 flux;
  bool pixelToUse;
  Float64 fluxError;
};

struct TCurve {
  TCurve(Int32 nMaxPixels) : nMaxPixels(nMaxPixels) {
    lambda.reserve(nMaxPixels);
    flux.reserve(nMaxPixels);
    pixelsToUse.reserve(nMaxPixels);
    fluxError.reserve(nMaxPixels);
  };

  Int32 nMaxPixels;
  TAxisSampleList lambda;
  TFloat64List flux;
  TBoolList pixelsToUse;
  TFloat64List fluxError;

  // Proxy class to handle assignment at a specific index
  // struct Proxy {
  //   TCurve &curve;
  //   Int32 idx;

  //   // Overload the assignment operator
  //   Proxy &operator=(const Proxy &other) {
  //     curve.lambda[idx] = other.curve.lambda[other.idx];
  //     curve.flux[idx] = other.curve.flux[other.idx];
  //     curve.pixelsToUse[idx] = other.curve.pixelsToUse[other.idx];
  //     curve.fluxError[idx] = other.curve.fluxError[other.idx];
  //     return *this;
  //   }

  //   Proxy &operator=(const Proxy &other) const {
  //     curve.lambda[idx] = other.curve.lambda[other.idx];
  //     curve.flux[idx] = other.curve.flux[other.idx];
  //     curve.pixelsToUse[idx] = other.curve.pixelsToUse[other.idx];
  //     curve.fluxError[idx] = other.curve.fluxError[other.idx];
  //     return *const_cast<Proxy *>(this);
  //   }
  // };

  // Overload the operator[] to return a proxy object
  // Proxy operator[](Int32 idx) { return Proxy{*this, idx}; }

  TCurveElement operator[](Int32 idx) {
    return {lambda[idx], flux[idx], pixelsToUse[idx], fluxError[idx]};
  }

  const TCurveElement operator[](Int32 idx) const {
    return {lambda[idx], flux[idx], pixelsToUse[idx], fluxError[idx]};
  }
  // const Proxy operator[](Int32 idx) const {
  //   return Proxy{*const_cast<TCurve *>(this), idx};
  // }

  void push_back(TCurveElement const &elem) {
    // TODO add an error if size too big

    lambda.push_back(elem.lambda);
    flux.push_back(elem.flux);
    pixelsToUse.push_back(elem.pixelToUse);
    fluxError.push_back(elem.fluxError);
  }

  Int32 size() const { return lambda.size(); }

  void setLambda(TFloat64List const &inputLambda) {
    if (inputLambda.size() > nMaxPixels) {
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "T3DCurve::setLambda Lambda size toolarge, should "
                            "be lower than"
                         << nMaxPixels << " but is " << inputLambda.size());
    }
    lambda = inputLambda;
  }

  template <typename T>
  void checkCurveElement(TList<T> const &inputElement,
                         std::string elementName) {
    if (inputElement.size() > nMaxPixels) {
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "Wrong curve " << elementName
                         << " size, should be lower than " << nMaxPixels
                         << " but is " << inputElement.size());
    }
  }

  void setFlux(TList<Float64> const &inputFlux) {
    checkCurveElement(inputFlux, "flux");
    flux = inputFlux;
  }

  void setPixelsToUse(TList<bool> const &inputPixelsToUse) {
    checkCurveElement(inputPixelsToUse, "pixelsToUse");
    pixelsToUse = inputPixelsToUse;
  }

  void setFluxError(TList<Float64> const &inputFluxError) {
    checkCurveElement(inputFluxError, "fluxError");
    fluxError = inputFluxError;
  }
};

struct T3DCurve {
  T3DCurve(Int32 npixels, Int32 nigm, Int32 nism)
      : nPixels(npixels), nIgm(nigm), nIsm(nism){};

  // rule of 5 defaults
  ~T3DCurve() = default;
  T3DCurve(const T3DCurve &) = default;
  T3DCurve(T3DCurve &&) = default;
  T3DCurve &operator=(const T3DCurve &) = default;
  T3DCurve &operator=(T3DCurve &&) = default;

  Int32 nPixels;
  Int32 nIgm;
  Int32 nIsm;
  TList<Float64> lambda = TList<Float64>(nPixels);
  T3DList<Float64> flux =
      T3DList<Float64>(nIgm, T2DList<Float64>(nIsm, TList<Float64>(nPixels)));
  T3DList<bool> pixelsToUse =
      T3DList<bool>(nIgm, T2DList<bool>(nIsm, TList<bool>(nPixels)));
  T3DList<Float64> fluxError =
      T3DList<Float64>(nIgm, T2DList<Float64>(nIsm, TList<Float64>(nPixels)));

  friend T3DCurve operator+(const T3DCurve &c1, const T3DCurve &c2) {
    // NB after addition lambda is unsorted
    if (c1.nIgm != c2.nIgm)
      THROWG(
          ErrorCode::INTERNAL_ERROR,
          Formatter() << "T3DCurve::+ : Impossible to add curves with nIgm = "
                      << c1.nIgm << "and " << c2.nIgm);
    if (c1.nIsm != c2.nIsm)
      THROWG(
          ErrorCode::INTERNAL_ERROR,
          Formatter() << "T3DCurve::+ : Impossible to add curves with nIsm = "
                      << c1.nIsm << "and " << c2.nIsm);
    Int32 newNPixels = c1.nPixels + c2.nPixels;

    T3DCurve newCurve(newNPixels, c1.nIgm, c1.nIsm);

    for (Int32 pixelIdx = 0; pixelIdx < c1.nPixels; pixelIdx++) {
      newCurve.lambda[pixelIdx] = c1.lambda[pixelIdx];
      for (Int16 igmIdx = 0; igmIdx < c1.nIgm; igmIdx++) {
        for (Int16 ismIdx = 0; ismIdx < c1.nIsm; ismIdx++) {
          newCurve.flux[igmIdx][ismIdx][pixelIdx] =
              c1.flux[igmIdx][ismIdx][pixelIdx];
          newCurve.pixelsToUse[igmIdx][ismIdx][pixelIdx] =
              c1.pixelsToUse[igmIdx][ismIdx][pixelIdx];
          newCurve.fluxError[igmIdx][ismIdx][pixelIdx] =
              c1.fluxError[igmIdx][ismIdx][pixelIdx];
        }
      }
    }
    for (Int32 pixelIdx = 0; pixelIdx < c2.nPixels; pixelIdx++) {
      newCurve.lambda[c1.nPixels + pixelIdx] = c2.lambda[pixelIdx];
      for (Int16 igmIdx = 0; igmIdx < c1.nIgm; igmIdx++) {
        for (Int16 ismIdx = 0; ismIdx < c1.nIsm; ismIdx++) {
          newCurve.flux[igmIdx][ismIdx][c1.nPixels + pixelIdx] =
              c2.flux[igmIdx][ismIdx][pixelIdx];
          newCurve.pixelsToUse[igmIdx][ismIdx][c1.nPixels + pixelIdx] =
              c2.pixelsToUse[igmIdx][ismIdx][pixelIdx];
          newCurve.fluxError[igmIdx][ismIdx][c1.nPixels + pixelIdx] =
              c2.fluxError[igmIdx][ismIdx][pixelIdx];
        }
      }
    }

    return newCurve;
  }

  Int32 size() const { return nPixels; }

  void setLambda(TFloat64List const &inputLambda) {
    if (inputLambda.size() != nPixels) {
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "T3DCurve::setLambda Wrong lambda size, should be "
                         << nPixels << " but is " << inputLambda.size());
    }
    lambda = inputLambda;
  }

  template <typename T>
  void checkCurveElement(T3DList<T> const &inputElement,
                         std::string elementName) {

    Int32 dim1 = inputElement.size();
    Int32 dim2;
    Int32 dim3;

    for (Int32 i = 0; i < dim1; i++) {
      dim2 = inputElement[i].size();
      for (Int32 j = 0; j < dim1; j++) {
        dim3 = inputElement[i][j].size();
        if (dim1 != nIgm || dim2 != nIsm || dim3 != nPixels) {
          THROWG(ErrorCode::INTERNAL_ERROR,
                 Formatter()
                     << "Wrong curve " << elementName << " size, should be "
                     << nIgm << ", " << nIsm << ", " << nPixels << " but is "
                     << dim1 << ", " << dim2 << ", " << dim3);
        }
      }
    }
  }

  void setFlux(T3DList<Float64> const &inputFlux) {
    checkCurveElement(inputFlux, "flux");
    flux = inputFlux;
  }

  void setPixelsToUse(T3DList<bool> const &inputPixelsToUse) {
    checkCurveElement(inputPixelsToUse, "pixelsToUse");
    pixelsToUse = inputPixelsToUse;
  }

  void setFluxError(T3DList<Float64> const &inputFluxError) {
    checkCurveElement(inputFluxError, "fluxError");
    fluxError = inputFluxError;
  }

  TCurve toCurve(Int16 igmIdx, Int16 ismIdx) const {
    if (ismIdx > nIsm)
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "T3DCurve::toCurve ismIdx = " << ismIdx << " < "
                         << nIsm);
    if (igmIdx > nIgm)
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "T3DCurve::toCurve ismIdx = " << igmIdx << " < "
                         << nIgm);
    TCurve curve(lambda.size());
    curve.setLambda(lambda);
    curve.setFlux(flux[igmIdx][ismIdx]);
    curve.setPixelsToUse(pixelsToUse[igmIdx][ismIdx]);
    curve.setFluxError(fluxError[igmIdx][ismIdx]);
    return curve;
  }
};
} // namespace NSEpic
#endif