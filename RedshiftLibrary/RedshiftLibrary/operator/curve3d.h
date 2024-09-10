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
#ifndef _REDSHIFT_CURVE3D_
#define _REDSHIFT_CURVE3D_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/operator/curve.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"

namespace NSEpic {

struct T3DCurve {
  T3DCurve(Int32 nigm, Int32 nism);

  // rule of 5 defaults
  ~T3DCurve() = default;
  T3DCurve(const T3DCurve &) = default;
  T3DCurve(T3DCurve &&) = default;
  T3DCurve &operator=(const T3DCurve &) = default;
  T3DCurve &operator=(T3DCurve &&) = default;

  Int32 size() const { return lambda.size(); }

  void copyCurveAtAllIgmIsm(TCurve &curve);

  void setLambda(TFloat64List inputLambda);
  template <typename T>
  void checkCurveElement(T3DList<T> const &inputElement,
                         std::string elementName);

  void setFlux(T3DList<Float64> inputFlux);
  void setFluxError(T3DList<Float64> inputFluxError);
  void setIsExtincted(T3DList<bool> inputIsExtincted);
  void setIsSnrCompliant(TList<bool> inputIsSnrCompliant);
  void setMask(TList<uint8_t> inputMask);

  void setFluxAt(Int16 igmIdx, Int16 ismIdx, Int32 pixelIdx, Float64 value);
  void setFluxErrorAt(Int16 igmIdx, Int16 ismIdx, Int32 pixelIdx,
                      Float64 value);

  TCurve toCurve(Int16 igmIdx, Int16 ismIdx) const;
  TCurve toCoefCurve(Int16 igmIdx, Int16 ismIdx) const;

  bool pixelIsCoefValid(Int16 igmIdx, Int16 ismIdx, Int32 pixelIdx) const;
  bool pixelIsChi2Valid(Int32 pixelIdx) const;

  void checkIgmIdx(Int16 igmIdx) const;
  void checkIsmIdx(Int16 ismIdx) const;
  void checkPixelIdx(Int16 pixelIdx) const;
  void checkIdxs(Int16 igmIdx, Int16 ismIdx, Int32 pixelIdx) const;

  const TList<Float64> &getLambda() const { return lambda; };
  const T3DList<Float64> &getFlux() const { return flux; };
  const T3DList<Float64> &getFluxError() const { return fluxError; };
  const T3DList<bool> &getIsExtincted() const { return isExtincted; };
  const TList<bool> &getIsSnrCompliant() const { return isSnrCompliant; };
  const TList<uint8_t> &getMask() const { return mask; };

  Float64 getLambdaAt(Int32 pixelIdx) const;
  Float64 getFluxAt(Int16 igmIdx, Int16 ismIdx, Int32 pixelIdx) const;
  Float64 getFluxErrorAt(Int16 igmIdx, Int16 ismIdx, Int32 pixelIdx) const;
  bool getIsExtinctedAt(Int16 igmIdx, Int16 ismIdx, Int32 pixelIdx) const;

  Int16 getNIgm() const { return nIgm; };
  Int16 getNIsm() const { return nIsm; };

private:
  Int32 nIgm;
  Int32 nIsm;
  TList<Float64> lambda;
  T3DList<Float64> flux;
  T3DList<Float64> fluxError;
  T3DList<bool> isExtincted;
  TList<bool> isSnrCompliant;
  TList<uint8_t> mask;
};

} // namespace NSEpic
#endif