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
#ifndef _REDSHIFT_CURVE_
#define _REDSHIFT_CURVE_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"

namespace NSEpic {

struct TCurveElement {
  Float64 lambda;
  Float64 flux;
  Float64 fluxError;
};

struct TCurve {
  TCurve();
  TCurve(TList<Float64> lambda, TList<Float64> flux, TList<Float64> fluxError);

  TCurveElement get_at_index(Int32 idx) const;

  void push_back(TCurveElement const &elem);
  Int32 size() const;

  void setLambda(TFloat64List inputLambda);

  void setFlux(TList<Float64> inputFlux);
  void setFluxError(TList<Float64> inputFluxError);
  void sort();
  void reserve(Int32 size);

  void checkIdx(Int32 pixelIdx) const;

  const TAxisSampleList &getLambda() const { return lambda; };
  const TFloat64List &getFlux() const { return flux; };
  const TFloat64List &getFluxError() const { return fluxError; };

  Float64 getLambdaAt(Int32 pixelIdx) const;
  Float64 getFluxAt(Int32 pixelIdx) const;
  Float64 getFluxErrorAt(Int32 pixelIdx) const;

private:
  TAxisSampleList lambda;
  TFloat64List flux;
  TFloat64List fluxError;
};

} // namespace NSEpic
#endif