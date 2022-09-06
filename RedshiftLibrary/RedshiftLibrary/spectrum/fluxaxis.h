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
#ifndef _REDSHIFT_SPECTRUM_FLUXAXIS_
#define _REDSHIFT_SPECTRUM_FLUXAXIS_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/spectrum/noiseaxis.h"

namespace FluxAxis_test { // boost_test_suite
// all boost_auto_test_case that use private method
class ComputeMeanAndSDev_test;
} // namespace FluxAxis_test
namespace NSEpic {

class CMask;

/**
 * \ingroup Redshift
 */
class CSpectrumFluxAxis : public CSpectrumAxis {

public:
  CSpectrumFluxAxis() = default;
  // value =0. is a default value for flux and not for error.
  explicit CSpectrumFluxAxis(Int32 n, Float64 value = 0.0);
  CSpectrumFluxAxis(CSpectrumAxis otherFlux, CSpectrumNoiseAxis otherError);
  CSpectrumFluxAxis(const Float64 *samples, Int32 n);
  CSpectrumFluxAxis(const TFloat64List &samples);
  CSpectrumFluxAxis(TFloat64List &&samples);
  CSpectrumFluxAxis(const Float64 *samples, Int32 n, const Float64 *error,
                    const Int32 m);

  const CSpectrumNoiseAxis &GetError() const;
  CSpectrumNoiseAxis &GetError();

  void setError(const CSpectrumNoiseAxis &otherError);
  void SetSize(Int32 s);
  void clear();
  bool ApplyMeanSmooth(Int32 kernelHalfWidth);
  bool ApplyMedianSmooth(Int32 kernelHalfWidth);

  bool ComputeMeanAndSDev(const CMask &mask, Float64 &mean,
                          Float64 &sdev) const;
  Float64 ComputeRMSDiff(const CSpectrumFluxAxis &other);
  bool Subtract(const CSpectrumFluxAxis &other);
  bool Invert();
  CSpectrumFluxAxis
  extract(Int32 startIdx,
          Int32 endIdx) const; // this is mainly applied on m_StdError

private:
  friend class FluxAxis_test::ComputeMeanAndSDev_test;

  bool ComputeMeanAndSDevWithoutError(const CMask &mask, Float64 &mean,
                                      Float64 &sdev) const;
  bool ComputeMeanAndSDevWithError(const CMask &mask, Float64 &mean,
                                   Float64 &sdev) const;

  CSpectrumNoiseAxis m_StdError; // STD
};

inline CSpectrumNoiseAxis &CSpectrumFluxAxis::GetError() { return m_StdError; }

inline const CSpectrumNoiseAxis &CSpectrumFluxAxis::GetError() const {
  return m_StdError;
}

inline CSpectrumFluxAxis CSpectrumFluxAxis::extract(Int32 startIdx,
                                                    Int32 endIdx) const {
  return CSpectrumFluxAxis(CSpectrumAxis::extract(startIdx, endIdx),
                           m_StdError.extract(startIdx, endIdx));
}
} // namespace NSEpic

#endif
