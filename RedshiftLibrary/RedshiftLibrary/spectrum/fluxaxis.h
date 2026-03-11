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
  using CSpectrumAxis::CSpectrumAxis;
  explicit CSpectrumFluxAxis(CSpectrumAxis otherFlux)
      : CSpectrumAxis(std::move(otherFlux)) {}
  CSpectrumFluxAxis(CSpectrumAxis otherFlux, CSpectrumNoiseAxis otherError);
  CSpectrumFluxAxis(const Float64 *samples, Int32 n, const Float64 *error,
                    const Int32 m);

  bool HasError() const;
  const CSpectrumNoiseAxis &GetError() const;
  Float64 GetWeight(Int32 idx, Float64 normFactor = 1.0) const;
  Float64 GetInverseWeight(Int32 Idx, Float64 normFactor = 1.0) const;
  std::pair<TAxisSampleList, TAxisSampleList> GetSamplesAndErrorVector() &&;
  void setSamplesVector(TAxisSampleList axisList) override;
  void setError(CSpectrumNoiseAxis otherError);
  void resize(Int32 s, Float64 valudDef = 0.0) override;
  void clear() override;
  bool ApplyMeanSmooth(Int32 kernelHalfWidth);
  bool ApplyMedianSmooth(Int32 kernelHalfWidth);
  Float64 computeMaxAbsValue(Int32 imin, Int32 imax) const;
  bool ComputeMeanAndSDev(const CMask &mask, Float64 &mean,
                          Float64 &sdev) const;
  Float64 ComputeRMSDiff(const CSpectrumFluxAxis &other);
  const TBoolList checkFlux() const;
  bool correctFluxAndNoiseAxis(Int32 iMin, Int32 iMax, Float64 coeffCorr);
  CSpectrumFluxAxis &operator*=(Float64 op) override;
  CSpectrumFluxAxis &operator/=(Float64 op) override;
  using CSpectrumAxis::operator+=, CSpectrumAxis::operator-=;
  CSpectrumFluxAxis &operator+=(CSpectrumFluxAxis const &other);
  CSpectrumFluxAxis &operator-=(CSpectrumFluxAxis const &);
  // Hidden friend symmetric operators (ie non-member, but here for ADL )
  friend CSpectrumFluxAxis operator*(CSpectrumFluxAxis axis, Float64 op);
  friend CSpectrumFluxAxis operator*(Float64 op, CSpectrumFluxAxis axis);
  friend CSpectrumFluxAxis operator+(CSpectrumFluxAxis const &axis1,
                                     CSpectrumFluxAxis const &axis2);
  friend CSpectrumFluxAxis operator+(CSpectrumAxis const &axis1,
                                     CSpectrumFluxAxis const &axis2);
  friend CSpectrumFluxAxis operator+(CSpectrumFluxAxis const &axis1,
                                     CSpectrumAxis const &axis2);
  friend CSpectrumFluxAxis operator-(CSpectrumFluxAxis const &axis1,
                                     CSpectrumFluxAxis const &axis2);
  friend CSpectrumFluxAxis operator-(CSpectrumAxis const &axis1,
                                     CSpectrumFluxAxis const &axis2);
  friend CSpectrumFluxAxis operator-(CSpectrumFluxAxis const &axis1,
                                     CSpectrumAxis const &axis2);
  // member assymetric operator
  CSpectrumFluxAxis operator/(Float64 op) const;

  CSpectrumFluxAxis extract(Int32 startIdx, Int32 endIdx) const;
  CSpectrumFluxAxis MaskAxis(const TMaskList &) const;
  void Invert() = delete;

private:
  friend class FluxAxis_test::ComputeMeanAndSDev_test;

  void checkSizes() const;

  CSpectrumNoiseAxis m_StdError; // STD
  bool m_hasStdError = false;
};

inline bool CSpectrumFluxAxis::HasError() const { return m_hasStdError; }

inline const CSpectrumNoiseAxis &CSpectrumFluxAxis::GetError() const {
  if (!m_hasStdError)
    THROWG(ErrorCode::INTERNAL_ERROR, "spectrum flux axis has no error vector");
  return m_StdError;
}

inline Float64 CSpectrumFluxAxis::GetWeight(Int32 idx,
                                            Float64 normFactor) const {
  if (m_hasStdError) {
    auto const &err = m_StdError[idx] * normFactor;
    return 1 / (err * err);
  }
  return 1 / (normFactor * normFactor);
}

inline Float64 CSpectrumFluxAxis::GetInverseWeight(Int32 idx,
                                                   Float64 normFactor) const {
  if (m_hasStdError) {
    auto const &err = m_StdError[idx] * normFactor;
    return err * err;
  }
  return 1;
}

inline std::pair<TAxisSampleList, TAxisSampleList>
CSpectrumFluxAxis::GetSamplesAndErrorVector() && {
  return {std::move(m_Samples), std::move(m_StdError).GetSamplesVector()};
}

inline CSpectrumFluxAxis CSpectrumFluxAxis::extract(Int32 startIdx,
                                                    Int32 endIdx) const {
  if (HasError())
    return CSpectrumFluxAxis(CSpectrumAxis::extract(startIdx, endIdx),
                             m_StdError.extract(startIdx, endIdx));
  else
    return CSpectrumFluxAxis(CSpectrumAxis::extract(startIdx, endIdx));
}

inline CSpectrumFluxAxis &CSpectrumFluxAxis::operator*=(Float64 op) {
  CSpectrumAxis::operator*=(op);
  if (m_hasStdError)
    m_StdError *= op;
  return *this;
}

inline CSpectrumFluxAxis &CSpectrumFluxAxis::operator/=(Float64 op) {
  operator*=(1 / op);
  return *this;
}

inline CSpectrumFluxAxis operator*(CSpectrumFluxAxis axis, Float64 op) {
  CSpectrumFluxAxis multipliedAxis(std::move(axis));
  multipliedAxis *= op;
  return multipliedAxis;
}

inline CSpectrumFluxAxis operator*(Float64 op, CSpectrumFluxAxis axis) {
  return std::move(axis) * op;
}

inline CSpectrumFluxAxis operator+(CSpectrumFluxAxis const &axis1,
                                   CSpectrumFluxAxis const &axis2) {
  CSpectrumFluxAxis sumaxis(axis1);
  sumaxis += axis2;
  return sumaxis;
}

inline CSpectrumFluxAxis operator+(CSpectrumAxis const &axis1,
                                   CSpectrumFluxAxis const &axis2) {
  CSpectrumFluxAxis sumaxis(axis1);
  sumaxis += axis2;
  return sumaxis;
}

inline CSpectrumFluxAxis operator+(CSpectrumFluxAxis const &axis1,
                                   CSpectrumAxis const &axis2) {
  return axis2 + axis1;
}

inline CSpectrumFluxAxis operator-(CSpectrumFluxAxis const &axis1,
                                   CSpectrumFluxAxis const &axis2) {
  CSpectrumFluxAxis diffaxis(axis1);
  diffaxis -= axis2;
  return diffaxis;
}

inline CSpectrumFluxAxis operator-(CSpectrumAxis const &axis1,
                                   CSpectrumFluxAxis const &axis2) {
  CSpectrumFluxAxis diffaxis(axis1);
  diffaxis -= axis2;
  return diffaxis;
}

inline CSpectrumFluxAxis operator-(CSpectrumFluxAxis const &axis1,
                                   CSpectrumAxis const &axis2) {
  return axis2 - axis1;
}

inline CSpectrumFluxAxis CSpectrumFluxAxis::operator/(Float64 op) const {
  CSpectrumFluxAxis dividedAxis(*this);
  dividedAxis /= op;
  return dividedAxis;
}

} // namespace NSEpic

#endif
