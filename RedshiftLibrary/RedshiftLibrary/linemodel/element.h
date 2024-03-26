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
#ifndef _REDSHIFT_LINEMODEL_ELEMENT_
#define _REDSHIFT_LINEMODEL_ELEMENT_

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/polynom.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/line/lineprofile.h"
#include "RedshiftLibrary/linemodel/elementparam.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

class RuleStrongHigherThanWeak_fixture;

namespace RuleStrongHigherThanWeak_test {
class Correct_test_no_change;
class Correct_test_one_high_weak;
} // namespace RuleStrongHigherThanWeak_test

namespace NSEpic {

using TLineModelElementParam_ptr = std::shared_ptr<TLineModelElementParam>;

/**
 * \ingroup Redshift
 */
class CLineModelElement {

public:
  CLineModelElement(const TLineModelElementParam_ptr elementParam);

  Float64 GetObservedPosition(Int32 line_index, Float64 redshift,
                              bool doAsymfitdelta = true) const;
  Float64 GetLineProfileAtRedshift(Int32 line_index, Float64 redshift,
                                   Float64 x) const;
  void getObservedPositionAndLineWidth(Int32 line_index, Float64 redshift,
                                       Float64 &mu, Float64 &sigma,
                                       bool doAsymfitdelta = true) const;
  void prepareSupport(const CSpectrumSpectralAxis &spectralAxis,
                      Float64 redshift, const TFloat64Range &lambdaRange,
                      Float64 max_offset = 0.0);
  TInt32RangeList getSupport() const;
  TInt32RangeList getTheoreticalSupport() const;
  void EstimateTheoreticalSupport(Int32 line_index,
                                  const CSpectrumSpectralAxis &spectralAxis,
                                  Float64 redshift,
                                  const TFloat64Range &lambdaRange,
                                  Float64 max_offset = 0.0);
  void computeOutsideLambdaRange();

  TInt32Range getSupportSubElt(Int32 line_index) const;
  TInt32Range getTheoreticalSupportSubElt(Int32 line_id) const;

  static TInt32Range
  EstimateIndexRange(const CSpectrumSpectralAxis &spectralAxis, Float64 mu,
                     const TFloat64Range &lambdaRange, Float64 winsizeAngstrom);

  Float64 GetContinuumAtCenterProfile(
      Int32 line_index, const CSpectrumSpectralAxis &spectralAxis,
      Float64 redshift, const CSpectrumFluxAxis &continuumfluxAxis,
      bool enableAmplitudeOffsets = false) const;

  Float64 getModelAtLambda(Float64 lambda, Float64 redshift,
                           Float64 continuumFlux,
                           Int32 line_index = undefIdx) const;
  Float64 GetModelDerivAmplitudeAtLambda(Float64 lambda, Float64 redshift,
                                         Float64 continuumFlux) const;
  Float64 GetModelDerivVelAtLambda(Float64 lambda, Float64 redshift,
                                   Float64 continuumFlux) const;
  Float64 GetModelDerivContinuumAmpAtLambda(Float64 lambda, Float64 redshift,
                                            Float64 continuumFluxUnscale) const;
  Float64 GetModelDerivZAtLambda(Float64 lambda, Float64 redshift,
                                 Float64 continuumFlux,
                                 Float64 continuumFluxDerivZ = 0.0) const;

  void addToSpectrumModel(const CSpectrumSpectralAxis &modelspectralAxis,
                          CSpectrumFluxAxis &modelfluxAxis,
                          const CSpectrumFluxAxis &continuumfluxAxis,
                          Float64 redshift, Int32 line_index = undefIdx) const;
  void
  addToSpectrumModelDerivVel(const CSpectrumSpectralAxis &modelspectralAxis,
                             CSpectrumFluxAxis &modelfluxAxis,
                             const CSpectrumFluxAxis &continuumFluxAxis,
                             Float64 redshift, bool emissionLine) const;

  void initSpectrumModel(CSpectrumFluxAxis &modelfluxAxis,
                         const CSpectrumFluxAxis &continuumfluxAxis,
                         Int32 line_index = undefIdx) const;
  void initSpectrumModelPolynomial(CSpectrumFluxAxis &modelfluxAxis,
                                   const CSpectrumSpectralAxis &spcAxis,
                                   Int32 line_index) const;

  void SetLSF(const std::shared_ptr<const CLSF> &lsf);

  Int32 GetSize() const { return m_size; };

  bool IsOutsideLambdaRange() const;
  Float64 GetLineWidth(Float64 lambda) const;
  bool IsOutsideLambdaRangeLine(Int32 line_index) const;
  void SetOutsideLambdaRangeList(Int32 line_index);

  void SetLineProfile(Int32 line_index, CLineProfile_ptr &&profile);

  bool isLineActiveOnSupport(Int32 line_indexA, Int32 line_indexB) const;
  Int32 getStartNoOverlap(Int32 line_index) const;
  Int32 getEndNoOverlap(Int32 line_index) const;

  void debug(std::ostream &os) const;
  void dumpElement(std::ostream &os) const;

  TLineModelElementParam_ptr getElementParam() { return m_ElementParam; }
  const TLineModelElementParam_ptr &getElementParam() const {
    return m_ElementParam;
  }

  Float64 GetElementAmplitude() const;
  Float64 GetElementError() const;
  void SetFittedAmplitude(Int32 line_index, Float64 A, Float64 AStd);

  Int32 computeCrossProducts(Float64 redshift,
                             const CSpectrumSpectralAxis &spectralAxis,
                             const CSpectrumFluxAxis &noContinuumfluxAxis,
                             const CSpectrumFluxAxis &continuumfluxAxis,
                             Int32 line_index);

  /*  void fitAmplitude(Float64 redshift, const CSpectrumSpectralAxis
     &spectralAxis, const CSpectrumFluxAxis &noContinuumfluxAxis, const
     CSpectrumFluxAxis &continuumfluxAxis, Int32 line_index = undefIdx);
  */
protected:
  friend ::RuleStrongHigherThanWeak_fixture;
  friend RuleStrongHigherThanWeak_test::Correct_test_no_change;
  friend RuleStrongHigherThanWeak_test::Correct_test_one_high_weak;

  const TLineModelElementParam_ptr m_ElementParam;

  const Float64 m_OutsideLambdaRangeOverlapThreshold;
  bool m_OutsideLambdaRange;

  std::shared_ptr<const CLSF> m_LSF;

  std::vector<TBoolList> m_LineIsActiveOnSupport;

  TInt32List m_StartNoOverlap;
  TInt32List m_EndNoOverlap;
  TInt32List m_StartTheoretical;
  TInt32List m_EndTheoretical;

  TBoolList m_OutsideLambdaRangeList;
  Int32 m_size;
};

inline bool CLineModelElement::IsOutsideLambdaRange() const {
  return m_OutsideLambdaRange;
}
/**
 * \brief Returns whether the line with id line_id is outside the lambda
 *range.
 **/
inline bool
CLineModelElement::IsOutsideLambdaRangeLine(Int32 line_index) const {
  return m_OutsideLambdaRangeList[line_index];
}

inline bool CLineModelElement::isLineActiveOnSupport(Int32 lineindexA,
                                                     Int32 lineindexB) const {
  return m_LineIsActiveOnSupport[lineindexA][lineindexB];
}

inline Int32 CLineModelElement::getStartNoOverlap(Int32 line_index) const {
  return m_StartNoOverlap[line_index];
}

inline Int32 CLineModelElement::getEndNoOverlap(Int32 line_index) const {
  return m_EndNoOverlap[line_index];
}

inline void CLineModelElement::SetLSF(const std::shared_ptr<const CLSF> &lsf) {
  m_LSF = lsf;
}

inline void CLineModelElement::SetOutsideLambdaRangeList(Int32 line_index) {
  m_OutsideLambdaRangeList[line_index] = true;
}

} // namespace NSEpic

#endif // _REDSHIFT_LINEMODEL_ELEMENT_
