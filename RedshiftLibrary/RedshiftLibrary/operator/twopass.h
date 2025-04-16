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
#ifndef _REDSHIFT_OPERATOR_TWOPASS_
#define _REDSHIFT_OPERATOR_TWOPASS_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/zgridparam.h"
#include "RedshiftLibrary/operator/extremaresult.h"
#include "RedshiftLibrary/operator/operator.h"
#include "RedshiftLibrary/operator/pass.h"
#include "RedshiftLibrary/operator/twopassresult.h"
#include "RedshiftLibrary/statistics/pdfcandidatesz.h"

namespace TwoPass {
class spanRedshift_test;
} // namespace TwoPass

namespace NSEpic {
struct TsecondPassIndices {
  Int32 insertionIdx = undefIdx;
  TInt32List overwrittenSourceIndices;
  std::size_t size = undefIdx;
  Int32 largeStepFactor;
};

template <class T> class COperatorTwoPass : public COperatorPass {
public:
  COperatorTwoPass() = default;
  virtual ~COperatorTwoPass() = default;
  COperatorTwoPass(const COperatorTwoPass &other) = default;
  COperatorTwoPass &operator=(const COperatorTwoPass &other) = default;
  COperatorTwoPass(COperatorTwoPass &&other) noexcept = default;
  COperatorTwoPass &operator=(COperatorTwoPass &&other) noexcept = default;

  void setTwoPassParameters(const Float64 halfWindowSize,
                            const bool zLogSampling, const Float64 fineStep,
                            Int32 twoPassZStepFactor);
  TFloat64List spanRedshiftWindow(const Float64 z) const;
  void buildExtendedRedshifts();
  TList<TsecondPassIndices> updateRedshiftGrid();
  void updateVectors(std::shared_ptr<CTwoPassResult> const &result,
                     TList<TsecondPassIndices> const &insertionIndicesList);
  void
  updateRedshiftGridAndResults(std::shared_ptr<CTwoPassResult> const &result);
  TZGridListParams getSPZGridParams();
  TCandidateZbyRank getFirstPassCandidatesZByRank() const {
    return m_firstpass_extremaResult->getCandidatesZByRank();
  }
  std::shared_ptr<const CExtremaResult<T>> getFirstPassExtremaResults() const {
    return m_firstpass_extremaResult;
  }
  std::vector<TFloat64List> getExtendedRedshifts() const {
    return m_extendedRedshifts;
  }
  TInt32Range getzIdxRangeToCompute(Int32 candidateIdx);

  void
  SetFirstPassExtremaResults(std::shared_ptr<ExtremaResult> const &result) {
    m_firstpass_extremaResult = result;
  };

protected:
  friend class TwoPass::spanRedshift_test;
  Float64 m_secondPass_halfwindowsize = NAN;
  bool m_zLogSampling = false;
  Float64 m_fineStep = NAN;
  Int32 m_twoPassZStepFactor = undefIdx;
  std::vector<TFloat64List> m_extendedRedshifts; // z range around extrema
  std::shared_ptr<CExtremaResult<T>> m_firstpass_extremaResult;
};

template <class T>
void COperatorTwoPass<T>::setTwoPassParameters(const Float64 halfWindowSize,
                                               const bool zLogSampling,
                                               const Float64 fineStep,
                                               Int32 twoPassZStepFactor) {
  m_secondPass_halfwindowsize = halfWindowSize;
  m_zLogSampling = zLogSampling;
  m_fineStep = fineStep;
  m_twoPassZStepFactor = twoPassZStepFactor;
}

template <class T>
TFloat64List COperatorTwoPass<T>::spanRedshiftWindow(const Float64 z) const {
  Float64 half_r = m_secondPass_halfwindowsize;
  Float64 half_l = m_secondPass_halfwindowsize;

  if (m_zLogSampling) {
    half_r = (exp(half_r) - 1.0) * (1. + z);
    half_l = (1.0 - exp(-half_l)) * (1. + z);
  }

  //
  TFloat64Range windowRange(z - half_l, z + half_r);
  windowRange.IntersectWith(m_redshifts);
  CZGridParam zparam(windowRange, m_fineStep, z);

  return zparam.getZGrid(m_zLogSampling);
};

template <class T> void COperatorTwoPass<T>::buildExtendedRedshifts() {
  // Refine redshift grid around extrema results redshifts

  m_extendedRedshifts.reserve(m_firstpass_extremaResult->size());
  for (const auto &[_, cand] : m_firstpass_extremaResult->m_ranked_candidates) {
    Log.LogInfo(Formatter() << "  Operator-TwoPass: Raw extr #" << cand->Rank
                            << ", z_e.X=" << cand->Redshift
                            << ", m_e.Y=" << cand->ValProba);
    m_extendedRedshifts.push_back(spanRedshiftWindow(cand->Redshift));
  }
}

template <class T>
void COperatorTwoPass<T>::updateRedshiftGridAndResults(
    std::shared_ptr<CTwoPassResult> const &result) {
  auto const insertionIndicesList = updateRedshiftGrid();
  updateVectors(result, insertionIndicesList);
}

template <class T>
TList<TsecondPassIndices> COperatorTwoPass<T>::updateRedshiftGrid() {
  TList<TsecondPassIndices> insertionIndicesList;
  for (auto &subgrid : m_extendedRedshifts) {
    auto const &[insertionIdx, overwrittenSourceIndices] =
        CZGridListParams::insertSubgrid(subgrid, m_redshifts);

    insertionIndicesList.push_back({insertionIdx,
                                    std::move(overwrittenSourceIndices),
                                    subgrid.size(), m_twoPassZStepFactor});
  }
  return insertionIndicesList;
}

template <class T>
void COperatorTwoPass<T>::updateVectors(
    std::shared_ptr<CTwoPassResult> const &result,
    TList<TsecondPassIndices> const &insertionIndicesList) {
  for (auto const &insertionIndices : insertionIndicesList)
    result->updateVectors(insertionIndices);

  result->Redshifts = m_redshifts;
}

// only for secondpass grid
template <class T> TZGridListParams COperatorTwoPass<T>::getSPZGridParams() {
  Int32 s = m_extendedRedshifts.size();
  TZGridListParams centeredZgrid_params(s);
  for (Int32 i = 0; i < s; i++) {
    const auto &extendedGrid = m_extendedRedshifts[i];
    centeredZgrid_params[i] =
        CZGridParam(TFloat64Range(extendedGrid), m_fineStep,
                    m_firstpass_extremaResult->Redshift(i));
  }
  return centeredZgrid_params;
}

template <class T>
TInt32Range COperatorTwoPass<T>::getzIdxRangeToCompute(Int32 candidateIdx) {
  auto const &extendedRedshifts = m_extendedRedshifts[candidateIdx];
  Int32 imin = undefIdx;
  Int32 imax = undefIdx;
  TFloat64Range(extendedRedshifts)
      .getClosedIntervalIndices(m_redshifts, imin, imax);
  return TInt32Range(imin, imax);
}

} // namespace NSEpic
#endif
