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
#include "RedshiftLibrary/operator/operator.h"
#include "RedshiftLibrary/operator/passextremaresult.h"
#include "RedshiftLibrary/operator/twopassresult.h"

namespace TwoPass {
class spanRedshift_test;
} // namespace TwoPass

namespace NSEpic {
template <class T> class COperatorTwoPass : public COperator {
public:
  void init(const Float64 halfWindowSize, const bool zLogSampling,
            const TFloat64List &redshifts, const Float64 fineStep);
  TFloat64List spanRedshiftWindow(const Float64 z) const;
  void buildExtendedRedshifts(CPassExtremaResult &passExtremaResult);
  void updateRedshiftGridAndResults(std::shared_ptr<CTwoPassResult> result);
  TZGridListParams getSPZGridParams();

  CPassExtremaResult &getFirstPassExtremaResult() {
    return m_firstpass_extremaResult;
  }

protected:
  friend class TwoPass::spanRedshift_test;
  Float64 m_secondPass_halfwindowsize;
  bool m_zLogSampling;
  TFloat64List m_Redshifts;
  Float64 m_fineStep;
  std::vector<TFloat64List> m_extendedRedshifts; // z range around extrema
  CPassExtremaResult m_firstpass_extremaResult;
  std::shared_ptr<CExtremaResult<T>> m_rebaseArg;
};

template <class T>
void COperatorTwoPass<T>::init(
    const Float64 halfWindowSize, const bool zLogSampling,
    const TFloat64List &redshifts, // TODO check if should be a shared_ptr /
                                   // check that is always synchro with operator
    const Float64 fineStep) {
  m_secondPass_halfwindowsize = halfWindowSize;
  m_zLogSampling = zLogSampling;
  m_Redshifts = redshifts;
  m_fineStep = fineStep;
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
  windowRange.IntersectWith(m_Redshifts);
  CZGridParam zparam(windowRange, m_fineStep, z);

  return zparam.getZGrid(m_zLogSampling);
};

template <class T>
void COperatorTwoPass<T>::buildExtendedRedshifts(
    CPassExtremaResult &passExtremaResult) {
  // Refine redshift grid around extrema results redshifts
  Int32 nExtremaResults = passExtremaResult.size();
  m_extendedRedshifts.reserve(nExtremaResults);

  for (Int32 candidateIdx = 0; candidateIdx < nExtremaResults; candidateIdx++) {
    const std::shared_ptr<const TCandidateZ> &cand =
        passExtremaResult.m_ranked_candidates[candidateIdx].second;

    Log.LogInfo(Formatter() << "  Operator-TwoPass: Raw extr #" << candidateIdx
                            << ", z_e.X=" << cand->Redshift
                            << ", m_e.Y=" << cand->ValProba);
    m_extendedRedshifts.push_back(spanRedshiftWindow(cand->Redshift));
  }
}

template <class T>
void COperatorTwoPass<T>::updateRedshiftGridAndResults(
    std::shared_ptr<CTwoPassResult> result) {

  for (auto &subgrid : m_extendedRedshifts) {
    Int32 imin, ndup;
    std::tie(imin, ndup) =
        CZGridListParams::insertSubgrid(subgrid, m_Redshifts);
    result->updateVectors(imin, ndup, subgrid.size());
  }
  result->Redshifts = m_Redshifts;
}

// only for secondpass grid
template <class T> TZGridListParams COperatorTwoPass<T>::getSPZGridParams() {
  Int32 s = m_extendedRedshifts.size();
  TZGridListParams centeredZgrid_params(s);
  for (Int32 i = 0; i < s; i++) {
    const auto &extendedGrid = m_extendedRedshifts[i];
    centeredZgrid_params[i] = CZGridParam(
        TFloat64Range(extendedGrid), m_fineStep,
        m_firstpass_extremaResult.m_ranked_candidates[i].second->Redshift);
  }
  return centeredZgrid_params;
}

} // namespace NSEpic

#endif
