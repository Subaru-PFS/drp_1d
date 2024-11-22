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
#ifndef _REDSHIFT_STATISTICS_PDFCANDIDATESZRESULT_
#define _REDSHIFT_STATISTICS_PDFCANDIDATESZRESULT_

#include <ostream>
#include <string>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/processflow/result.h"
#include "RedshiftLibrary/statistics/pdfcandidatesz.h"

namespace NSEpic {

template <typename T>
using TRankedCandidates =
    std::vector<std::pair<std::string, std::shared_ptr<T>>>;
template <class T> class CPdfCandidateszResult : public COperatorResult {

public:
  CPdfCandidateszResult(Int32 optMethod = 0)
      : COperatorResult("PdfCandidatesZResult"), m_optMethod(optMethod){};

  // rule of 5 defaults
  CPdfCandidateszResult(const CPdfCandidateszResult &) = default;
  CPdfCandidateszResult(CPdfCandidateszResult &&) = default;
  CPdfCandidateszResult &operator=(const CPdfCandidateszResult &) = default;
  CPdfCandidateszResult &operator=(CPdfCandidateszResult &&) = default;
  virtual ~CPdfCandidateszResult() = default;
  Int32 m_optMethod; // 0: direct integration, 1:gaussian fit

  Int32 size() const { return m_ranked_candidates.size(); }

  std::shared_ptr<const T> getRankedCandidateCPtr(int rank) const {
    return std::dynamic_pointer_cast<const T>(m_ranked_candidates[rank].second);
  }

  const std::shared_ptr<T> &getRankedCandidatePtr(int rank) const {
    return m_ranked_candidates[rank].second;
  }

  TRankedCandidates<T> m_ranked_candidates;

  TCandidateZbyRank getCandidatesZByRank() {
    TCandidateZbyRank ret;
    for (auto &cand : m_ranked_candidates) {
      ret.push_back(std::make_pair(
          cand.first, std::dynamic_pointer_cast<TCandidateZ>(cand.second)));
    }
    return ret;
  }

  TStringList GetIDs() const {
    TStringList ids;
    ids.reserve(m_ranked_candidates.size());
    for (auto c : m_ranked_candidates)
      ids.push_back(c.first);
    return ids;
  }

  std::string ID(Int32 i) const { return m_ranked_candidates[i].first; }

  Float64 Redshift(Int32 i) const {
    return m_ranked_candidates[i].second->Redshift;
  }
  Float64 ValProba(Int32 i) const {
    return m_ranked_candidates[i].second->ValProba;
  }
  Float64 ValSumProba(Int32 i) const {
    return m_ranked_candidates[i].second->ValSumProba;
  }
  Float64 DeltaZ(Int32 i) const {
    return m_ranked_candidates[i].second->Deltaz;
  }
};

typedef CPdfCandidateszResult<TCandidateZ> PdfCandidatesZResult;
} // namespace NSEpic

#endif
