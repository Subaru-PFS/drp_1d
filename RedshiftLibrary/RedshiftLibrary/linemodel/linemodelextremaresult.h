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
#ifndef _REDSHIFT_LINEMODEL_LINEMODELEXTREMARESULT_
#define _REDSHIFT_LINEMODEL_LINEMODELEXTREMARESULT_

#include <unordered_set>

#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/linemodel/continuummodelsolution.h"
#include "RedshiftLibrary/linemodel/linemodelsolution.h"
#include "RedshiftLibrary/operator/extremaresult.h"
#include "RedshiftLibrary/operator/modelphotvalueresult.h"
#include "RedshiftLibrary/processflow/result.h"

namespace NSEpic {
class CModelSpectrumResult;
class CLineModelSolution;
class CModelRulesResult;
class CModelPhotValueResult;
class CLineModelFitting;
class CTplratioManager;
class CLineModelResult;
class TCandidateZ;
#include "RedshiftLibrary/linemodel/linemodelextremaresult.i"

template <>
class CExtremaResult<TLineModelResult>
    : public CPdfCandidateszResult<TLineModelResult> {
public:
  std::vector<TFloat64List> ExtendedRedshifts; // z range around extrema

  std::vector<std::shared_ptr<const CLineModelSolution>>
      m_savedModelFittingResults;
  std::vector<std::shared_ptr<const CModelRulesResult>>
      m_savedModelRulesResults;
  std::vector<std::shared_ptr<const CModelSpectrumResult>>
      m_savedModelContinuumSpectrumResults;
  std::vector<std::shared_ptr<const CModelSpectrumResult>>
      m_savedModelSpectrumResults;
  std::vector<std::shared_ptr<const CModelPhotValueResult>> m_modelPhotValues;

  CExtremaResult<TLineModelResult>(const TCandidateZbyRank &zCandidates)
      : m_savedModelFittingResults(zCandidates.size()),
        m_savedModelRulesResults(zCandidates.size()),
        m_savedModelContinuumSpectrumResults(zCandidates.size()),
        m_savedModelSpectrumResults(zCandidates.size()),
        m_modelPhotValues(zCandidates.size()) {
    m_type = "LineModelExtremaResult";
    SetRankedCandidates(zCandidates);
  };

  CExtremaResult<TLineModelResult>() { m_type = "LineModelExtremaResult"; }
  void SetRankedCandidates(const TCandidateZbyRank &zCandidates) {
    m_ranked_candidates.clear();
    for (const auto &cand : zCandidates) {
      m_ranked_candidates.push_back(
          std::make_pair(std::string(cand.first),
                         std::make_shared<TLineModelResult>(*cand.second)));
    }
  }

  std::shared_ptr<const COperatorResult>
  getCandidate(const int &rank, const std::string &dataset,
               bool firstpassResults = false) const override;

  const std::string &
  getCandidateDatasetType(const std::string &dataset) const override;

  bool HasCandidateDataset(const std::string &dataset) const override;

  std::shared_ptr<const COperatorResult>
  getCandidateParent(const int &rank, const std::string &dataset) const;

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

typedef CExtremaResult<TLineModelResult> LineModelExtremaResult;

} // namespace NSEpic

#endif
