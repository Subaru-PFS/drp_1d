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

#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/linemodel/continuummodelsolution.h"
#include "RedshiftLibrary/linemodel/linemodelsolution.h"
#include "RedshiftLibrary/operator/extremaresult.h"
#include "RedshiftLibrary/processflow/result.h"

namespace NSEpic {
class CModelSpectrumResult;
class CLineModelSolution;
class CModelRulesResult;
class CSpectraFluxResult;
class CLineModelFitting;

class CLineModelResult;
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
  std::vector<std::shared_ptr<const CSpectraFluxResult>>
      m_savedModelContinuumSpectrumResults;
  std::vector<std::shared_ptr<const CModelSpectrumResult>>
      m_savedModelSpectrumResults;

  CExtremaResult<TLineModelResult>(const TCandidateZbyRank &zCandidates) {
    this->m_type = "LineModelExtremaResult";
    Int32 i = 0;
    for (std::pair<std::string, const std::shared_ptr<TCandidateZ> &> cand :
         zCandidates) {
      this->m_ranked_candidates.push_back(
          std::make_pair<std::string, std::shared_ptr<TLineModelResult>>(
              std::string(cand.first),
              std::make_shared<TLineModelResult>(*cand.second)));
      i++;
    }
    this->Resize(zCandidates.size());
  };

  // mainly for saving parentObject
  CExtremaResult<TLineModelResult>(
      const std::pair<std::string, std::shared_ptr<TCandidateZ>> cand) {
    this->m_type = "LineModelExtremaResult";
    this->m_ranked_candidates.push_back(
        std::make_pair<std::string, std::shared_ptr<TLineModelResult>>(
            std::string(cand.first),
            std::make_shared<TLineModelResult>(*cand.second)));
  };

  void Resize(Int32 size) {
    // CExtremaResult::Resize(size);

    m_savedModelFittingResults.resize(size);
    m_savedModelRulesResults.resize(size);
    m_savedModelContinuumSpectrumResults.resize(size);
    m_savedModelSpectrumResults.resize(size);
  }

  void setCandidateFromContinuumSolution(int rank,
                                         CContinuumModelSolution cms) {
    m_ranked_candidates[rank].second =
        std::make_shared<TLineModelResult>(TLineModelResult(cms));
  }

  std::shared_ptr<const COperatorResult>
  getCandidate(const int &rank, const std::string &dataset,
               bool firstpassResults = false) const override;

  const std::string &
  getCandidateDatasetType(const std::string &dataset) const override;

  bool HasCandidateDataset(const std::string &dataset) const override;

  std::shared_ptr<const COperatorResult>
  getCandidateParent(const int &rank, const std::string &dataset) const;
};

typedef CExtremaResult<TLineModelResult> LineModelExtremaResult;

} // namespace NSEpic

#endif
