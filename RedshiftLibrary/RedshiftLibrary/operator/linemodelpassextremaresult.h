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
#ifndef _REDSHIFT_OPERATORRESULT_LINEMODEL_
#define _REDSHIFT_OPERATORRESULT_LINEMODEL_

#include <unordered_set>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/continuum/indexes.h"
#include "RedshiftLibrary/linemodel/linemodelextremaresult.h"

namespace NSEpic {

class CLineModelPassExtremaResult {

public:
  CLineModelPassExtremaResult(Int32 n);

  CLineModelPassExtremaResult() = default;

  void Resize(Int32 size);
  TInt32List
  getUniqueCandidates(const CLineModelPassExtremaResult &results_b) const;
  TFloat64List GetRedshifts() const;
  void fillWithContinuumModelSolutionAtIndex(
      Int32 i, const CContinuumModelSolution &contModelSol);

  TCandidateZbyRank m_ranked_candidates;

  std::vector<CContinuumModelSolution> m_fittedContinuum;

  std::vector<TFloat64List> ExtendedRedshifts; // z range around extrema

  // line width
  TFloat64List Elv;                    // emission line width
  TFloat64List Alv;                    // absorption line width
  std::vector<TFloat64List> GroupsELv; // per fitting group line width , EL
  std::vector<TFloat64List> GroupsALv; // per fitting group line width , AL

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
  Int32 size() const { return m_ranked_candidates.size(); }
};

} // namespace NSEpic
#endif
