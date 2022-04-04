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

#include "RedshiftLibrary/operator/linemodelpassextremaresult.h"
#include <cfloat>
using namespace NSEpic;
using namespace std;

CLineModelPassExtremaResult::CLineModelPassExtremaResult(Int32 n) { Resize(n); }

void CLineModelPassExtremaResult::Resize(Int32 size) {
  m_ranked_candidates.resize(size);
  FittedTplName.resize(size);
  FittedTplAmplitude.resize(size);
  FittedTplAmplitudeError.resize(size);
  FittedTplMerit.resize(size);
  FittedTplMeritPhot.resize(size);
  FittedTplEbmvCoeff.resize(size);
  FittedTplMeiksinIdx.resize(size);
  FittedTplDtm.resize(size);
  FittedTplMtm.resize(size);
  FittedTplLogPrior.resize(size);
  FittedTplSNR.resize(size);

  MeritContinuum.resize(size);

  mTransposeM.resize(size);
  CorrScaleMarg.resize(size);
  NDof.resize(size);
  Redshift_lmfit.resize(size);
  snrHa.resize(size);
  lfHa.resize(size);
  snrOII.resize(size);
  lfOII.resize(size);

  ExtendedRedshifts.resize(size);
  NLinesOverThreshold.resize(size);
  LogArea.resize(size);
  LogAreaCorrectedExtrema.resize(size);
  SigmaZ.resize(size);

  StrongELSNR.resize(size);
  StrongELSNRAboveCut.resize(size);
  bic.resize(size);
  ContinuumIndexes.resize(size);
  OutsideLinesMask.resize(size);
  OutsideLinesSTDFlux.resize(size);
  OutsideLinesSTDError.resize(size);

  Elv.resize(size);
  Alv.resize(size);
  GroupsELv.resize(size);
  GroupsALv.resize(size);
  for (Int32 ke = 0; ke < size; ke++) {
    GroupsELv[ke] = TFloat64List(250, -1); // WARNING: hardcoded ngroups max
    GroupsALv[ke] = TFloat64List(250, -1); // WARNING: hardcoded ngroups max
  }

  FittedTplRedshift.resize(size);
  FittedTplpCoeffs.resize(size);

  FittedTplratioName.resize(size);
  FittedTplratioAmplitude.resize(size);
  FittedTplratioDtm.resize(size);
  FittedTplratioMtm.resize(size);
  FittedTplratioIsmCoeff.resize(size);
}

TFloat64List CLineModelPassExtremaResult::GetRedshifts() const {
  TFloat64List redshifts;
  redshifts.reserve(m_ranked_candidates.size());
  for (auto c : m_ranked_candidates)
    redshifts.push_back(c.second.Redshift);
  return redshifts;
}

TInt32List CLineModelPassExtremaResult::getUniqueCandidates(
    std::shared_ptr<const CLineModelPassExtremaResult> results_b) {
  TInt32List uniqueIndices;
  Float64 skip_thres_absdiffz =
      5e-4; // threshold to remove duplicate extrema/candidates
  for (Int32 keb = 0; keb < results_b->m_ranked_candidates.size(); keb++) {
    const Float64 &z_fpb = results_b->m_ranked_candidates[keb].second.Redshift;
    // skip if z_fpb is nearly the same as any z_fp
    Float64 minAbsDiffz = DBL_MAX;
    bool duplicate = false;
    for (Int32 ke = 0; ke < m_ranked_candidates.size(); ke++) {
      Float64 z_diff = z_fpb - m_ranked_candidates[ke].second.Redshift;
      if (std::abs(z_diff) < skip_thres_absdiffz) {
        duplicate = true;
        break; // no need to look for more
      }
    }
    if (!duplicate)
      uniqueIndices.push_back(keb);
  }
  return uniqueIndices;
}

TStringList CLineModelPassExtremaResult::GetIDs() const {
  TStringList ids;
  ids.reserve(m_ranked_candidates.size());
  for (auto c : m_ranked_candidates)
    ids.push_back(c.first);
  return ids;
}
