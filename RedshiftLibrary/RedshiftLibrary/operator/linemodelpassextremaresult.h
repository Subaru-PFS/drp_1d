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

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/continuum/indexes.h"
#include "RedshiftLibrary/linemodel/linemodelextremaresult.h"

#include <unordered_set>

namespace NSEpic {

class CLineModelPassExtremaResult {

public:
  CLineModelPassExtremaResult(Int32 n);

  CLineModelPassExtremaResult() = default;

  void Resize(Int32 size);
  TInt32List getUniqueCandidates(
      std::shared_ptr<const CLineModelPassExtremaResult> results_b);
  TFloat64List GetRedshifts() const;

  Int32 m_optMethod; // 0: direct integration, 1:gaussian fit

  TCandidateZbyRank m_ranked_candidates;

  // template continuum
  TStringList FittedTplName; // Name of the best template fitted for continuum
  TFloat64List FittedTplAmplitude; // Amplitude for the best template fitted for
                                   // continuum
  TFloat64List FittedTplAmplitudeError; // Amplitude error for the best template
                                        // fitted for continuum
  TFloat64List
      FittedTplMerit; // Chisquare for the best template fitted for continuum
  TFloat64List FittedTplMeritPhot; // extra chisquare term due to photometry
  TFloat64List FittedTplEbmvCoeff; // Calzetti ebmvcoeff for the best template
                                   // fitted for continuum
  TInt32List FittedTplMeiksinIdx;  // Meiksin igm index for the best template
                                   // fitted for continuum
  TFloat64List FittedTplDtm; // DTM for the best template fitted for continuum
  TFloat64List FittedTplMtm; // MTM for the best template fitted for continuum
  TFloat64List
      FittedTplLogPrior; // log prior for the best template fitted for continuum
  TFloat64List FittedTplSNR;

  // Extrema results
  TFloat64List MeritContinuum; // extrema merit for continuum

  TFloat64List mTransposeM;    // extrema model norm
  TFloat64List CorrScaleMarg;  // extrema scale marg. correction
  TInt32List NDof;             // non zero elements in the lambdarange
  TFloat64List Redshift_lmfit; // z found with lmfit
  TFloat64List snrHa;
  TFloat64List lfHa;
  TFloat64List snrOII;
  TFloat64List lfOII;

  std::vector<TFloat64List> ExtendedRedshifts; // z range around extrema
  TFloat64List NLinesOverThreshold;
  TFloat64List LogArea;                 // log area for each extrema
  TFloat64List LogAreaCorrectedExtrema; // corrected z for each extrema
  TFloat64List SigmaZ;                  // sigmaz for each extrema

  TFloat64List StrongELSNR;
  std::vector<std::unordered_set<std::string>> StrongELSNRAboveCut;
  TFloat64List bic; // bayesian information criterion for each extrema
  std::vector<CContinuumIndexes::TContinuumIndexList>
      ContinuumIndexes; // continuum indexes for each extrema
  std::vector<CMask>
      OutsideLinesMask; // Mask with 0 under the lines and 1 anywhere else
  TFloat64List OutsideLinesSTDFlux; // STD measured on the spectrum continuum
                                    // substracted outside lines
  TFloat64List
      OutsideLinesSTDError; // STD measured on the error spectrum outside lines

  // line width
  TFloat64List Elv;                    // emission line width
  TFloat64List Alv;                    // absorption line width
  std::vector<TFloat64List> GroupsELv; // per fitting group line width , EL
  std::vector<TFloat64List> GroupsALv; // per fitting group line width , AL

  // template continuum (+ base class)
  TFloat64List
      FittedTplRedshift; // Redshift for the best template fitted for continuum
  std::vector<TFloat64List> FittedTplpCoeffs; // poly coeffs for the best
                                              // template fitted for continuum

  // template ratio
  TStringList FittedTplratioName;       // Name of the best template fitted for
                                        // tplcorr/tplratio
  TFloat64List FittedTplratioAmplitude; // amp of the best template fitted for
                                        // tplcorr/tplratio
  TFloat64List
      FittedTplratioDtm; // dtm of the best template fitted for tplcorr/tplratio
  TFloat64List
      FittedTplratioMtm; // mtm of the best template fitted for tplcorr/tplratio
  TFloat64List FittedTplratioIsmCoeff; // IsmCoeff/EBMV of the best template
                                       // fitted for tplcorr/tplratio

  mutable std::map<int, TFloat64List> continuumIndexesColorCopy;
  mutable std::map<int, TFloat64List> continuumIndexesBreakCopy;

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