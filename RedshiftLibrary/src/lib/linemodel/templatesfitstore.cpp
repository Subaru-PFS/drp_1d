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
#include <climits>

#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/linemodel/linemodelfitting.h"
#include "RedshiftLibrary/linemodel/templatesfitstore.h"

using namespace NSEpic;

CTemplatesFitStore::CTemplatesFitStore(const TFloat64List &redshifts)
    : redshiftgrid(redshifts) {
  initFitValues();
  m_fitMaxValues = std::make_shared<fitMaxValues>();
}

/**
 * @brief CTemplatesFitStore::initFitValues
 * allocate the [nz][n continuum candidates values] structure
 */
void CTemplatesFitStore::initFitValues() {
  std::vector<CTplModelSolution> zfitvals;
  m_fitValues = std::vector<std::vector<CTplModelSolution>>(redshiftgrid.size(),
                                                            zfitvals);
}

const TFloat64List &CTemplatesFitStore::GetRedshiftList() const {
  return redshiftgrid;
}

Int32 CTemplatesFitStore::GetRedshiftIndex(Float64 z) const {
  auto it = std::find(redshiftgrid.begin(), redshiftgrid.end(), z);
  if (it != redshiftgrid.end())
    return std::distance(redshiftgrid.begin(), it);
  else
    return -1;
}

Int32 CTemplatesFitStore::getClosestLowerRedshiftIndex(Float64 z) const {
  Int32 idx = -1;
  CIndexing<Float64>::getClosestLowerIndex(redshiftgrid, z, idx);
  return idx;
}

/**
 * @brief CTemplatesFitStore::Add
 * @param tplName
 * @param ismEbmvCoeff
 * @param igmMeiksinIdx
 * @param redshift
 * @param merit
 * @param fitAmplitude
 * @param fitAmplitudeError
 * @param fitDtM
 * @param fitMtM
 *
 * brief: try to insert the fit values into the mFitValues table at the correct
 * idxz:
 *   - no insertion if redshift can't be found in the redshiftgrid
 *   - no insertion if the merit is higher than the highest rank continuum
 * candidate
 *   - insertion is done at a given continuum_candidate_rank position wrt merit
 * value
 * @return False if there was a problem.
 */
bool CTemplatesFitStore::Add(std::string tplName, Float64 ismEbmvCoeff,
                             Int32 igmMeiksinIdx, Float64 redshift,
                             Float64 merit, Float64 chiSquare_phot,
                             Float64 fitAmplitude, Float64 fitAmplitudeError,
                             Float64 fitAmplitudeSigma, Float64 fitDtM,
                             Float64 fitMtM, Float64 logprior, Float64 snr) {
  CTplModelSolution tmpCContinuumModelSolution;
  tmpCContinuumModelSolution.tplMerit = merit;
  tmpCContinuumModelSolution.tplMeritPhot = chiSquare_phot;
  tmpCContinuumModelSolution.tplAmplitude = fitAmplitude;
  tmpCContinuumModelSolution.tplAmplitudeError = fitAmplitudeError;
  tmpCContinuumModelSolution.tplAmplitudeSigma = fitAmplitudeSigma;
  tmpCContinuumModelSolution.tplDtM = fitDtM;
  tmpCContinuumModelSolution.tplMtM = fitMtM;
  tmpCContinuumModelSolution.tplLogPrior = logprior;
  tmpCContinuumModelSolution.tplEbmvCoeff = ismEbmvCoeff;
  tmpCContinuumModelSolution.tplMeiksinIdx = igmMeiksinIdx;
  tmpCContinuumModelSolution.tplName = tplName;
  tmpCContinuumModelSolution.tplRedshift = redshift;
  tmpCContinuumModelSolution.tplSNR = snr;

  //
  Int32 idxz = GetRedshiftIndex(redshift);
  if (idxz < 0) {
    Log.LogDebug("CTemplatesFitStore::Unable to find z index for redshift=%f",
                 redshift);
    return false;
  }

  // if chi2 val is the lowest, and condition on tplName, insert at position
  // ipos
  auto ipos =
      std::upper_bound(m_fitValues[idxz].begin(), m_fitValues[idxz].end(),
                       tmpCContinuumModelSolution.tplMerit,
                       [](Float64 val, const CTplModelSolution &lhs) {
                         return (val < lhs.tplMerit);
                       });

  Log.LogDebug("CTemplatesFitStore::Add iz=%d (z=%f) - adding at pos=%d "
               "(merit=%e, ebmv=%e, imeiksin=%d)",
               idxz, redshift, std::distance(m_fitValues[idxz].begin(), ipos),
               tmpCContinuumModelSolution.tplMerit,
               tmpCContinuumModelSolution.tplEbmvCoeff,
               tmpCContinuumModelSolution.tplMeiksinIdx);

  // insert the new SValue and move all the older candidates position according
  // to ipos found
  m_fitValues[idxz].insert(ipos, tmpCContinuumModelSolution);

  // this is not very secure. it should be checked that all redshifts have the
  // same fitValues count
  if (n_continuum_candidates < m_fitValues[idxz].size()) {
    n_continuum_candidates = m_fitValues[idxz].size();
    Log.LogDebug("CTemplatesFitStore::n_continuum_candidates set to %d)",
                 n_continuum_candidates);
  }

  return true;
}

Int32 CTemplatesFitStore::GetContinuumCount() const {
  return n_continuum_candidates;
}

const CTplModelSolution &
CTemplatesFitStore::GetFitValues(Int32 idxz,
                                 Int32 continuumCandidateRank) const {
  if (continuumCandidateRank > n_continuum_candidates - 1)
    THROWG(INTERNAL_ERROR,
           Formatter() << "Cannot find the "
                          "correct pre-computed continuum: candidateRank ("
                       << continuumCandidateRank
                       << ") >= n_continuum_candidates ("
                       << n_continuum_candidates << ")");
  if (continuumCandidateRank < 0)
    THROWG(INTERNAL_ERROR,
           Formatter() << "Cannot find the "
                          "correct pre-computed continuum: candidateRank ("
                       << continuumCandidateRank << ") <0");

  if ((idxz < 0) || (idxz > redshiftgrid.size() - 1)) {
    THROWG(INTERNAL_ERROR,
           Formatter() << "redshift idx " << idxz << " is outside range");
  }

  return m_fitValues[idxz][continuumCandidateRank];
}

const CTplModelSolution &
CTemplatesFitStore::GetFitValues(Float64 redshiftVal,
                                 Int32 continuumCandidateRank) const {
  if (continuumCandidateRank > n_continuum_candidates - 1)
    THROWG(INTERNAL_ERROR,
           Formatter() << "Cannot find the "
                          "correct pre-computed continuum: candidateRank ("
                       << continuumCandidateRank
                       << ") >= n_continuum_candidates ("
                       << n_continuum_candidates << ")");

  if (continuumCandidateRank < 0)
    THROWG(INTERNAL_ERROR,
           Formatter() << "Cannot find the "
                          "correct pre-computed continuum: candidateRank("
                       << continuumCandidateRank << ") <0");

  if (redshiftVal < redshiftgrid.front() || redshiftVal > redshiftgrid.back())
    THROWG(INTERNAL_ERROR, Formatter()
                               << "Looking for redshiftVal=" << redshiftVal
                               << " outside range [" << redshiftgrid.front()
                               << ", " << redshiftgrid.back() << "]");

  // find the idxz
  Int32 idxz = GetRedshiftIndex(redshiftVal);

  if (idxz < 0) {
    Log.LogDetail("CTemplatesFitStore::GetFitValues - cannot find the correct "
                  "pre-computed continuum.");
    Log.LogDetail("CTemplatesFitStore::GetFitValues - redshiftVal=%e, but lt "
                  "m_fitValues[0].redshift=%e",
                  redshiftVal, redshiftgrid[0]);
    Log.LogDetail("CTemplatesFitStore::GetFitValues - redshiftVal=%e, but ht "
                  "m_fitValues[m_fitValues.size()-1].redshift=%e",
                  redshiftVal, redshiftgrid[redshiftgrid.size() - 1]);

    for (Int32 k = 0; k < 10; k++)
      Log.LogDebug("CTemplatesFitStore::GetFitValues - redshiftVal=%e, lt "
                   "m_fitValues[%d].redshift=%e",
                   redshiftVal, k, redshiftgrid[k]);

    THROWG(INTERNAL_ERROR, "Cannot find redshiftVal");
  }
  // CTplModelSolution cms = m_fitValues[idxz][continuumCandidateRank];
  //   cms.tplRedshift = redshiftVal;
  return m_fitValues[idxz][continuumCandidateRank];
}

Float64
CTemplatesFitStore::FindMaxAmplitudeSigma(Float64 &z,
                                          CTplModelSolution &fitValues) {
  Int32 icontinuum = 0;
  m_fitMaxValues->fitAmplitudeSigmaMAX = -INFINITY;
  // TemplateFitValues fitValues;
  for (Int32 i = 0; i < redshiftgrid.size(); i++) {
    const CTplModelSolution &thisfitValues = m_fitValues[i][icontinuum];
    if (thisfitValues.tplAmplitudeSigma >
        m_fitMaxValues->fitAmplitudeSigmaMAX) {
      m_fitMaxValues->fitAmplitudeSigmaMAX = thisfitValues.tplAmplitudeSigma;
      z = redshiftgrid[i];
      fitValues = thisfitValues;
    }
  }

  return m_fitMaxValues->fitAmplitudeSigmaMAX;
}
