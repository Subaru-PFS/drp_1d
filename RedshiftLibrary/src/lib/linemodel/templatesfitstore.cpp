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
#include "RedshiftLibrary/linemodel/templatesfitstore.h"
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/linemodel/linemodelfitting.h"

#include <float.h>

using namespace NSEpic;

CTemplatesFitStore::CTemplatesFitStore(const TFloat64List &redshifts)
    : redshiftgrid(redshifts) {
  initFitValues();
}

/**
 * @brief CTemplatesFitStore::initFitValues
 * allocate the [nz][n continuum candidates values] structure
 */
void CTemplatesFitStore::initFitValues() {
  std::vector<SValues> zfitvals;
  m_fitValues =
      std::vector<std::vector<SValues>>(redshiftgrid.size(), zfitvals);
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
                             Float64 fitMtM, Float64 logprior) {
  SValues tmpSValues;
  tmpSValues.merit = merit;
  tmpSValues.chiSquare_phot = chiSquare_phot;
  tmpSValues.fitAmplitude = fitAmplitude;
  tmpSValues.fitAmplitudeError = fitAmplitudeError;
  tmpSValues.fitAmplitudeSigma = fitAmplitudeSigma;
  tmpSValues.fitDtM = fitDtM;
  tmpSValues.fitMtM = fitMtM;
  tmpSValues.logprior = logprior;
  tmpSValues.ismEbmvCoeff = ismEbmvCoeff;
  tmpSValues.igmMeiksinIdx = igmMeiksinIdx;
  tmpSValues.tplName = tplName;

  //
  Int32 idxz = GetRedshiftIndex(redshift);
  if (idxz < 0) {
    Log.LogDebug("CTemplatesFitStore::Unable to find z index for redshift=%f",
                 redshift);
    return false;
  }

  // if chi2 val is the lowest, and condition on tplName, insert at position
  // ipos
  auto ipos = std::upper_bound(
      m_fitValues[idxz].begin(), m_fitValues[idxz].end(), tmpSValues.merit,
      [](Float64 val, const SValues &lhs) { return (val < lhs.merit); });

  Log.LogDebug("CTemplatesFitStore::Add iz=%d (z=%f) - adding at pos=%d "
               "(merit=%e, ebmv=%e, imeiksin=%d)",
               idxz, redshift, std::distance(m_fitValues[idxz].begin(), ipos),
               tmpSValues.merit, tmpSValues.ismEbmvCoeff,
               tmpSValues.igmMeiksinIdx);

  // insert the new SValue and move all the older candidates position according
  // to ipos found
  m_fitValues[idxz].insert(ipos, tmpSValues);

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

CTemplatesFitStore::TemplateFitValues
CTemplatesFitStore::GetFitValues(Int32 idxz,
                                 Int32 continuumCandidateRank) const {
  if (continuumCandidateRank > n_continuum_candidates - 1) {
    Log.LogError("CTemplatesFitStore::GetFitValues - cannot find the correct "
                 "pre-computed continuum: candidateRank (%d) >= "
                 "n_continuum_candidates (%d)",
                 continuumCandidateRank, n_continuum_candidates);
    THROWG(INTERNAL_ERROR, "Cannot find the "
                           "correct pre-computed continuum");
  } else if (continuumCandidateRank < 0) {
    Log.LogError("CTemplatesFitStore::GetFitValues - cannot find the correct "
                 "pre-computed continuum: candidateRank (%d) <0",
                 continuumCandidateRank);
    THROWG(INTERNAL_ERROR, "Cannot find the "
                           "correct pre-computed continuum");
  }

  if ((idxz < 0) || (idxz > redshiftgrid.size() - 1)) {
    THROWG(INTERNAL_ERROR,
           Formatter() << "redshift idx " << idxz << " is outside range");
  }

  return m_fitValues[idxz][continuumCandidateRank];
}

CTemplatesFitStore::TemplateFitValues
CTemplatesFitStore::GetFitValues(Float64 redshiftVal,
                                 Int32 continuumCandidateRank) const {
  if (continuumCandidateRank > n_continuum_candidates - 1) {
    Log.LogError("CTemplatesFitStore::GetFitValues - cannot find the correct "
                 "pre-computed continuum: candidateRank (%d) >= "
                 "n_continuum_candidates (%d)",
                 continuumCandidateRank, n_continuum_candidates);
    THROWG(INTERNAL_ERROR, "Cannot find the "
                           "correct pre-computed continuum");

  } else if (continuumCandidateRank < 0) {
    Log.LogError("CTemplatesFitStore::GetFitValues - cannot find the correct "
                 "pre-computed continuum: candidateRank (%d) <0",
                 continuumCandidateRank);
    THROWG(INTERNAL_ERROR, "Cannot find the "
                           "correct pre-computed continuum");
  }

  if (redshiftVal < redshiftgrid[0]) {
    Log.LogError("CTemplatesFitStore - GetFitValues, looking for "
                 "redshiftVal=%f, but lt redshiftgrid[0]=%f",
                 redshiftVal, redshiftgrid[0]);
    THROWG(INTERNAL_ERROR, "Looking for "
                           "outside range redshiftVal");
  }
  if (redshiftVal > redshiftgrid[redshiftgrid.size() - 1]) {
    Log.LogError(
        "CTemplatesFitStore - GetFitValues, looking for redshiftVal=%f, but ht "
        "redshiftgrid[redshiftgrid.size()-1]=%f",
        redshiftVal, redshiftgrid[redshiftgrid.size() - 1]);
    THROWG(INTERNAL_ERROR, "Looking for "
                           "outside range redshiftVal");
  }

  // find the idxz
  Int32 idxz = GetRedshiftIndex(redshiftVal);

  if (idxz < 0) {
    Log.LogError("CTemplatesFitStore::GetFitValues - cannot find the correct "
                 "pre-computed continuum.");
    Log.LogError("CTemplatesFitStore::GetFitValues - redshiftVal=%e, but lt "
                 "m_fitValues[0].redshift=%e",
                 redshiftVal, redshiftgrid[0]);
    Log.LogError("CTemplatesFitStore::GetFitValues - redshiftVal=%e, but ht "
                 "m_fitValues[m_fitValues.size()-1].redshift=%e",
                 redshiftVal, redshiftgrid[redshiftgrid.size() - 1]);

    for (Int32 k = 0; k < 10; k++) {
      Log.LogDebug("CTemplatesFitStore::GetFitValues - redshiftVal=%e, lt "
                   "m_fitValues[%d].redshift=%e",
                   redshiftVal, k, redshiftgrid[k]);
    }

    THROWG(INTERNAL_ERROR, "Cannot find redshiftVal");
  }

  return m_fitValues[idxz][continuumCandidateRank];
}

Float64
CTemplatesFitStore::FindMaxAmplitudeSigma(Float64 &z,
                                          TemplateFitValues &fitValues) {
  Int32 icontinuum = 0;
  m_fitContinuum_fitAmplitudeSigmaMAX = -INFINITY;
  // TemplateFitValues fitValues;
  for (Int32 i = 0; i < redshiftgrid.size(); i++) {
    const TemplateFitValues &thisfitValues = m_fitValues[i][icontinuum];
    if (thisfitValues.fitAmplitudeSigma > m_fitContinuum_fitAmplitudeSigmaMAX) {
      m_fitContinuum_fitAmplitudeSigmaMAX = thisfitValues.fitAmplitudeSigma;
      z = redshiftgrid[i];
      fitValues = thisfitValues;
    }
  }

  return m_fitContinuum_fitAmplitudeSigmaMAX;
}
