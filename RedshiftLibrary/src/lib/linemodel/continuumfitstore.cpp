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

CContinuumFitStore::CContinuumFitStore(const TFloat64List &redshifts)
    : m_redshiftgrid(redshifts) {}

Int32 CContinuumFitStore::GetRedshiftIndex(Float64 z) const {
  auto it = std::find(m_redshiftgrid.begin(), m_redshiftgrid.end(), z);
  if (it != m_redshiftgrid.end())
    return std::distance(m_redshiftgrid.begin(), it);
  else
    return -1;
}

Int32 CContinuumFitStore::getClosestLowerRedshiftIndex(Float64 z) const {
  Int32 idx = -1;
  CIndexing<Float64>::getClosestLowerIndex(m_redshiftgrid, z, idx);
  return idx;
}

const TFloat64List &CContinuumFitStore::GetRedshiftList() const {
  return m_redshiftgrid;
}

void CContinuumFitStore::Add(std::string tplName, Float64 ismEbmvCoeff,
                             Int32 igmMeiksinIdx, Float64 redshift,
                             Float64 merit, Float64 chiSquare_phot,
                             Float64 fitAmplitude, Float64 fitAmplitudeError,
                             Float64 fitAmplitudeSigma, Float64 fitDtM,
                             Float64 fitMtM, Float64 logprior, Float64 snr) {
  THROWG(ErrorCode::INTERNAL_ERROR, Formatter() << "Wrong method");
}

void CContinuumFitStore::Add(
    Float64 ismEbmvCoeff, Int32 igmMeiksinIdx, Float64 redshift,
    Float64 chi2, // TODO see if chi2 and merit is the same
    Float64 a1, Float64 a2, Float64 b1, Float64 b2, Float64 snr) {
  THROWG(ErrorCode::INTERNAL_ERROR, Formatter() << "Wrong method");
}

const CContinuumModelSolution &
CContinuumFitStore::GetFitValues(Int32 idxz,
                                 Int32 continuumCandidateRank) const {
  if (continuumCandidateRank > n_continuum_candidates - 1)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Cannot find the "
                          "correct pre-computed continuum: candidateRank ("
                       << continuumCandidateRank
                       << ") >= n_continuum_candidates ("
                       << n_continuum_candidates << ")");
  if (continuumCandidateRank < 0)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Cannot find the "
                          "correct pre-computed continuum: candidateRank ("
                       << continuumCandidateRank << ") <0");

  if ((idxz < 0) || (idxz > m_redshiftgrid.size() - 1)) {
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "redshift idx " << idxz << " is outside range");
  }

  return m_fitValues[idxz][continuumCandidateRank];
}

const CContinuumModelSolution &
CContinuumFitStore::GetFitValues(Float64 redshiftVal,
                                 Int32 continuumCandidateRank) const {
  if (continuumCandidateRank > n_continuum_candidates - 1)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Cannot find the "
                          "correct pre-computed continuum: candidateRank ("
                       << continuumCandidateRank
                       << ") >= n_continuum_candidates ("
                       << n_continuum_candidates << ")");

  if (continuumCandidateRank < 0)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Cannot find the "
                          "correct pre-computed continuum: candidateRank("
                       << continuumCandidateRank << ") <0");

  if (redshiftVal < m_redshiftgrid.front() ||
      redshiftVal > m_redshiftgrid.back())
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Looking for redshiftVal=" << redshiftVal
                       << " outside range [" << m_redshiftgrid.front() << ", "
                       << m_redshiftgrid.back() << "]");

  // find the idxz
  Int32 idxz = GetRedshiftIndex(redshiftVal);

  if (idxz < 0) {
    Log.LogDetail(
        Formatter()
        << "CTemplatesFitStore::GetFitValues - cannot find the correct "
           "pre-computed continuum.");
    Log.LogDetail(Formatter()
                  << "CTemplatesFitStore::GetFitValues - redshiftVal="
                  << redshiftVal
                  << ", but lt "
                     "m_fitValues[0].redshift="
                  << m_redshiftgrid[0]);
    Log.LogDetail(Formatter()
                  << "CTemplatesFitStore::GetFitValues - redshiftVal="
                  << redshiftVal
                  << ", but ht "
                     "m_fitValues[m_fitValues.size()-1].redshift="
                  << m_redshiftgrid[m_redshiftgrid.size() - 1]);

    for (Int32 k = 0; k < 10; k++)
      Log.LogDebug(Formatter()
                   << "CTemplatesFitStore::GetFitValues - redshiftVal="
                   << redshiftVal
                   << ", lt "
                      "m_fitValues["
                   << k << "].redshift=" << m_redshiftgrid[k]);

    THROWG(ErrorCode::INTERNAL_ERROR, "Cannot find redshiftVal");
  }
  // CContinuumModelSolution cms = m_fitValues[idxz][continuumCandidateRank];
  //   cms.tplRedshift = redshiftVal;
  return m_fitValues[idxz][continuumCandidateRank];
}

/**
 * @brief CTemplatesFitStore::initFitValues
 * allocate the [nz][n continuum candidates values] structure
 */
void CContinuumFitStore::initFitValues() {
  std::vector<CContinuumModelSolution> zfitvals;
  m_fitValues = std::vector<std::vector<CContinuumModelSolution>>(
      m_redshiftgrid.size(), zfitvals);
}