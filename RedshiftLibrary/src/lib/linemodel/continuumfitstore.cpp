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
#include "RedshiftLibrary/linemodel/continuumfitstore.h"
#include "RedshiftLibrary/linemodel/linemodelfitting.h"

using namespace NSEpic;

CContinuumFitStore::CContinuumFitStore(const TFloat64List &redshifts)
    : m_redshiftgrid(redshifts) {
  initFitValues();
}

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

const CContinuumModelSolution &
CContinuumFitStore::GetFitValues(Int32 idxz,
                                 Int32 continuumCandidateRank) const {
  Int32 continuumCount = getContinuumCount();
  if (continuumCandidateRank > continuumCount - 1)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Cannot find the "
                          "correct pre-computed continuum: candidateRank ("
                       << continuumCandidateRank << ") >= continuumCount ("
                       << continuumCount << ")");
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

std::vector<std::vector<CContinuumModelSolution>>
CContinuumFitStore::GetFitValuesVector() {
  return m_fitValues;
}

const CContinuumModelSolution &
CContinuumFitStore::GetFitValues(Float64 redshiftVal,
                                 Int32 continuumCandidateRank) const {
  Int32 continuumCount = getContinuumCount();
  if (continuumCandidateRank > continuumCount - 1)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Cannot find the "
                          "correct pre-computed continuum: candidateRank ("
                       << continuumCandidateRank << ") >= continuumCount ("
                       << continuumCount << ")");

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

std::pair<Float64, CContinuumModelSolution const>
CContinuumFitStore::FindMaxAmplitudeSigma() const {
  Int32 const icontinuum = 0;
  Float64 fitAmplitudeSigmaMAX = -INFINITY;
  Int32 i_max = undefIdx;
  for (Int32 i = 0; i < m_redshiftgrid.size(); i++) {
    const CContinuumModelSolution &continuum = m_fitValues[i][icontinuum];
    Float64 fracAmplitudeSigma = getFracAmplitudeSigma(continuum);
    if (fracAmplitudeSigma > fitAmplitudeSigmaMAX) {
      fitAmplitudeSigmaMAX = fracAmplitudeSigma;
      i_max = i;
    }
  }

  return std::make_pair(fitAmplitudeSigmaMAX,
                        std::cref(m_fitValues[i_max][icontinuum]));
}

CContinuumModelSolution const &CContinuumFitStore::FindMinReducedChi2() const {
  Int32 const icontinum = 0;
  Float64 max_reduced_chi2 = -INFINITY;
  auto const max_elt = std::min_element(
      m_fitValues.cbegin(), m_fitValues.cend(),
      [this](std::vector<CContinuumModelSolution> const &left_fitValues,
             std::vector<CContinuumModelSolution> const &right_fitValues) {
        return left_fitValues[icontinum].reducedChi2 <
               right_fitValues[icontinum].reducedChi2;
      });
  return (*max_elt)[icontinum];
}
