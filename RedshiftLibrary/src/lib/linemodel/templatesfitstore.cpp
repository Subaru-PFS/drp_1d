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
 *   - no insertion if redshift can't be found in the m_redshiftgrid
 *   - no insertion if the merit is higher than the highest rank continuum
 * candidate
 *   - insertion is done at a given continuum_candidate_rank position wrt merit
 * value
 * @return False if there was a problem.
 */
void CTemplatesFitStore::Add(std::string tplName, Float64 ismEbmvCoeff,
                             Int32 igmMeiksinIdx, Float64 redshift,
                             Float64 merit, Float64 chiSquare_phot,
                             Float64 fitAmplitude, Float64 fitAmplitudeError,
                             Float64 fitAmplitudeSigma, Float64 fitDtM,
                             Float64 fitMtM, Float64 logprior, Float64 snr) {
  CContinuumModelSolution tmpCContinuumModelSolution;
  tmpCContinuumModelSolution.merit = merit;
  tmpCContinuumModelSolution.tplMeritPhot = chiSquare_phot;
  tmpCContinuumModelSolution.tplAmplitude = fitAmplitude;
  tmpCContinuumModelSolution.tplAmplitudeError = fitAmplitudeError;
  tmpCContinuumModelSolution.tplAmplitudeSigma = fitAmplitudeSigma;
  tmpCContinuumModelSolution.tplDtM = fitDtM;
  tmpCContinuumModelSolution.tplMtM = fitMtM;
  tmpCContinuumModelSolution.tplLogPrior = logprior;
  tmpCContinuumModelSolution.ebmvCoef = ismEbmvCoeff;
  tmpCContinuumModelSolution.meiksinIdx = igmMeiksinIdx;
  tmpCContinuumModelSolution.tplName = tplName;
  tmpCContinuumModelSolution.redshift = redshift;
  tmpCContinuumModelSolution.SNR = snr;

  Int32 idxz = GetRedshiftIndex(redshift);
  if (idxz < 0)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Unable to find z index for redshift=" << redshift);

  // if chi2 val is the lowest, and condition on tplName, insert at position
  // ipos
  auto ipos = std::upper_bound(
      m_fitValues[idxz].begin(), m_fitValues[idxz].end(),
      std::make_pair(tmpCContinuumModelSolution.merit,
                     std::abs(tmpCContinuumModelSolution.tplAmplitudeSigma)),
      [](std::pair<Float64, Float64> const &val,
         CContinuumModelSolution const &lhs) {
        if (val.first != lhs.merit)
          return (val.first < lhs.merit);
        return val.second > std::abs(lhs.tplAmplitudeSigma);
      });

  Log.LogDebug(Formatter() << "CTemplatesFitStore::Add iz=" << idxz
                           << " (z=" << redshift << ") - adding at pos="
                           << std::distance(m_fitValues[idxz].begin(), ipos)
                           << "(merit=" << tmpCContinuumModelSolution.merit
                           << ", ebmv=" << tmpCContinuumModelSolution.ebmvCoef
                           << ", imeiksin="
                           << tmpCContinuumModelSolution.meiksinIdx << ")");

  // insert the new SValue and move all the older candidates position according
  // to ipos found
  m_fitValues[idxz].insert(ipos, tmpCContinuumModelSolution);

  if (getContinuumCount() < m_fitValues[idxz].size()) {
    m_nContinuumCandidates = m_fitValues[idxz].size();
    Log.LogDebug(Formatter()
                 << "CTemplatesFitStore::m_nContinuumCandidates set to "
                 << getContinuumCount() << ")");
  }
}

Int32 CTemplatesFitStore::getContinuumCount() const {
  return m_nContinuumCandidates;
}

Float64 CTemplatesFitStore::getFracAmplitudeSigma(CContinuumModelSolution const &continuum) const {
  return continuum.tplAmplitudeSigma;
}