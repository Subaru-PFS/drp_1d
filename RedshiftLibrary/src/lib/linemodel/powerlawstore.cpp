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
#include "RedshiftLibrary/linemodel/powerlawstore.h"

using namespace NSEpic;

CPowerLawStore::CPowerLawStore(const TFloat64List &redshifts)
    : CContinuumFitStore(redshifts) {
  initFitValues();
  m_fitMaxValues = std::make_shared<fitMaxValues>();
}

void CPowerLawStore::Add(Float64 ismEbmvCoeff, Int32 igmMeiksinIdx,
                         Float64 redshift,
                         Float64 chi2, // TODO see if chi2 and merit is the same
                         Float64 a1, Float64 a2, Float64 b1, Float64 b2,
                         Float64 snr) {
  CContinuumModelSolution tmpCContinuumModelSolution;
  tmpCContinuumModelSolution.a1 = a1;
  tmpCContinuumModelSolution.a2 = a2;
  tmpCContinuumModelSolution.b1 = b1;
  tmpCContinuumModelSolution.b2 = b2;
  tmpCContinuumModelSolution.EbmvCoeff = ismEbmvCoeff;
  tmpCContinuumModelSolution.MeiksinIdx = igmMeiksinIdx;
  tmpCContinuumModelSolution.redshift = redshift;
  tmpCContinuumModelSolution.tplMerit = chi2;
  tmpCContinuumModelSolution.SNR = snr;

  Int32 idxz = GetRedshiftIndex(redshift);
  if (idxz < 0)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Unable to find z index for redshift=" << redshift);

  // if chi2 val is the lowest, and condition on tplName, insert at position
  // ipos
  // TODO : snr is incorrect should find an equivalent of amplitude sigma
  auto ipos =
      std::upper_bound(m_fitValues[idxz].begin(), m_fitValues[idxz].end(),
                       std::make_pair(tmpCContinuumModelSolution.tplMerit,
                                      std::abs(tmpCContinuumModelSolution.SNR)),
                       [](std::pair<Float64, Float64> const &val,
                          CContinuumModelSolution const &lhs) {
                         if (val.first != lhs.tplMerit)
                           return (val.first < lhs.tplMerit);
                         return val.second > std::abs(lhs.SNR);
                       });

  Log.LogDebug(
      Formatter() << "CTemplatesFitStore::Add iz=" << idxz << " (z=" << redshift
                  << ") - adding at pos="
                  << std::distance(m_fitValues[idxz].begin(), ipos)
                  << "(merit=" << tmpCContinuumModelSolution.tplMerit
                  << ", ebmv=" << tmpCContinuumModelSolution.tplEbmvCoeff
                  << ", imeiksin=" << tmpCContinuumModelSolution.tplMeiksinIdx
                  << ")");

  // insert the new SValue and move all the older candidates position according
  // to ipos found
  m_fitValues[idxz].insert(ipos, tmpCContinuumModelSolution);

  // this is not very secure. it should be checked that all redshifts have the
  // same fitValues count
  if (n_continuum_candidates < m_fitValues[idxz].size()) {
    n_continuum_candidates = m_fitValues[idxz].size();
    Log.LogDebug(Formatter()
                 << "CTemplatesFitStore::n_continuum_candidates set to "
                 << n_continuum_candidates << ")");
  }
};
