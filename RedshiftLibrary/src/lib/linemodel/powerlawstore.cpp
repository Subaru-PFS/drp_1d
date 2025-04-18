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

void CPowerLawStore::Add(Float64 ismEbmvCoeff, Int32 igmMeiksinIdx,
                         Float64 redshift, Float64 chi2, Float64 reducedChi2,
                         Float64 pValue, TPowerLawCoefsPair powerLawCoefs,
                         Float64 snr) {
  CContinuumModelSolution tmpCContinuumModelSolution;
  tmpCContinuumModelSolution.name = "powerLaw";
  tmpCContinuumModelSolution.a1 = powerLawCoefs.first.a;
  tmpCContinuumModelSolution.a2 = powerLawCoefs.second.a;
  tmpCContinuumModelSolution.b1 = powerLawCoefs.first.b;
  tmpCContinuumModelSolution.b2 = powerLawCoefs.second.b;
  tmpCContinuumModelSolution.a1std = powerLawCoefs.first.stda;
  tmpCContinuumModelSolution.a2std = powerLawCoefs.second.stda;
  tmpCContinuumModelSolution.b1std = powerLawCoefs.first.stdb;
  tmpCContinuumModelSolution.b2std = powerLawCoefs.second.stdb;
  tmpCContinuumModelSolution.ebmvCoef = ismEbmvCoeff;
  tmpCContinuumModelSolution.meiksinIdx = igmMeiksinIdx;
  tmpCContinuumModelSolution.redshift = redshift;
  tmpCContinuumModelSolution.merit = chi2;
  tmpCContinuumModelSolution.reducedChi2 = reducedChi2;
  tmpCContinuumModelSolution.pValue = pValue;
  tmpCContinuumModelSolution.SNR = snr;

  Int32 idxz = GetRedshiftIndex(redshift);
  if (idxz < 0)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "Unable to find z index for redshift=" << redshift);

  m_fitValues[idxz].push_back(std::move(tmpCContinuumModelSolution));
};

Float64 CPowerLawStore::getFracAmplitudeSigma(
    CContinuumModelSolution const &continuum) const {
  return std::max(continuum.a1 / continuum.a1std,
                  continuum.a2 / continuum.a2std);
};
