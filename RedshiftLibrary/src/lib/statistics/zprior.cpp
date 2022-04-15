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
#include "RedshiftLibrary/statistics/zprior.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/operator/pdfz.h"

#include "RedshiftLibrary/log/log.h"
#include <float.h>

using namespace NSEpic;
using namespace std;

void CZPrior::NormalizePrior(TFloat64List &logzPrior) const {
  Float64 logSum = COperatorPdfz::logSumExpTrick(logzPrior, m_redshifts);
  for (Float64 &zp : logzPrior)
    zp -= logSum;
}

// setting cte priors for all redshift values
TFloat64List CZPrior::GetConstantLogZPrior(Int32 nredshifts) const {
  TFloat64List logzPrior(nredshifts, 0.0);

  if (m_normalizePrior)
    NormalizePrior(logzPrior);

  return logzPrior;
}

TFloat64List CZPrior::GetStrongLinePresenceLogZPrior(
    const TBoolList &linePresence, const Float64 penalization_factor) const {
  Float64 logprobaPresent = 0.0;
  Float64 logprobaAbsent = log(penalization_factor);
  TFloat64List logzPrior(linePresence.size(), logprobaAbsent);
  for (Int32 kz = 0; kz < linePresence.size(); kz++) {
    if (linePresence[kz]) {
      logzPrior[kz] = logprobaPresent;
    }
  }

  if (m_normalizePrior)
    NormalizePrior(logzPrior);

  return logzPrior;
}

TFloat64List CZPrior::GetNLinesSNRAboveCutLogZPrior(
    const TInt32List &nlinesAboveSNR, const Float64 penalization_factor) const {
  Int32 nz = nlinesAboveSNR.size();
  Int32 nlinesThres = 2;
  Float64 logprobaPresent = 0.0;
  Float64 logprobaAbsent = log(penalization_factor);
  TFloat64List logzPrior(nz, logprobaAbsent);
  for (Int32 kz = 0; kz < nz; kz++) {
    if (nlinesAboveSNR[kz] >= nlinesThres) {
      logzPrior[kz] = logprobaPresent;
      // Log.LogDetail("ZPrior: Prior: nlinesAboveSNR[kz] >= nlinesThres for
      // kz=%d", kz);
    }
  }

  if (m_normalizePrior)
    NormalizePrior(logzPrior);

  return logzPrior;
}

TFloat64List CZPrior::GetEuclidNhaLogZPrior(const TRedshiftList &redshifts,
                                            const Float64 aCoeff) const {
  if (aCoeff <= 0.0) {
    throw GlobalException(ErrorCode::INTERNAL_ERROR,
                          Formatter() << "CZPrior::GetEuclidNhaLogZPrior: "
                                         "problem found aCoeff<=0: aCoeff="
                                      << aCoeff);
  }

  TFloat64List zPrior(redshifts.size(), 0.0);
  TFloat64List logzPrior(redshifts.size(), -INFINITY);

  Float64 maxP = -DBL_MAX;
  Float64 minP = DBL_MAX;
  for (Int32 kz = 0; kz < redshifts.size(); kz++) {
    Float64 z = redshifts[kz];
    Float64 z2 = z * z;
    Float64 z3 = z2 * z;
    Float64 z4 = z3 * z;
    Float64 z5 = z4 * z;
    Float64 z6 = z5 * z;

    // poly reg pozzetti model 1 at FHa=1e-16
    zPrior[kz] = (-54.7422088727874 * z6 + 1203.94994364807 * z5 -
                  10409.6716744981 * z4 + 44240.3837462642 * z3 -
                  92914.84430357 * z2 + 79004.76406 * z - 2288.98457865);

    // shape prior at low z, left of the bell
    bool enable_low_z_flat = false;
    if (enable_low_z_flat && z < 0.7204452872044528) {
      zPrior[kz] = 20367.877916402278;
    } else if (zPrior[kz] < 0) {
      zPrior[kz] = DBL_MIN;
    }

    // apply strength & switch to log
    logzPrior[kz] = log(zPrior[kz]) * aCoeff;

    if (logzPrior[kz] > maxP) {
      maxP = logzPrior[kz];
    }
    if (logzPrior[kz] < minP) {
      minP = logzPrior[kz];
    }
  }

  Log.LogDebug("Pdfz: log zPrior: using HalphaZPrior min=%e", minP);
  Log.LogDebug("Pdfz: log zPrior: using HalphaZPrior max=%e", maxP);

  Float64 log_dynamicCut = -log(1e12);
  if (maxP == maxP && maxP != -INFINITY) // test NAN & -inf
  {
    for (Float64 &zP : logzPrior) {
      zP -= maxP;
      if (zP < log_dynamicCut)
        zP = log_dynamicCut;
    }
  }

  if (m_normalizePrior)
    NormalizePrior(logzPrior);

  return logzPrior;
}

/**
 * @brief CombineLogZPrior
 * returns a vector with log(prior1*prior2) for each z, with normalization so
 * that sumPrior=1
 * @param logprior1
 * @param logprior2
 * @return
 */
TFloat64List CZPrior::CombineLogZPrior(const TFloat64List &logprior1,
                                       const TFloat64List &logprior2) const {
  if (logprior1.size() != logprior2.size()) {
    throw GlobalException(
        ErrorCode::INTERNAL_ERROR,
        "CZPrior::CombineLogZPrior, the 2 zpriors have not the same size");
  }
  Int32 n = logprior1.size();

  TFloat64List logzPriorCombined(n);

  for (Int32 k = 0; k < n; k++) {
    Float64 val = logprior1[k] + logprior2[k];
    logzPriorCombined[k] = val;
  }

  if (m_normalizePrior)
    NormalizePrior(logzPriorCombined);

  return logzPriorCombined;
}
