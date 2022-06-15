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
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/statistics/zprior.h"
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <math.h>

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(Statistics_zprior)

TFloat64List genRedshifts(Int32 size) {
  TFloat64List redshifts(size);
  for (Int32 i = 0; i < redshifts.size(); i++) {
    redshifts[i] = i * 0.1;
  }
  return redshifts;
}

Float64 precision = 1e-12;

BOOST_AUTO_TEST_CASE(NormalizePrior_test) {
  // Range < 2 -> Global exception in COperatorPdfz::logSumExpTrick
  TFloat64List redshifts = genRedshifts(1);
  CZPrior zprior = CZPrior(true, redshifts);
  TFloat64List logzPrior(redshifts.size(), 0.0);
  BOOST_CHECK_THROW(zprior.NormalizePrior(logzPrior), GlobalException);

  // Range >= 2 -> OK
  redshifts = genRedshifts(5);
  TFloat64List logzPrior_2 = {0., 0., 1., 0., 0.};
  BOOST_CHECK_NO_THROW(zprior.NormalizePrior(logzPrior_2));

  // Check that normalization is OK
  Float64 logSum = COperatorPdfz::logSumExpTrick(logzPrior_2, redshifts);
  BOOST_CHECK_CLOSE(exp(logSum), 1., precision);
}

BOOST_AUTO_TEST_CASE(GetConstantLogZPrior_test) {
  // GetConstantLogZPrior without normalization
  TFloat64List redshifts = genRedshifts(9);
  CZPrior zprior = CZPrior(false, redshifts);
  TFloat64List logzPrior;

  logzPrior = zprior.GetConstantLogZPrior(redshifts.size());
  for (Int32 i = 1; i < logzPrior.size(); i++) {
    BOOST_CHECK_CLOSE((logzPrior[i] - logzPrior[i - 1]), 0.0, precision);
  }

  // GetConstantLogZPrior with normalization
  CZPrior zprior_2 = CZPrior(true, redshifts);
  logzPrior = zprior_2.GetConstantLogZPrior(redshifts.size());
  for (Int32 i = 1; i < logzPrior.size(); i++) {
    BOOST_CHECK_CLOSE((logzPrior[i] - logzPrior[i - 1]), 0.0, precision);
  }

  // Check that normalization is OK
  Float64 logSum = COperatorPdfz::logSumExpTrick(logzPrior, redshifts);
  BOOST_CHECK_CLOSE(exp(logSum), 1., precision);
}

BOOST_AUTO_TEST_CASE(GetStrongLinePresenceLogZPrior_test) {

  TFloat64List redshifts = genRedshifts(9);
  CZPrior zprior = CZPrior(false, redshifts);
  TBoolList linePresence = {0, 0, 0, 1, 1, 1, 0, 0, 0};
  Float64 penalization_factor = 0.5;

  // GetStrongLinePresenceLogZPrior
  // without normalization
  TFloat64List logzPrior =
      zprior.GetStrongLinePresenceLogZPrior(linePresence, penalization_factor);
  for (Int32 i = 0; i < logzPrior.size(); i++) {
    if (i == 3 || i == 4 || i == 5)
      BOOST_CHECK(logzPrior[i] == 0.);
    else
      BOOST_CHECK_CLOSE(logzPrior[i], log(penalization_factor), precision);
  }
  zprior.NormalizePrior(logzPrior);

  // GetStrongLinePresenceLogZPrior
  // with normalization
  CZPrior zprior_2 = CZPrior(true, redshifts);
  TFloat64List logzPrior_2 = zprior_2.GetStrongLinePresenceLogZPrior(
      linePresence, penalization_factor);
  BOOST_CHECK(logzPrior_2 == logzPrior);
}

BOOST_AUTO_TEST_CASE(CombineLogZPrior_test) {
  // GetConstantLogZPrior
  TFloat64List redshifts = genRedshifts(5);
  CZPrior zprior = CZPrior(true, redshifts);
  TFloat64List logzPrior;
  logzPrior = zprior.GetConstantLogZPrior(redshifts.size());

  // size of logzPrior_2 different from size of logzPrior_1
  TFloat64List logzPrior_2 = {0, 0, 0, 0};
  BOOST_CHECK_THROW(zprior.CombineLogZPrior(logzPrior, logzPrior_2),
                    GlobalException);

  // size of logzPrior_2 equal to size of logzPrior_1
  logzPrior_2 = zprior.GetConstantLogZPrior(redshifts.size());
  TFloat64List logzPriorCombined =
      zprior.CombineLogZPrior(logzPrior, logzPrior_2);
  BOOST_CHECK(logzPriorCombined == logzPrior &&
              logzPriorCombined == logzPrior_2);
}

BOOST_AUTO_TEST_CASE(GetNLinesSNRAboveCutLogZPrior_test) {
  CZPrior zprior;
  TInt32List nlinesAboveSNR = {0, 0, 0, 5, 5, 5, 0, 0, 0};
  Float64 penalization_factor = 0.5;

  // GetNLinesSNRAboveCutLogZPrior_test
  // without normalization
  TFloat64List logzPrior =
      zprior.GetNLinesSNRAboveCutLogZPrior(nlinesAboveSNR, penalization_factor);
  for (Int32 i = 0; i < logzPrior.size(); i++) {
    if (i == 3 || i == 4 || i == 5)
      BOOST_CHECK(logzPrior[i] == 0.);
    else
      BOOST_CHECK_CLOSE(logzPrior[i], log(penalization_factor), precision);
  }

  // GetNLinesSNRAboveCutLogZPrior_test
  // with normalization
  TFloat64List redshifts = genRedshifts(9);
  CZPrior zprior_ref = CZPrior(true, redshifts);
  zprior_ref.NormalizePrior(logzPrior);

  CZPrior zprior_2 = CZPrior(true, redshifts);
  TFloat64List logzPrior_2 = zprior_2.GetNLinesSNRAboveCutLogZPrior(
      nlinesAboveSNR, penalization_factor);
  BOOST_CHECK(logzPrior_2 == logzPrior);
}

BOOST_AUTO_TEST_CASE(GetEuclidNhaLogZPrior_test) {
  TFloat64List redshifts = genRedshifts(5);
  CZPrior zprior(true, redshifts);
  TFloat64List logzprior_ref = {-26.325572589439748, 0.6646422550761828,
                                1.0460496580236389, 1.2134721088424685,
                                1.3054485264887992};

  Float64 aCoeff = 0.;
  BOOST_CHECK_THROW(zprior.GetEuclidNhaLogZPrior(redshifts, aCoeff),
                    GlobalException);

  aCoeff = 0.5;
  TFloat64List logzPrior = zprior.GetEuclidNhaLogZPrior(redshifts, aCoeff);
  for (Int32 i = 1; i < logzPrior.size(); i++) {
    BOOST_CHECK_CLOSE(logzPrior[i], logzprior_ref[i], precision);
  }
}

BOOST_AUTO_TEST_SUITE_END()