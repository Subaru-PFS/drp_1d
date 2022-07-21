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
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/statistics/pdfcandidatesz.h"
#include "RedshiftLibrary/statistics/pdfcandidateszresult.h"
#include "RedshiftLibrary/statistics/zprior.h"

#include <boost/test/unit_test.hpp>
using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(Statistics_pdfcandidatesz)

TFloat64Range redshiftRange = TFloat64Range(0, 5);
Float64 redshiftStep = 1E-4;
TRedshiftList pdfz = redshiftRange.SpreadOverLogZplusOne(redshiftStep);
Float64 precision = 1e-10;

// Gaussian definition
TFloat64List generateGaussian(TFloat64List z, Float64 sigma, Float64 mu) {
  TFloat64List gaussY(z.size());
  for (Int32 i = 0; i < z.size(); i++) {
    gaussY[i] = -log(sigma * std::sqrt(2 * M_PI)) -
                0.5 * ((z[i] - mu) / sigma) * ((z[i] - mu) / sigma);
  }
  return gaussY;
}

TFloat64Range gaussRange = TFloat64Range(1e-4, 2.);
TRedshiftList gaussZ = gaussRange.SpreadOver(redshiftStep);
Float64 sigma = 0.1;
Float64 mu = 1.;
TFloat64List gaussY = generateGaussian(gaussZ, sigma, mu);

// Constructor
BOOST_AUTO_TEST_CASE(TCandidateZ_test) {

  // Default
  TCandidateZ candidate;
  BOOST_CHECK(candidate.ValSumProba == 0.);

  // copy
  candidate.ValSumProba = 2.1;
  TCandidateZ candidate_2(candidate);
  BOOST_CHECK(candidate_2.ValSumProba == 2.1);

  // copy assignement
  TCandidateZ candidate_3;
  candidate_3 = candidate;
  BOOST_CHECK(candidate_3.ValSumProba == 2.1);

  // move
  TCandidateZ candidate_4(std::move(candidate));
  BOOST_CHECK(candidate_4.ValSumProba == 2.1);

  // move assignement
  TCandidateZ candidate_5;
  candidate_5 = std::move(candidate_2);
  BOOST_CHECK(candidate_5.ValSumProba == 2.1);
}

// Constructor
BOOST_AUTO_TEST_CASE(CPdfCandidatesZ_test) {

  // const with redshifts list
  TRedshiftList redshifts = {1.0, 2.0};
  CPdfCandidatesZ zcand_redshift = CPdfCandidatesZ(redshifts);

  BOOST_CHECK(zcand_redshift.m_candidates.size() == 2);
  BOOST_CHECK(zcand_redshift.m_candidates["EXT0"]->Redshift == 1.0);
  BOOST_CHECK(zcand_redshift.m_candidates["EXT1"]->Redshift == 2.0);

  // const with candidates list
  TCandidateZbyID zcandidates = {{"EXT0", std::make_shared<TCandidateZ>()},
                                 {"EXT1", std::make_shared<TCandidateZ>()}};
  zcandidates["EXT0"]->Redshift = 1.0;
  zcandidates["EXT1"]->Redshift = 2.0;
  CPdfCandidatesZ zcand_candList = CPdfCandidatesZ(zcandidates);

  BOOST_CHECK(zcand_candList.m_candidates.size() == 2);
  BOOST_CHECK(zcand_candList.m_candidates["EXT0"]->Redshift == 1.0);
  BOOST_CHECK(zcand_candList.m_candidates["EXT1"]->Redshift == 2.0);
}

// both redshifts belong to overlapping range
BOOST_AUTO_TEST_CASE(Deltaz_overlapping1) {
  TRedshiftList center_redshifts = {1.0, 1.5};
  TRedshiftList deltaz = {0.5 / 3, 0.5 / 3};
  TCandidateZRangebyID ranges;
  TCandidateZRangebyID correct_ranges = {{"EXT0", {0.5, 1.24990}},
                                         {"EXT1", {1.25, 2}}};

  TCandidateZbyID zcandidates = {{"EXT0", std::make_shared<TCandidateZ>()},
                                 {"EXT1", std::make_shared<TCandidateZ>()}};
  zcandidates["EXT0"]->Redshift = center_redshifts[0];
  zcandidates["EXT0"]->Deltaz = deltaz[0];
  zcandidates["EXT1"]->Redshift = center_redshifts[1];
  zcandidates["EXT1"]->Deltaz = deltaz[1];
  CPdfCandidatesZ zcand_op = CPdfCandidatesZ(zcandidates);
  zcand_op.SetIntegrationWindows(TFloat64Range(pdfz), ranges);

  BOOST_CHECK_CLOSE(ranges["EXT0"].GetEnd(), correct_ranges["EXT0"].GetEnd(),
                    1E-4);
  BOOST_CHECK_CLOSE(ranges["EXT1"].GetEnd(), correct_ranges["EXT1"].GetEnd(),
                    1E-4);
  BOOST_CHECK_CLOSE(ranges["EXT0"].GetBegin(),
                    correct_ranges["EXT0"].GetBegin(), 1E-4);
  BOOST_CHECK_CLOSE(ranges["EXT1"].GetBegin(),
                    correct_ranges["EXT1"].GetBegin(), 1E-4);

  // deltaz no set : m_dzDefault * (1 + cand.Redshift)
  zcandidates["EXT0"]->Deltaz = -1;
  CPdfCandidatesZ zcand_op_2 = CPdfCandidatesZ(zcandidates);
  zcand_op_2.SetIntegrationWindows(TFloat64Range(pdfz), ranges);
  BOOST_CHECK_CLOSE(zcand_op_2.m_candidates["EXT0"]->Deltaz,
                    1e-3 * (1 + zcandidates["EXT0"]->Redshift), precision);
}

// no overlapping
BOOST_AUTO_TEST_CASE(Deltaz_nooverlapping) {
  TRedshiftList center_redshifts = {1.0, 4.0};
  TRedshiftList deltaz = {0.5 / 3, 0.5 / 3};
  TCandidateZRangebyID ranges;
  TCandidateZRangebyID correct_ranges = {{"EXT0", {0.5, 1.5}},
                                         {"EXT1", {3.5, 4.5}}};

  TCandidateZbyID zcandidates = {{"EXT0", std::make_shared<TCandidateZ>()},
                                 {"EXT1", std::make_shared<TCandidateZ>()}};
  zcandidates["EXT0"]->Redshift = center_redshifts[0];
  zcandidates["EXT0"]->Deltaz = deltaz[0];
  zcandidates["EXT1"]->Redshift = center_redshifts[1];
  zcandidates["EXT1"]->Deltaz = deltaz[1];
  CPdfCandidatesZ zcand_op = CPdfCandidatesZ(zcandidates);
  zcand_op.SetIntegrationWindows(TFloat64Range(pdfz), ranges);

  BOOST_CHECK_CLOSE(ranges["EXT0"].GetEnd(), correct_ranges["EXT0"].GetEnd(),
                    1E-4);
  BOOST_CHECK_CLOSE(ranges["EXT1"].GetEnd(), correct_ranges["EXT1"].GetEnd(),
                    1E-4);
  BOOST_CHECK_CLOSE(ranges["EXT0"].GetBegin(),
                    correct_ranges["EXT0"].GetBegin(), 1E-4);
  BOOST_CHECK_CLOSE(ranges["EXT1"].GetBegin(),
                    correct_ranges["EXT1"].GetBegin(), 1E-4);
}
// only the smallest redshift belongs to the overlapping region
BOOST_AUTO_TEST_CASE(Deltaz_overlapping_3) {
  TRedshiftList center_redshifts = {1.0, 4.0};
  TRedshiftList deltaz = {1. / 3, 2.5 / 3};
  TCandidateZRangebyID ranges;
  TCandidateZRangebyID correct_ranges = {{"EXT0", {0., 1.7499}},
                                         {"EXT1", {1.75, pdfz.back()}}};

  TCandidateZbyID zcandidates = {{"EXT0", std::make_shared<TCandidateZ>()},
                                 {"EXT1", std::make_shared<TCandidateZ>()}};
  zcandidates["EXT0"]->Redshift = center_redshifts[0];
  zcandidates["EXT0"]->Deltaz = deltaz[0];
  zcandidates["EXT1"]->Redshift = center_redshifts[1];
  zcandidates["EXT1"]->Deltaz = deltaz[1];
  CPdfCandidatesZ zcand_op = CPdfCandidatesZ(zcandidates);
  zcand_op.SetIntegrationWindows(TFloat64Range(pdfz), ranges);

  BOOST_CHECK_CLOSE(ranges["EXT0"].GetEnd(), correct_ranges["EXT0"].GetEnd(),
                    1E-4);
  BOOST_CHECK_CLOSE(ranges["EXT1"].GetEnd(), correct_ranges["EXT1"].GetEnd(),
                    1E-4);
  BOOST_CHECK_CLOSE(ranges["EXT0"].GetBegin(),
                    correct_ranges["EXT0"].GetBegin(), 1E-4);
  BOOST_CHECK_CLOSE(ranges["EXT1"].GetBegin(),
                    correct_ranges["EXT1"].GetBegin(), 1E-4);
}

// only the greatest redshift belong to overlapping region
BOOST_AUTO_TEST_CASE(Deltaz_overlapping_4) {
  TRedshiftList center_redshifts = {1.0, 4.0};
  TRedshiftList deltaz = {2.0 / 3, 2.5 / 3};
  TCandidateZRangebyID ranges;
  TCandidateZRangebyID correct_ranges = {{"EXT0", {0., 2.2499}},
                                         {"EXT1", {2.25, pdfz.back()}}};

  TCandidateZbyID zcandidates = {{"EXT0", std::make_shared<TCandidateZ>()},
                                 {"EXT1", std::make_shared<TCandidateZ>()}};
  zcandidates["EXT0"]->Redshift = center_redshifts[0];
  zcandidates["EXT0"]->Deltaz = deltaz[0];
  zcandidates["EXT1"]->Redshift = center_redshifts[1];
  zcandidates["EXT1"]->Deltaz = deltaz[1];
  CPdfCandidatesZ zcand_op = CPdfCandidatesZ(zcandidates);
  zcand_op.SetIntegrationWindows(TFloat64Range(pdfz), ranges);

  BOOST_CHECK_CLOSE(ranges["EXT0"].GetEnd(), correct_ranges["EXT0"].GetEnd(),
                    1E-4);
  BOOST_CHECK_CLOSE(ranges["EXT1"].GetEnd(), correct_ranges["EXT1"].GetEnd(),
                    1E-4);
  BOOST_CHECK_CLOSE(ranges["EXT0"].GetBegin(),
                    correct_ranges["EXT0"].GetBegin(), 1E-4);
  BOOST_CHECK_CLOSE(ranges["EXT1"].GetBegin(),
                    correct_ranges["EXT1"].GetBegin(), 1E-4);
}

// candidates too closed
BOOST_AUTO_TEST_CASE(Deltaz_overlapping_5) {
  TRedshiftList center_redshifts = {1.0, 1.00005};
  TRedshiftList deltaz = {2.0 / 3, 2.5 / 3};
  TCandidateZRangebyID ranges;
  TCandidateZRangebyID correct_ranges = {{"EXT1", {pdfz.front(), 3.50005}}};

  TCandidateZbyID zcandidates = {{"EXT0", std::make_shared<TCandidateZ>()},
                                 {"EXT1", std::make_shared<TCandidateZ>()}};
  zcandidates["EXT0"]->Redshift = center_redshifts[0];
  zcandidates["EXT0"]->Deltaz = deltaz[0];
  zcandidates["EXT1"]->Redshift = center_redshifts[1];
  zcandidates["EXT1"]->Deltaz = deltaz[1];
  CPdfCandidatesZ zcand_op = CPdfCandidatesZ(zcandidates);
  TStringList b = zcand_op.SetIntegrationWindows(TFloat64Range(pdfz), ranges);

  BOOST_CHECK(b[0] == "EXT0");
  BOOST_CHECK_CLOSE(ranges["EXT1"].GetEnd(), correct_ranges["EXT1"].GetEnd(),
                    1E-4);
  BOOST_CHECK_CLOSE(ranges["EXT1"].GetBegin(),
                    correct_ranges["EXT1"].GetBegin(), 1E-4);
}

// no range found including the candidate
BOOST_AUTO_TEST_CASE(Deltaz_overlapping_6) {
  TRedshiftList center_redshifts = {5.5};
  TRedshiftList deltaz = {2.5 / 3};
  TCandidateZRangebyID ranges;

  TCandidateZbyID zcandidates = {{"EXT0", std::make_shared<TCandidateZ>()}};
  zcandidates["EXT0"]->Redshift = center_redshifts[0];
  zcandidates["EXT0"]->Deltaz = deltaz[0];
  CPdfCandidatesZ zcand_op = CPdfCandidatesZ(zcandidates);

  BOOST_CHECK_THROW(zcand_op.SetIntegrationWindows(TFloat64Range(pdfz), ranges),
                    GlobalException);
}

BOOST_AUTO_TEST_CASE(SortByValSumProbaInt_test) {
  TCandidateZbyID zcandidates = {{"EXT0", std::make_shared<TCandidateZ>()},
                                 {"EXT1", std::make_shared<TCandidateZ>()},
                                 {"EXT2", std::make_shared<TCandidateZ>()},
                                 {"EXT3", std::make_shared<TCandidateZ>()}};
  zcandidates["EXT0"]->ValSumProba = 0.;
  zcandidates["EXT1"]->ValSumProba = 1.5;
  zcandidates["EXT2"]->ValSumProba = 1.;
  TCandidateZbyRank ranked_candidates;

  CPdfCandidatesZ zcand_op = CPdfCandidatesZ(zcandidates);
  zcand_op.SortByValSumProbaInt(ranked_candidates);

  // Check first is OK :
  BOOST_CHECK(ranked_candidates[0].first == "EXT1");
  BOOST_CHECK(ranked_candidates[1].first == "EXT2");
  BOOST_CHECK(ranked_candidates[2].first == "EXT0");

  // Check second valSumProba is OK :
  BOOST_CHECK(ranked_candidates[0].second->ValSumProba ==
              zcandidates["EXT1"]->ValSumProba);
  BOOST_CHECK(ranked_candidates[1].second->ValSumProba ==
              zcandidates["EXT2"]->ValSumProba);
  BOOST_CHECK(ranked_candidates[2].second->ValSumProba ==
              zcandidates["EXT0"]->ValSumProba);

  // Check order is preserved :
  zcandidates["EXT0"]->ValSumProba = 0.;
  zcandidates["EXT1"]->ValSumProba = 0.;
  zcandidates["EXT2"]->ValSumProba = 0.;
  zcandidates["EXT3"]->ValSumProba = 0.;
  CPdfCandidatesZ zcand_op_2 = CPdfCandidatesZ(zcandidates);
  zcand_op_2.SortByValSumProbaInt(ranked_candidates);

  BOOST_CHECK(ranked_candidates[0].first == "EXT0");
  BOOST_CHECK(ranked_candidates[1].first == "EXT1");
  BOOST_CHECK(ranked_candidates[2].first == "EXT2");
  BOOST_CHECK(ranked_candidates[3].first == "EXT3");

  // Test with one candidate
  TCandidateZbyID zcandidates_alone = {
      {"EXT0", std::make_shared<TCandidateZ>()}};
  zcandidates_alone["EXT0"]->ValSumProba = 0.;
  CPdfCandidatesZ zcand_alone = CPdfCandidatesZ(zcandidates_alone);
  zcand_alone.SortByValSumProbaInt(ranked_candidates);
  BOOST_CHECK(ranked_candidates[0].second->ValSumProba ==
              zcandidates["EXT0"]->ValSumProba);
}

BOOST_AUTO_TEST_CASE(getCandidateSumTrapez_test) {
  TFloat64List valprobalog = {0.2, 0.8, 0.2};
  TFloat64Range zrange(2.1, 5.1);
  std::shared_ptr<TCandidateZ> candidate = std::make_shared<TCandidateZ>();
  candidate->Redshift = 1.;

  // redshifts not sorted
  TRedshiftList center_redshifts = {5.0, 2.0, 3.0};
  CPdfCandidatesZ zcand_op = CPdfCandidatesZ(center_redshifts);
  BOOST_CHECK_THROW(zcand_op.getCandidateSumTrapez(
                        center_redshifts, valprobalog, zrange, candidate),
                    GlobalException);

  // no enclosing interval
  center_redshifts = {2.0, 3.0, 5.0};
  BOOST_CHECK_THROW(zcand_op.getCandidateSumTrapez(
                        center_redshifts, valprobalog, zrange, candidate),
                    GlobalException);

  // Test OK
  TFloat64Range zrange_2(gaussZ.front() + 1e-6, gaussZ.back() - 1e-6);
  CPdfCandidatesZ zcand_op_2 = CPdfCandidatesZ(gaussZ);
  bool ret =
      zcand_op_2.getCandidateSumTrapez(gaussZ, gaussY, zrange_2, candidate);
  BOOST_CHECK(ret == 1);
  BOOST_CHECK_CLOSE(candidate->ValSumProbaZmin, gaussZ.front(), precision);
  BOOST_CHECK_CLOSE(candidate->ValSumProbaZmax, gaussZ.back(), precision);

  Float64 logSum = COperatorPdfz::logSumExpTrick(gaussY, gaussZ);
  BOOST_CHECK_CLOSE(candidate->ValSumProba, exp(logSum), precision);
  BOOST_CHECK_CLOSE(candidate->ValSumProba, 1.0, precision);
}

BOOST_AUTO_TEST_CASE(getCandidateRobustGaussFit_test) {
  std::shared_ptr<TCandidateZ> candidate = std::make_shared<TCandidateZ>();
  candidate->Redshift = 1.;

  TFloat64Range zrange_2(gaussZ.front() + 1e-6, gaussZ.back() - 1e-6);
  CPdfCandidatesZ zcand_op = CPdfCandidatesZ(gaussZ);
  bool fitSuccessful =
      zcand_op.getCandidateRobustGaussFit(gaussZ, gaussY, zrange_2, candidate);

  BOOST_CHECK_CLOSE(candidate->GaussAmp * candidate->GaussSigma *
                        sqrt(2 * M_PI),
                    1.0, precision);
  BOOST_CHECK_CLOSE(candidate->GaussSigma, sigma, precision);
  BOOST_CHECK_CLOSE(candidate->GaussAmp, 1 / (sigma * std::sqrt(2 * M_PI)),
                    precision);
}

BOOST_AUTO_TEST_CASE(Compute_test) {
  TCandidateZbyID candidate_ref = {{"EXT0", std::make_shared<TCandidateZ>()}};
  candidate_ref["EXT0"]->Redshift = 1.;

  TFloat64List gaussY_2 = generateGaussian(gaussZ, 0.3, mu);
  CPdfCandidatesZ zcand_op = CPdfCandidatesZ(candidate_ref);

  // compute with direct integration
  std::shared_ptr<PdfCandidatesZResult> result =
      zcand_op.Compute(gaussZ, gaussY_2);

  TCandidateZbyRank candidate = result->m_ranked_candidates;

  BOOST_CHECK_CLOSE(candidate[0].second->ValSumProbaZmin, gaussZ.front(),
                    precision);
  BOOST_CHECK_CLOSE(candidate[0].second->ValSumProbaZmax, gaussZ.back(),
                    precision);

  Float64 logSum = COperatorPdfz::logSumExpTrick(gaussY_2, gaussZ);
  BOOST_CHECK_CLOSE(candidate[0].second->ValSumProba, exp(logSum), precision);

  // compute with gaussian fit
  zcand_op.m_optMethod = 1;
  result = zcand_op.Compute(gaussZ, gaussY_2);
  candidate = result->m_ranked_candidates;

  BOOST_CHECK_CLOSE(candidate[0].second->GaussAmp *
                        candidate[0].second->GaussSigma * sqrt(2 * M_PI),
                    1.0, precision);
  BOOST_CHECK_CLOSE(candidate[0].second->GaussSigma, 0.3, precision);
  BOOST_CHECK_CLOSE(candidate[0].second->GaussAmp,
                    1 / (0.3 * std::sqrt(2 * M_PI)), precision);
}
BOOST_AUTO_TEST_SUITE_END()
