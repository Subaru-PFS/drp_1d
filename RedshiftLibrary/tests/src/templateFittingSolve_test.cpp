
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
#include <boost/test/unit_test.hpp>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/method/templatefittingsolve.h"
#include "RedshiftLibrary/method/templatefittingsolveresult.h"
#include "RedshiftLibrary/operator/extremaresult.h"
#include "RedshiftLibrary/operator/templatefittinglog.h"
#include "RedshiftLibrary/processflow/context.h"
#include "tests/src/templatefittingfortests.h"
#include "tests/src/tool/inputContextLight.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(templateFittingSolve_test)

BOOST_FIXTURE_TEST_CASE(computeNoFFT_test,
                        fixture_TemplateFittingSolveTestNoFFT) {

  CAutoScope spectrumModel_autoscope(Context.m_ScopeStack, "galaxy",
                                     ScopeType::SPECTRUMMODEL);
  CAutoScope stage_autoscope(Context.m_ScopeStack, "redshiftSolver",
                             ScopeType::STAGE);
  CTemplateFittingSolve templateFittingSolve;
  BOOST_REQUIRE_NO_THROW(templateFittingSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("galaxy", "redshiftSolver",
                                               "templateFittingSolve");
  BOOST_CHECK(result_out.lock()->getType() == "CTemplateFittingSolveResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "redshiftSolver", "templateFittingSolve", "pdf");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "redshiftSolver", "templateFittingSolve", "pdf_params");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  std::string resType = Context.GetResultStore()->GetCandidateResultType(
      "galaxy", "redshiftSolver", "templateFittingSolve", "extrema_results",
      "model_parameters");
  BOOST_CHECK(resType == "TExtremaResult");

  std::shared_ptr<const TExtremaResult> res =
      Context.GetResultStore()->GetExtremaResult(
          "galaxy", "redshiftSolver", "templateFittingSolve", "extrema_results",
          "model_parameters", 0);
  Float64 z = res->Redshift;
  BOOST_CHECK_CLOSE(z, 2.8770415147926256, 1e-6);

  // Checks that fit quality indicators are correctly set
  Float64 chi2 = res->fittedContinuum.merit;
  BOOST_CHECK_CLOSE(chi2, 366.77585884307723, 1e-4);
  Float64 chi2r = res->fittedContinuum.reducedChi2;
  BOOST_CHECK_CLOSE(chi2r, 6.6686519789650403, 1e-4);
  Float64 pValue = res->fittedContinuum.pValue;
  BOOST_CHECK_CLOSE(pValue, 4.6089808815878171e-48, 1e-4);

  Context.reset();
}

BOOST_FIXTURE_TEST_CASE(compute2Pass_test,
                        fixture_TemplateFittingSolve2PassTest) {
  CAutoScope spectrumModel_autoscope(Context.m_ScopeStack, "galaxy",
                                     ScopeType::SPECTRUMMODEL);
  CAutoScope stage_autoscope(Context.m_ScopeStack, "redshiftSolver",
                             ScopeType::STAGE);

  CTemplateFittingSolve templateFittingSolve;
  BOOST_REQUIRE_NO_THROW(templateFittingSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("galaxy", "redshiftSolver",
                                               "templateFittingSolve");
  BOOST_CHECK(result_out.lock()->getType() == "CTemplateFittingSolveResult");

  // First pass pdf results
  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "redshiftSolver", "templateFittingSolve", "firstpass_pdf");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");
  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "redshiftSolver", "templateFittingSolve",
      "firstpass_pdf_params");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "redshiftSolver", "templateFittingSolve", "pdf");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "redshiftSolver", "templateFittingSolve", "pdf_params");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  std::string resType = Context.GetResultStore()->GetCandidateResultType(
      "galaxy", "redshiftSolver", "templateFittingSolve", "extrema_results",
      "model_parameters");
  BOOST_CHECK(resType == "TExtremaResult");

  std::shared_ptr<const TExtremaResult> res =
      Context.GetResultStore()->GetExtremaResult(
          "galaxy", "redshiftSolver", "templateFittingSolve", "extrema_results",
          "model_parameters", 0);
  Float64 z = res->Redshift;

  // Checks that fit quality indicators are correctly set
  Float64 chi2 = res->fittedContinuum.merit;
  BOOST_CHECK_CLOSE(chi2, 332.39827811512555, 1e-4);
  Float64 chi2r = res->fittedContinuum.reducedChi2;
  BOOST_CHECK_CLOSE(chi2r, 6.2716656248136893, 1e-4);
  Float64 pValue = res->fittedContinuum.pValue;
  BOOST_CHECK_CLOSE(pValue, 1.6442425635585607e-42, 1e-4);
}

BOOST_FIXTURE_TEST_CASE(computeFFT_test, fixture_TemplateFittingSolveTestFFT) {
  CAutoScope spectrumModel_autoscope(Context.m_ScopeStack, "galaxy",
                                     ScopeType::SPECTRUMMODEL);
  CAutoScope stage_autoscope(Context.m_ScopeStack, "redshiftSolver",
                             ScopeType::STAGE);

  CTemplateFittingSolve templateFittingSolve;
  BOOST_REQUIRE_NO_THROW(templateFittingSolve.Compute());

  std::weak_ptr<const COperatorResult> result_out =
      Context.GetResultStore()->GetSolveResult("galaxy", "redshiftSolver",
                                               "templateFittingSolve");
  BOOST_CHECK(result_out.lock()->getType() == "CTemplateFittingSolveResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "redshiftSolver", "templateFittingSolve", "pdf");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  result_out = Context.GetResultStore()->GetLogZPdfResult(
      "galaxy", "redshiftSolver", "templateFittingSolve", "pdf_params");
  BOOST_CHECK(result_out.lock()->getType() == "CLogZPdfResult");

  std::string resType = Context.GetResultStore()->GetCandidateResultType(
      "galaxy", "redshiftSolver", "templateFittingSolve", "extrema_results",
      "model_parameters");
  BOOST_CHECK(resType == "TExtremaResult");

  std::shared_ptr<const TExtremaResult> res =
      Context.GetResultStore()->GetExtremaResult(
          "galaxy", "redshiftSolver", "templateFittingSolve", "extrema_results",
          "model_parameters", 0);

  Float64 z = res->Redshift;
  BOOST_CHECK_CLOSE(z, 2.880219830862035, 1e-6);

  // Checks that fit quality indicators are correctly set
  Float64 chi2 = res->fittedContinuum.merit;
  BOOST_CHECK_CLOSE(chi2, 331.98611446809548, 1e-4);
  Float64 chi2r = res->fittedContinuum.reducedChi2;
  BOOST_CHECK_CLOSE(chi2r, 6.1478910086684344, 1e-4);
  Float64 pValue = res->fittedContinuum.pValue;
  BOOST_CHECK_CLOSE(pValue, 4.9918699914782802e-42, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
