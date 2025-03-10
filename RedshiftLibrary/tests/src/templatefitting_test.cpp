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

#include "RedshiftLibrary/operator/powerlaw.h"
#include "RedshiftLibrary/processflow/context.h"
#include "tests/src/templatefittingfortests.h"
#include "tests/src/tool/inputContextLight.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(templateFitting_test)

BOOST_FIXTURE_TEST_CASE(fitQuality_test,
                        fixture_TemplateFittingSolveTestNoFFT) {

  // Prepares contexte for BasicFit
  CAutoScope spectrumModel_autoscope(Context.m_ScopeStack, "galaxy",
                                     ScopeType::SPECTRUMMODEL);
  CAutoScope stage_autoscope(Context.m_ScopeStack, "redshiftSolver",
                             ScopeType::STAGE);
  CTemplateFittingSolve templateFittingSolve;

  auto const &inputContext = *Context.GetInputContext();
  templateFittingSolve.InitRanges(inputContext);

  auto templateFittingOperator =
      COperatorTemplateFitting(templateFittingSolve.m_redshifts);

  // Prepares arguments for BasicFit
  std::shared_ptr<const CTemplate> tpl = fixture_SharedGalaxyTemplate().tpl;
  Float64 redshift = templateFittingSolve.m_redshifts[0];
  Float64 overlapThreshold = 1;
  bool opt_extinction = false;
  bool opt_dustFitting = false;
  CPriorHelper::TPriorEList logpriore = CPriorHelper::TPriorEList();
  auto MeiksinList = {0};
  auto EbmvList = {0};

  auto result = templateFittingOperator.BasicFit(
      tpl, redshift, overlapThreshold, opt_extinction, opt_dustFitting,
      logpriore, MeiksinList, EbmvList);

  // Checks that fit quality values are set
  BOOST_CHECK_CLOSE(result.chiSquare, 335.74708870710344, 1e-6);
  BOOST_CHECK_CLOSE(result.reducedChiSquare, 6.3348507303227066, 1e-6);
  BOOST_CHECK_CLOSE(result.pValue, 3.9524304409427213e-43, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()