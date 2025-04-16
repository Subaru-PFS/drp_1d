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

#include "RedshiftLibrary/processflow/context.h"
#include "tests/src/templatefittingfortests.h"
#include "tests/src/tool/inputContextLight.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(templateFittingLog_test)

BOOST_FIXTURE_TEST_CASE(EstimateXtY_test, fixture_TemplateFittingSolveTestFFT) {
  CAutoScope spectrumModel_autoscope(Context.m_ScopeStack, "galaxy",
                                     ScopeType::SPECTRUMMODEL);
  CAutoScope stage_autoscope(Context.m_ScopeStack, "redshiftSolver",
                             ScopeType::STAGE);

  // Creation of useful objects
  CTemplateFittingSolve templateFittingSolve;
  BOOST_REQUIRE_NO_THROW(templateFittingSolve.Compute());

  Float64 precision = 1e-12;

  TFloat64Range lbdaR;
  TFloat64List XtY;
  TFloat64List XtYres;
  TFloat64List redshifts = {2.8399999999999999, 2.8404879658869557,
                            2.8409759937819086};
  COperatorTemplateFittingLog tplFittingLog(redshifts);

  // X size even, Y size even
  lbdaR.Set(1, 10);
  TFloat64List X = lbdaR.SpreadOver(1);
  lbdaR.Set(1, 14);
  TFloat64List Y = lbdaR.SpreadOver(1);

  EstimateXtYSlow(X, Y, XtY);

  tplFittingLog.m_nPaddedSamples = ceil(Y.size() / 2.0) * 2;
  tplFittingLog.InitFFT(tplFittingLog.m_nPaddedSamples);
  tplFittingLog.EstimateXtY(X, Y, XtYres, 0);

  for (std::size_t i = 0; i < XtY.size(); i++)
    BOOST_CHECK_CLOSE(XtY[i], XtYres[i], precision);

  // X size even, Y size odd
  lbdaR.Set(1, 15);
  Y = lbdaR.SpreadOver(1);

  EstimateXtYSlow(X, Y, XtY);

  tplFittingLog.m_nPaddedSamples = ceil(Y.size() / 2.0) * 2;
  tplFittingLog.InitFFT(tplFittingLog.m_nPaddedSamples);
  tplFittingLog.EstimateXtY(X, Y, XtYres, 0);

  for (std::size_t i = 0; i < XtY.size(); i++)
    BOOST_CHECK_CLOSE(XtY[i], XtYres[i], precision);

  // X size odd, Y size odd
  lbdaR.Set(1, 11);
  X = lbdaR.SpreadOver(1);

  EstimateXtYSlow(X, Y, XtY);

  tplFittingLog.m_nPaddedSamples = ceil(Y.size() / 2.0) * 2;
  tplFittingLog.InitFFT(tplFittingLog.m_nPaddedSamples);
  tplFittingLog.EstimateXtY(X, Y, XtYres, 0);

  for (std::size_t i = 0; i < XtY.size(); i++)
    BOOST_CHECK_CLOSE(XtY[i], XtYres[i], precision);

  // X size odd, Y size even
  lbdaR.Set(1, 14);
  Y = lbdaR.SpreadOver(1);

  EstimateXtYSlow(X, Y, XtY);

  tplFittingLog.m_nPaddedSamples = ceil(Y.size() / 2.0) * 2;
  tplFittingLog.InitFFT(tplFittingLog.m_nPaddedSamples);
  tplFittingLog.EstimateXtY(X, Y, XtYres, 0);

  for (std::size_t i = 0; i < XtY.size(); i++)
    BOOST_CHECK_CLOSE(XtY[i], XtYres[i], precision);
}

BOOST_AUTO_TEST_SUITE_END()
