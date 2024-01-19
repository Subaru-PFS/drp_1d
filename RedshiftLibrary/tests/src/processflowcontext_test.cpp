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
#include "RedshiftLibrary/method/templatefittingsolve.h"
#include "RedshiftLibrary/method/templatefittingsolveresult.h"
#include "RedshiftLibrary/processflow/context.h"
#include "tests/src/tool/inputContextLight.h"
#include <boost/test/unit_test.hpp>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_multifit_nlin.h>

using namespace NSEpic;

const std::string jsonString =
    "{\"lambdaRange\" : [ 4630, 4815 ],"
    "\"smoothWidth\" : 0.0,"
    "\"templateCatalog\" : {"
    "\"continuumRemoval\" : {"
    "\"method\" : \"zero\","
    "\"medianKernelWidth\" : 75,"
    "\"medianEvenReflection\" : true}},"
    "\"ebmv\" : {\"start\" : 0, \"step\" : 0.1, \"count\" : 10},"
    "\"continuumRemoval\" : {"
    "\"method\" : \"irregularSamplingMedian\","
    "\"medianKernelWidth\" : 400,"
    "\"medianEvenReflection\" : true,"
    "\"decompScales\" : 9},"
    "\"lsf\" : {\"lsfType\" : \"gaussianConstantResolution\", \"resolution\" : "
    "4300},"
    "\"extremaRedshiftSeparation\" : 0.01,"
    "\"spectrumModels\" : [\"galaxy\"],"
    "\"autoCorrectInput\" : false,"
    "\"galaxy\" : {"
    "\"stages\": [\"redshiftSolver\"],"
    "\"redshiftRange\" : [ 2.84, 2.88 ],"
    "\"redshiftStep\" : 0.0001,"
    "\"redshiftSampling\" : \"log\","
    "\"templateDir\" : \"templates/BC03_sdss_tremonti21\","
    "\"redshiftSolver\": {"
    "\"method\" : \"templateFittingSolve\","
    "\"linemeas_method\" : null,"
    "\"lineModelSolve\" : {"
    "\"lineModel\" : {"
    "\"lineTypeFilter\" : \"no\","
    "\"lineForceFilter\" : \"no\"}},"
    "\"templateFittingSolve\" : {"
    "\"extremaCount\" : 5,"
    "\"overlapThreshold\" : 1,"
    "\"spectrum\" : {\"component\" : \"raw\"},"
    "\"fftProcessing\" : true,"
    "\"interpolation\" : \"preComputedFineGrid\","
    "\"igmFit\" : true,"
    "\"ismFit\" : true,"
    "\"pdfCombination\" : \"marg\","
    "\"enablePhotometry\" : false}},"
    "\"airVacuumMethod\" : \"default\"}}";

class fixture_processflowcontextTest {
public:
  TScopeStack scopeStack;
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin =
      fixture_MeiskinCorrection().igmCorrectionMeiksin;
  std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti =
      fixture_CalzettiCorrection().ismCorrectionCalzetti;
  std::shared_ptr<CLSF> LSF =
      fixture_LSFGaussianConstantResolution(scopeStack).LSF;
  std::shared_ptr<CSpectrum> spc = fixture_SharedSpectrumExtended().spc;
  std::shared_ptr<CTemplateCatalog> catalog =
      fixture_sharedTemplateCatalog().catalog;
  std::shared_ptr<CPhotBandCatalog> photoBandCatalog =
      fixture_PhotoBandCatalog().photoBandCatalog;
  std::shared_ptr<CLineCatalog> lineCatalog = fixture_LineCatalog().lineCatalog;
  std::shared_ptr<CLineCatalogsTplRatio> lineRatioTplCatalog =
      fixture_LineRatioTplCatalog().lineRatioTplCatalog;
};

BOOST_FIXTURE_TEST_SUITE(processflowcontext_test,
                         fixture_processflowcontextTest)

BOOST_AUTO_TEST_CASE(context_test) {

  Context.LoadParameterStore(jsonString);
  std::shared_ptr<const CParameterStore> paramStore =
      Context.GetParameterStore();
  BOOST_CHECK(paramStore->Get<TFloat64Range>("lambdaRange") ==
              TFloat64Range(4630, 4815));

  std::shared_ptr<const CInputContext> inputCtx = Context.GetInputContext();
  BOOST_CHECK(inputCtx->getSpectra().size() == 0);

  std::shared_ptr<COperatorResultStore> resultStore = Context.GetResultStore();
  BOOST_CHECK(resultStore->HasDataset("galaxy", "redshiftSolver", "templateFittingSolve",
                                      "solveResult") == false);

  spc->SetLSF(LSF);
  Context.addSpectrum(spc);
  BOOST_CHECK(Context.GetSpectrum() == spc);

  BOOST_CHECK(Context.getRebinnedSpectra().size() == 0);

  catalog->Add(fixture_SharedGalaxyTemplate().tpl);
  catalog->m_logsampling = 1;
  catalog->Add(fixture_SharedGalaxyTemplate().tpl);
  Context.setTemplateCatalog(catalog);
  BOOST_CHECK(Context.GetTemplateCatalog() == catalog);

  Context.setPhotBandCatalog(photoBandCatalog);
  BOOST_CHECK(Context.GetPhotBandCatalog() == photoBandCatalog);

  Context.setLineRatioCatalogCatalog("galaxy", lineRatioTplCatalog);
  Context.m_ScopeStack.push_back("galaxy");
  BOOST_CHECK(Context.GetTplRatioCatalog() == lineRatioTplCatalog);

  Context.setLineCatalog("galaxy", "lineModelSolve", lineCatalog);
  BOOST_CHECK(Context.GetLineCatalog("galaxy", "lineModelSolve") ==
              lineCatalog);
  Context.m_ScopeStack.push_back("redshiftSolver");
  Context.m_ScopeStack.push_back("lineModelSolve");
  BOOST_CHECK(Context.getCLineMap().size() ==
              fixture_LineCatalog().lineCatalogSize);

  Context.setfluxCorrectionCalzetti(ismCorrectionCalzetti);
  Context.setfluxCorrectionMeiksin(igmCorrectionMeiksin);

  Context.Init();
  BOOST_CHECK(Context.GetRebinnedSpectrum() = spc);

  std::shared_ptr<const TFloat64Range> lbdaRange =
      std::make_shared<const TFloat64Range>(4630, 4815);
  BOOST_CHECK(Context.GetLambdaRange()->GetBegin() == lbdaRange->GetBegin());
  BOOST_CHECK(Context.GetLambdaRange()->GetEnd() == lbdaRange->GetEnd());

  lbdaRange = std::make_shared<const TFloat64Range>(4630.478, 4815);
  BOOST_CHECK(Context.GetClampedLambdaRange(false)->GetBegin() ==
              lbdaRange->GetBegin());
  BOOST_CHECK(Context.GetClampedLambdaRange(false)->GetEnd() ==
              lbdaRange->GetEnd());

  BOOST_CHECK(Context.getRebinnedClampedLambdaRanges()[0]->GetBegin() ==
              Context.GetRebinnedSpectrum()
                  ->GetSpectralAxis()
                  .GetSamplesVector()
                  .front());
  BOOST_CHECK(Context.getRebinnedClampedLambdaRanges()[0]->GetEnd() ==
              Context.GetRebinnedSpectrum()
                  ->GetSpectralAxis()
                  .GetSamplesVector()
                  .back());

  Context.m_ScopeStack.pop_back();
  BOOST_CHECK_THROW(Context.GetCurrentMethod(), GlobalException);

  Context.m_ScopeStack.pop_back();
  Context.m_ScopeStack.pop_back();
  BOOST_CHECK_THROW(Context.GetCurrentCategory(), GlobalException);

  CTemplateFittingSolve templateFittingSolve(Context.m_ScopeStack, "galaxy");
  templateFittingSolve.Compute();

  BOOST_CHECK(resultStore->HasDataset("galaxy", "redshiftSolver", "templateFittingSolve",
                                      "solveResult") == true);
  Context.reset();
  BOOST_CHECK(Context.getSpectra().size() == 0);
  BOOST_CHECK(resultStore->HasDataset("galaxy", "redshiftSolver", "templateFittingSolve",
                                      "solveResult") == false);

  // GSL error
  Int32 dim = 1;
  gsl_matrix *covarMatrix = gsl_matrix_alloc(dim, dim);
  BOOST_CHECK_THROW(gsl_matrix_set(covarMatrix, 0, 1, 1.0), GlobalException);
  gsl_matrix_free(covarMatrix);
}

BOOST_AUTO_TEST_SUITE_END()