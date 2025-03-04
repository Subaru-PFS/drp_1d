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
#ifndef TEMPLATEFITTING_TEST_H
#define TEMPLATEFITTING_TEST_H

#include <boost/test/unit_test.hpp>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/method/templatefittingsolve.h"
#include "RedshiftLibrary/method/templatefittingsolveresult.h"
#include "RedshiftLibrary/operator/extremaresult.h"
#include "RedshiftLibrary/operator/templatefittinglog.h"
#include "RedshiftLibrary/processflow/context.h"
#include "tests/src/tool/inputContextLight.h"

namespace NSEpic {

const std::string jsonString =
    "{\"multiObsMethod\" : \"\","
    "\"lambdaRange\" : [ 4680, 4712 ],"
    "\"smoothWidth\" : 0.0,"
    "\"nbSamplesMin\" : 1,"
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
    "\"airVacuumMethod\" : \"default\",";

const std::string jsonStringFFT = {
    "\"galaxy\" : {"
    "\"redshiftRange\" : [ 2.84, 2.88 ],"
    "\"redshiftStep\" : 0.0001,"
    "\"redshiftSampling\" : \"log\","
    "\"linemeas_method\" : null,"
    "\"redshiftSolver\": {"
    "\"method\" : \"templateFittingSolve\","
    "\"templateFittingSolve\" : {"
    "\"singlePass\" : true,"
    "\"extremaCount\" : 5,"
    "\"overlapThreshold\" : 1,"
    "\"spectrum\" : {\"component\" : \"raw\"},"
    "\"fftProcessing\" : true,"
    "\"interpolation\" : \"preComputedFineGrid\","
    "\"igmFit\" : true,"
    "\"ismFit\" : true,"
    "\"pdfCombination\" : \"marg\","
    "\"enablePhotometry\" : false}}}}"};

const std::string jsonStringNoFFT = {
    "\"galaxy\" : {"
    "\"redshiftRange\" : [ 2.84, 2.88 ],"
    "\"redshiftStep\" : 0.0001,"
    "\"redshiftSampling\" : \"log\","
    "\"linemeas_method\" : null,"
    "\"redshiftSolver\": {"
    "\"method\" : \"templateFittingSolve\","
    "\"templateFittingSolve\" : {"
    "\"singlePass\" : true,"
    "\"extremaCount\" : 5,"
    "\"overlapThreshold\" : 1,"
    "\"spectrum\" : {\"component\" : \"raw\"},"
    "\"fftProcessing\" : false,"
    "\"interpolation\" : \"preComputedFineGrid\","
    "\"igmFit\" : true,"
    "\"ismFit\" : true,"
    "\"pdfCombination\" : \"marg\","
    "\"enablePhotometry\" : true,"
    "\"photometry\": {\"weight\" : 1.0},"
    "\"extremaCutProbaThreshold\" : -1,"
    "\"firstPass\": {\"largeGridStepRatio\": 10, \"extremaCount\" : 5},"
    "\"secondPass\": {\"continuumFit\": \"fromFirstPass\", \"halfWindowSize\": "
    "0.001}"
    "}}}}"};
// Question: here on halfwindowsize : should it be < redshiftStep ?

class fixture_TemplateFittingCommon {
public:
  fixture_Context ctx;
  void Init(std::string fullJsonString) {
    fillCatalog();
    ctx.reset();
    ctx.loadParameterStore(fullJsonString);
    ctx.setCorrections(igmCorrectionMeiksin, ismCorrectionCalzetti);
    ctx.setCatalog(catalog);
    ctx.addSpectrum(spc, LSF);
    ctx.initContext();
    ctx.setPhotoBandCatalog(photoBandCatalog);
    spc->SetPhotData(photoData);
  }

  std::shared_ptr<CScopeStack> scopeStack = std::make_shared<CScopeStack>();
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin =
      fixture_MeiskinCorrection().igmCorrectionMeiksin;
  std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti =
      fixture_CalzettiCorrection().ismCorrectionCalzetti;
  std::shared_ptr<CLSF> LSF =
      fixture_LSFGaussianConstantResolution(scopeStack).LSF;
  std::shared_ptr<CSpectrum> spc = fixture_SharedSpectrum().spc;
  std::shared_ptr<CTemplateCatalog> catalog =
      fixture_sharedTemplateCatalog().catalog;
  std::shared_ptr<CPhotBandCatalog> photoBandCatalog =
      fixture_PhotoBandCatalog().photoBandCatalog;
  std::shared_ptr<CPhotometricData> photoData = fixture_PhotoData().photoData;

  void fillCatalog() {
    catalog->Add(fixture_SharedGalaxyTemplate().tpl);
    catalog->m_logsampling = 1;
    catalog->Add(fixture_SharedGalaxyTemplate().tpl);
  }
};

class fixture_TemplateFittingSolveTestNoFFT
    : public fixture_TemplateFittingCommon {
public:
  fixture_TemplateFittingSolveTestNoFFT() {
    Init(jsonString + jsonStringNoFFT);
  };
};

class fixture_TemplateFittingSolve2PassTest
    : public fixture_TemplateFittingCommon {
public:
  fixture_TemplateFittingSolve2PassTest() {
    std::string json2Pass = jsonStringNoFFT;
    std::string target = "\"singlePass\" : true,";
    size_t pos = json2Pass.find(target);
    json2Pass.replace(pos, target.length(), "\"singlePass\" : false,");

    // TODO do with photometry in a second time
    target = "\"enablePhotometry\" : true,";
    pos = json2Pass.find(target);
    json2Pass.replace(pos, target.length(), "\"enablePhotometry\" : false,");

    spc = fixture_SharedSpectrumExtended().spc;
    Init(jsonString + json2Pass);
  }
};

class fixture_TemplateFittingSolveTestFFT
    : public fixture_TemplateFittingCommon {
public:
  fixture_Context ctx;
  fixture_TemplateFittingSolveTestFFT() {
    spc = fixture_SharedSpectrumExtended().spc;
    Init(jsonString + jsonStringFFT);
  }
};

Int32 EstimateXtYSlow(const TFloat64List &X, const TFloat64List &Y,
                      TFloat64List &XtY) {
  Int32 nShifts = Y.size() - X.size() + 1;
  XtY.resize(nShifts);

  Int32 nX = X.size();
  Float64 xty = 0.0;
  for (std::size_t k = 0; k < nShifts; k++) {
    xty = 0.0;
    for (std::size_t j = 0; j < nX; j++) {
      xty += X[j] * Y[j + k];
    }
    XtY[k] = xty;
  }
  return 0;
}
} // namespace NSEpic

#endif