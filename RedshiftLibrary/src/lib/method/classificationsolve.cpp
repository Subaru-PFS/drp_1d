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
#include "RedshiftLibrary/method/classificationsolve.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/method/classificationresult.h"
#include "RedshiftLibrary/method/solveresult.h"
#include "RedshiftLibrary/processflow/parameterstore.h"

using namespace NSEpic;

CClassificationSolve::CClassificationSolve() : CSolve("classification") {}

std::shared_ptr<CSolveResult> CClassificationSolve::compute() {

  auto const &inputContext = Context.GetInputContext();
  auto const &resultStore = Context.GetResultStore();
  std::map<std::string, std::weak_ptr<const CPdfSolveResult>> results;
  std::map<std::string, Float64> logEvidences;
  std::map<std::string, bool> hasResult;
  std::string stage = "redshiftSolver";
  for (const std::string &spectrumModel : inputContext->m_categories) {
    const std::string &method =
        inputContext->GetParameterStore()->Get<std::string>(
            spectrumModel + "." + stage + ".method");

    if (resultStore->hasSolveResult(spectrumModel, stage, method)) {
      results[spectrumModel] = std::dynamic_pointer_cast<const CPdfSolveResult>(
          resultStore->GetSolveResult(spectrumModel, stage, method).lock());
      hasResult[spectrumModel] = true;
    } else
      hasResult[spectrumModel] = false;

    logEvidences[spectrumModel] = -INFINITY;
  }
  std::shared_ptr<CClassificationResult> classifResult =
      std::make_shared<CClassificationResult>();
  Float64 MaxLogEvidence = -DBL_MAX;
  for (const std::string &spectrumModel : inputContext->m_categories) {
    if (hasResult[spectrumModel]) {
      logEvidences[spectrumModel] =
          results[spectrumModel].lock()->getEvidence();
      if (logEvidences[spectrumModel] > MaxLogEvidence) {
        MaxLogEvidence = logEvidences[spectrumModel];
        typeLabel = spectrumModel;
      }
    }
  }
  Log.LogInfo("Setting object type: %s", typeLabel.c_str());
  Float64 sum = 0.;

  for (const std::string &spectrumModel : inputContext->m_categories) {
    if (hasResult[spectrumModel]) {
      Float64 Proba = exp(logEvidences[spectrumModel] - MaxLogEvidence);
      sum += Proba;
    }
  }
  if (sum <= 0) {
    THROWG(NO_CLASSIFICATION,
           "Classification failed, all probabilities undefined");
  }
  for (const std::string &spectrumModel : inputContext->m_categories) {
    if (hasResult[spectrumModel]) {
      Float64 proba = exp(logEvidences[spectrumModel] - MaxLogEvidence);
      proba /= sum;
      classifResult->SetProba(spectrumModel, proba);
      classifResult->SetEvidence(spectrumModel, logEvidences[spectrumModel]);
    } else {
      classifResult->SetProba(spectrumModel, 0.);
      classifResult->SetEvidence(spectrumModel, -INFINITY);
    }
  }
  classifResult->SetTypeLabel(typeLabel);

  return classifResult;
}
