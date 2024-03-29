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

CClassificationSolve::CClassificationSolve(TScopeStack &scope,
                                           std::string objectType)
    : CSolve("classification", scope, objectType) {}

std::shared_ptr<CSolveResult>
CClassificationSolve::compute(std::shared_ptr<const CInputContext> inputContext,
                              std::shared_ptr<COperatorResultStore> resultStore,
                              TScopeStack &scope) {

  std::map<std::string, std::weak_ptr<const CPdfSolveResult>> results;
  std::map<std::string, Float64> logEvidences;
  std::map<std::string, bool> hasResult;
  for (const std::string &category : inputContext->m_categories) {
    const std::string &method =
        inputContext->GetParameterStore()->Get<std::string>(category +
                                                            ".method");

    if (resultStore->hasSolveResult(category, method)) {
      results[category] = std::dynamic_pointer_cast<const CPdfSolveResult>(
          resultStore->GetSolveResult(category, method).lock());
      hasResult[category] = true;
    } else
      hasResult[category] = false;

    logEvidences[category] = -INFINITY;
  }
  std::shared_ptr<CClassificationResult> classifResult =
      std::make_shared<CClassificationResult>();
  Float64 MaxLogEvidence = -DBL_MAX;
  for (const std::string &category : inputContext->m_categories) {
    if (hasResult[category]) {
      logEvidences[category] = results[category].lock()->getEvidence();
      if (logEvidences[category] > MaxLogEvidence) {
        MaxLogEvidence = logEvidences[category];
        typeLabel = category;
      }
    }
  }
  Log.LogInfo("Setting object type: %s", typeLabel.c_str());
  Float64 sum = 0.;

  for (const std::string &category : inputContext->m_categories) {
    if (hasResult[category]) {
      Float64 Proba = exp(logEvidences[category] - MaxLogEvidence);
      sum += Proba;
    }
  }
  if (sum <= 0) {
    THROWG(NO_CLASSIFICATION,
           "Classification failed, all probabilities undefined");
  }
  for (const std::string &category : inputContext->m_categories) {
    if (hasResult[category]) {
      Float64 proba = exp(logEvidences[category] - MaxLogEvidence);
      proba /= sum;
      classifResult->SetProba(category, proba);
      classifResult->SetEvidence(category, logEvidences[category]);
    } else {
      classifResult->SetProba(category, 0.);
      classifResult->SetEvidence(category, -INFINITY);
    }
  }
  classifResult->SetTypeLabel(typeLabel);

  return classifResult;
}
