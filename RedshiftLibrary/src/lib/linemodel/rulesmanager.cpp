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

#include "RedshiftLibrary/linemodel/rulesmanager.h"
#include "RedshiftLibrary/linemodel/spectrummodel.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;

CRulesManager::CRulesManager(
    const CLMEltListVectorPtr &elementsVector, const CSpcModelVectorPtr &models,
    const CCSpectrumVectorPtr &inputSpcs,
    const CTLambdaRangePtrVector &lambdaRanges,
    std::shared_ptr<CContinuumManager> continuumManager,
    const TLineVector &restLineList)
    : CLineRatioManager(elementsVector, models, inputSpcs, lambdaRanges,
                        continuumManager, restLineList) {}

Float64 CRulesManager::computeMerit(Int32 iratio) {
  applyRules(true);
  getModel().refreshModel();
  return getLeastSquareMerit();
}

/**
 * /brief Calls the rules' methods depending on the JSON options.
 * If m_rulesoption is "no", do nothing.
 * If either "balmer" or "all" is in the rules string, call
 *ApplyBalmerRuleLinSolve. If "all" or "ratiorange" is in the rules string,
 *call ApplyAmplitudeRatioRangeRule parameterized for OII. If "all" or
 *"strongweak" is in the rules string, call ApplyStrongHigherWeakRule for
 *emission and then for absorption.
 **/
void CRulesManager::applyRules(bool enableLogs) {
  if (m_rulesoption == "no") {
    return;
  }

  m_Regulament.EnableLogs(enableLogs);
  m_Regulament.Apply(getElementList());
}

const TStringList &CRulesManager::GetModelRulesLog() const {
  return m_Regulament.GetLogs();
}

void CRulesManager::setRulesOption(std::string rulesOption) {

  m_Regulament.CreateRulesFromJSONFiles();
  if (rulesOption == "") {
    CAutoScope autoscope(Context.m_ScopeStack, "linemodel");
    std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();
    m_rulesoption = ps->GetScoped<std::string>("rules");
  } else
    m_rulesoption = rulesOption;

  m_Regulament.EnableRulesAccordingToParameters(m_rulesoption);
}
