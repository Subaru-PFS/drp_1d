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
#include <iostream>
#include <memory>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/line/regulament.h"
#include "RedshiftLibrary/log/log.h"

// To be removed once JSON code is in <--
#include "RedshiftLibrary/line/rule2singlelinesamplitude.h"
#include "RedshiftLibrary/line/ruleOIIRatioRange.h"
#include "RedshiftLibrary/line/ruleStrongHigherThanWeak.h"
#include "RedshiftLibrary/line/ruleSuperStrongHighest.h"
// -->

using namespace NSEpic;
using namespace std;

void CRegulament::Apply(CLineModelElementList &LineModelElementList) {
  for (unique_ptr<CRule> &rule : m_RulesVector) {
    rule->Apply(LineModelElementList);
    std::string logRule = rule->GetLogs();
    if (logRule.size() > 0 && m_LogsEnabled) {
      m_RulesLog.push_back(logRule);
    }
  }
}
/*
void CRegulament::ApplyWithRedshift( Float64 Redshift )
{

}
*/
bool CRegulament::CreateRulesFromJSONFiles(void) {
  // To be removed once JSON code is in <--
  bool True = true;

  // m_RulesVector.push_back( std::make_unique<CRuleBalmerLinearSolver>() );
  // m_RulesVector.back()->SetUp( True );

  m_RulesVector.push_back(
      unique_ptr<CRule2SingleLinesAmplitude>(new CRule2SingleLinesAmplitude()));
  m_RulesVector.back()->SetUp(True, CLine::nType_Emission, linetags::halpha_em,
                              linetags::hbeta_em, 1.0 / 2.86 * 1.1);

  m_RulesVector.push_back(
      unique_ptr<CRule2SingleLinesAmplitude>(new CRule2SingleLinesAmplitude()));
  m_RulesVector.back()->SetUp(True, CLine::nType_Emission, linetags::hbeta_em,
                              linetags::hgamma_em, 0.47 * 1.1);

  m_RulesVector.push_back(
      unique_ptr<CRule2SingleLinesAmplitude>(new CRule2SingleLinesAmplitude()));
  m_RulesVector.back()->SetUp(True, CLine::nType_Emission, linetags::hgamma_em,
                              linetags::hdelta_em, 1.1);

  m_RulesVector.push_back(
      unique_ptr<CRule2SingleLinesAmplitude>(new CRule2SingleLinesAmplitude()));
  m_RulesVector.back()->SetUp(True, CLine::nType_Emission, linetags::hdelta_em,
                              linetags::h8_em, 1.1);

  m_RulesVector.push_back(
      unique_ptr<CRule2SingleLinesAmplitude>(new CRule2SingleLinesAmplitude()));
  m_RulesVector.back()->SetUp(True, CLine::nType_Emission, linetags::h8_em,
                              linetags::h9_em, 1.1);

  m_RulesVector.push_back(
      unique_ptr<CRule2SingleLinesAmplitude>(new CRule2SingleLinesAmplitude()));
  m_RulesVector.back()->SetUp(True, CLine::nType_Emission, linetags::h9_em,
                              linetags::h10_em, 1.1);

  m_RulesVector.push_back(
      unique_ptr<CRule2SingleLinesAmplitude>(new CRule2SingleLinesAmplitude()));
  m_RulesVector.back()->SetUp(True, CLine::nType_Emission, linetags::h10_em,
                              linetags::h11_em, 1.1);

  m_RulesVector.push_back(
      unique_ptr<CRule2SingleLinesAmplitude>(new CRule2SingleLinesAmplitude()));
  m_RulesVector.back()->SetUp(True, CLine::nType_Emission, linetags::hbeta_abs,
                              linetags::hgamma_abs, 1.1);

  m_RulesVector.push_back(
      unique_ptr<CRule2SingleLinesAmplitude>(new CRule2SingleLinesAmplitude()));
  m_RulesVector.back()->SetUp(True, CLine::nType_Emission, linetags::hgamma_abs,
                              linetags::hdelta_abs, 1.1);

  m_RulesVector.push_back(
      unique_ptr<CRule2SingleLinesAmplitude>(new CRule2SingleLinesAmplitude()));
  m_RulesVector.back()->SetUp(True, CLine::nType_Emission, linetags::hdelta_abs,
                              linetags::h8_abs, 1.1);

  m_RulesVector.push_back(
      unique_ptr<CRule2SingleLinesAmplitude>(new CRule2SingleLinesAmplitude()));
  m_RulesVector.back()->SetUp(True, CLine::nType_Emission, linetags::h8_abs,
                              linetags::h9_abs, 1.1);

  m_RulesVector.push_back(
      unique_ptr<CRule2SingleLinesAmplitude>(new CRule2SingleLinesAmplitude()));
  m_RulesVector.back()->SetUp(True, CLine::nType_Emission, linetags::h9_abs,
                              linetags::h10_abs, 1.1);

  m_RulesVector.push_back(
      unique_ptr<CRule2SingleLinesAmplitude>(new CRule2SingleLinesAmplitude()));
  m_RulesVector.back()->SetUp(True, CLine::nType_Emission, linetags::h10_abs,
                              linetags::h11_abs, 1.1);

  m_RulesVector.push_back(unique_ptr<CRuleRatioRange>(new CRuleRatioRange()));
  m_RulesVector.back()->SetUp(True, CLine::nType_Emission, linetags::oII3726_em,
                              linetags::oII3729_em, 2.5);

  m_RulesVector.push_back(
      unique_ptr<CRuleStrongHigherThanWeak>(new CRuleStrongHigherThanWeak()));
  m_RulesVector.back()->SetUp(True, CLine::nType_Emission);

  m_RulesVector.push_back(
      unique_ptr<CRuleStrongHigherThanWeak>(new CRuleStrongHigherThanWeak()));
  m_RulesVector.back()->SetUp(True, CLine::nType_Absorption);

  m_RulesVector.push_back(unique_ptr<CRuleRatioRange>(new CRuleRatioRange()));
  m_RulesVector.back()->SetUp(True, CLine::nType_Emission,
                              linetags::cIII1907_em, linetags::cIII1909_em,
                              2.0);

  // OII and Halpha Super Strong
  // m_RulesVector.push_back( unique_ptr<CRuleSuperStrong>(new
  // CRuleSuperStrong()) ); m_RulesVector.back()->SetUp( True,
  // CLine::nType_Emission, linetags::oII3726_em, linetags::oII3729_em,
  // linetags::halpha_em, 1.1 );

  if (m_LogsEnabled) {
    m_RulesLog.push_back("Linemodel-Regulament: rules creation");
  }
  return true;
}

void CRegulament::EnableRulesAccordingToParameters(std::string Parameters) {
  if (Parameters == "no") {
    return;
  }
  for (unique_ptr<CRule> &rule : m_RulesVector) {
    bool enableRule = Parameters.find(rule->Name) != std::string::npos;
    if (Parameters == "all" || enableRule) {
      Log.LogDebug("Enabling rule %s.", rule->Name.c_str());
      rule->Enabled = true;
    } else {
      Log.LogDebug("Disabling rule %s.", rule->Name.c_str());
      rule->Enabled = false;
    }
  }
}

void CRegulament::EnableLogs(bool enable) {
  m_LogsEnabled = enable;
  m_RulesLog.clear();
}

const std::vector<string> &CRegulament::GetLogs() const { return m_RulesLog; }
