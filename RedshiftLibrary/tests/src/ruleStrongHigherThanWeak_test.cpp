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
#include <cmath>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include "RedshiftLibrary/line/ruleStrongHigherThanWeak.h"
#include "RedshiftLibrary/linemodel/obsiterator.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace std;
using namespace boost;

using namespace NSEpic;

class RuleStrongHigherThanWeak_fixture {
public:
  void setLineModelElementsAmplitudes(TFloat64List amps);

  CLineCatalog line_catalog = makeLineCatalog();
  CLMEltListVector element_list_vector = makeElementListVector();

private:
  CLine createLine(std::string const &name, CLine::EForce force, Int32 id,
                   std::string const &group);
  CLineCatalog makeLineCatalog();
  CLMEltListVector makeElementListVector();
};

CLine RuleStrongHigherThanWeak_fixture::createLine(std::string const &name,
                                                   CLine::EForce force,
                                                   Int32 id,
                                                   std::string const &group) {
  return CLine(name, -1., CLine::EType::nType_Emission,
               std::make_unique<CLineProfileSYM>(), force, 0., false, group,
               1.0, "Em1", id, std::to_string(id));
};

CLineCatalog RuleStrongHigherThanWeak_fixture::makeLineCatalog() {
  CLineCatalog line_catalog;
  CLineProfile_ptr profilesym = std::make_unique<CLineProfileSYM>();
  line_catalog.Add(
      createLine("LineWeak1", CLine::EForce::nForce_Weak, 0, "groupWeak"));
  line_catalog.Add(
      createLine("LineWeak2", CLine::EForce::nForce_Weak, 1, "groupWeak"));
  line_catalog.Add(createLine("LineStrong1", CLine::EForce::nForce_Strong, 2,
                              "groupStrong"));
  line_catalog.Add(createLine("LineStrong2", CLine::EForce::nForce_Strong, 3,
                              "groupStrong"));
  return line_catalog;
}

CLMEltListVector RuleStrongHigherThanWeak_fixture::makeElementListVector() {
  CLineMap const &lm = line_catalog.GetList();
  CSpectraGlobalIndex spc_index = CSpectraGlobalIndex(1);
  CAutoScope autoscope1(Context.m_ScopeStack, "model",
                        ScopeType::SPECTRUMMODEL);
  CAutoScope autoscope2(Context.m_ScopeStack, "stage", ScopeType::STAGE);
  CAutoScope autoscope3(Context.m_ScopeStack, "methodSolve", ScopeType::METHOD);
  std::string const jsonString = {
      "{\"model\" : {\"stage\": {\"methodSolve\": {\"lineModel\": {"
      "\"velocityEmission\": 100,"
      "\"lineWidthType\": \"instrumentDriven\""
      "}}}}}"};
  Context.LoadParameterStore(jsonString);
  CLMEltListVector element_list_vector =
      CLMEltListVector(spc_index, lm, ElementComposition::Default);
  spc_index.setAtBegining();
  for (auto &elt : element_list_vector.getElementList()) {
    elt->m_OutsideLambdaRangeList.assign(2, false);
    elt->computeOutsideLambdaRange();
  }
  return element_list_vector;
}

void RuleStrongHigherThanWeak_fixture::setLineModelElementsAmplitudes(
    TFloat64List amps) {
  auto &param1 = *(element_list_vector.getElementParam()[0]);
  param1.m_FittedAmplitudes = {amps[2], amps[3]};
  param1.m_NominalAmplitudes = {amps[2], amps[3]};
  param1.m_FittedAmplitudesStd = TFloat64List(param1.size());

  auto &param2 = *(element_list_vector.getElementParam()[1]);
  param2.m_FittedAmplitudes = {amps[0], amps[1]};
  param2.m_NominalAmplitudes = {amps[0], amps[1]};
  param2.m_FittedAmplitudesStd = TFloat64List(param2.size());
}

BOOST_FIXTURE_TEST_SUITE(RuleStrongHigherThanWeak_test,
                         RuleStrongHigherThanWeak_fixture)

BOOST_AUTO_TEST_CASE(Correct_test_no_change) {
  // Checks that if strongest weak line initial amplitude is lower than weakest
  // strong line, amplitudes are not changed
  Float64 lineWeak1InitialAmplitude = 0.8;
  Float64 lineWeak2InitialAmplitude = 1.;
  Float64 lineStrong1InitialAmplitude = 1.8;
  Float64 lineStrong2InitialAmplitude = 2.;

  setLineModelElementsAmplitudes(
      {lineWeak1InitialAmplitude, lineWeak2InitialAmplitude,
       lineStrong1InitialAmplitude, lineStrong2InitialAmplitude});

  CRuleStrongHigherThanWeak rule;
  rule.SetUp(true, CLine::EType::nType_Emission);
  rule.Correct(element_list_vector);

  auto &params = element_list_vector.getElementParam();
  Float64 correctedAmpWeak1 = params[1]->m_FittedAmplitudes[0];
  Float64 correctedAmpWeak2 = params[1]->m_FittedAmplitudes[1];
  Float64 correctedAmpStrong1 = params[0]->m_FittedAmplitudes[0];
  Float64 correctedAmpStrong2 = params[0]->m_FittedAmplitudes[1];

  BOOST_CHECK(correctedAmpWeak1 == lineWeak1InitialAmplitude);
  BOOST_CHECK(correctedAmpWeak2 == lineWeak2InitialAmplitude);
  BOOST_CHECK(correctedAmpStrong1 == lineStrong1InitialAmplitude);
  BOOST_CHECK(correctedAmpStrong2 == lineStrong2InitialAmplitude);
}

BOOST_AUTO_TEST_CASE(Correct_test_one_high_weak) {
  // Checks that if strongest weak line initial amplitude is higher than weakest
  // strong line, amplitudes of both weak line element are reduced
  Float64 lineWeak1InitialAmplitude = 0.8;
  Float64 lineWeak2InitialAmplitude = 4.;

  Float64 lineStrong1InitialAmplitude = 2;
  Float64 lineStrong2InitialAmplitude = 2.2;
  setLineModelElementsAmplitudes(
      {lineWeak1InitialAmplitude, lineWeak2InitialAmplitude,
       lineStrong1InitialAmplitude, lineStrong2InitialAmplitude});

  CRuleStrongHigherThanWeak rule;

  rule.SetUp(true, CLine::EType::nType_Emission);
  rule.Correct(element_list_vector);

  auto &params = element_list_vector.getElementParam();
  Float64 correctedAmpWeak1 = params[1]->m_FittedAmplitudes[0];
  Float64 correctedAmpWeak2 = params[1]->m_FittedAmplitudes[1];
  Float64 correctedAmpStrong1 = params[0]->m_FittedAmplitudes[0];
  Float64 correctedAmpStrong2 = params[0]->m_FittedAmplitudes[1];

  Float64 reductionCoef =
      lineStrong1InitialAmplitude / lineWeak2InitialAmplitude;
  cout << "correctedAmpWeak1 " << correctedAmpWeak1
       << "lineWeak1InitialAmplitude * reductionCoef "
       << lineWeak1InitialAmplitude * reductionCoef;
  BOOST_CHECK(correctedAmpWeak1 == lineWeak1InitialAmplitude * reductionCoef);
  BOOST_CHECK(correctedAmpWeak2 == lineWeak2InitialAmplitude * reductionCoef);
  BOOST_CHECK(correctedAmpStrong1 == lineStrong1InitialAmplitude);
  BOOST_CHECK(correctedAmpStrong2 == lineStrong2InitialAmplitude);
}

BOOST_AUTO_TEST_SUITE_END()
