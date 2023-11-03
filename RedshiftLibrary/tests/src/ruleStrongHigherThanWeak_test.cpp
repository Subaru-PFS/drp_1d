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
#include "RedshiftLibrary/line/ruleStrongHigherThanWeak.h"
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <iostream>

using namespace std;
using namespace boost;

using namespace NSEpic;

class RuleStrongHigherThanWeak_fixture {
public:
  static CLine createLine(std::string name, CLine::EForce force, Int32 id);
  static CLineModelElement createCLineModelElement(CLineVector lines,
                                                   std::vector<Float64> amps);
  static CLineModelElementList
  createElementList(vector<CLineModelElement> elements, Int32 nElements);
  CLineModelElementList makeSomeElementList(Float64 weak1Amp, Float64 weak2Amp,
                                            Float64 strong1Amp,
                                            Float64 strong2Amp);
};

CLine RuleStrongHigherThanWeak_fixture::createLine(std::string name,
                                                   CLine::EForce force,
                                                   Int32 id) {
  return CLine(name, -1., CLine::EType::nType_Emission,
               std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM()), force,
               0., false, undefStr, 1.0, "Em1", id, std::to_string(id));
};

CLineModelElement RuleStrongHigherThanWeak_fixture::createCLineModelElement(
    CLineVector lines, std::vector<Float64> amps) {
  TLineModelElementParam_ptr elementParam_ptr =
      std::make_shared<TLineModelElementParam>(lines, 1.0, 1.1);
  CLineModelElement cLineModelElement(elementParam_ptr, "instrumentDriven");

  vector<Float64> fittedAmplitudesErrorSigmas(lines.size());
  cLineModelElement.m_ElementParam->m_FittedAmplitudes = amps;
  cLineModelElement.m_ElementParam->m_NominalAmplitudes = amps;
  cLineModelElement.m_ElementParam->m_FittedAmplitudeErrorSigmas =
      fittedAmplitudesErrorSigmas;

  // Forces model elements to be included in lambda range
  cLineModelElement.m_OutsideLambdaRangeList.assign(lines.size(), false);

  return cLineModelElement;
}

CLineModelElementList RuleStrongHigherThanWeak_fixture::createElementList(
    vector<CLineModelElement> elements, Int32 nElements) {
  CLineModelElementList elementList;
  int i;
  for (i = 0; i < nElements; i++) {
    elementList.push_back(std::make_shared<CLineModelElement>(elements[i]));
  }
  return elementList;
}

CLineModelElementList RuleStrongHigherThanWeak_fixture::makeSomeElementList(
    Float64 weak1Amp, Float64 weak2Amp, Float64 strong1Amp,
    Float64 strong2Amp) {
  CLineProfile_ptr profilesym =
      std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM());

  CLine lineWeak1 = RuleStrongHigherThanWeak_fixture::createLine(
      "LineWeak1", CLine::EForce::nForce_Weak, 0);
  CLine lineWeak2 = RuleStrongHigherThanWeak_fixture::createLine(
      "LineWeak2", CLine::EForce::nForce_Weak, 1);
  CLine lineStrong1 = RuleStrongHigherThanWeak_fixture::createLine(
      "LineStrong1", CLine::EForce::nForce_Strong, 2);
  CLine lineStrong2 = RuleStrongHigherThanWeak_fixture::createLine(
      "LineStrong2", CLine::EForce::nForce_Strong, 3);
  CLineMap lm;
  CLineModelElement weakElement1 =
      RuleStrongHigherThanWeak_fixture::createCLineModelElement(
          {lineWeak1, lineWeak2}, {weak1Amp, weak2Amp});
  CLineModelElement strongElement1 =
      RuleStrongHigherThanWeak_fixture::createCLineModelElement(
          {lineStrong1, lineStrong2}, {strong1Amp, strong2Amp});

  CLineModelElementList elementList =
      RuleStrongHigherThanWeak_fixture::createElementList(
          {weakElement1, strongElement1}, 2);

  CRuleStrongHigherThanWeak rule;
  CLMEltListVector lmeltlistv = CLMEltListVector(elementList, lm);
  rule.SetUp(true, CLine::EType::nType_Emission);
  rule.Correct(lmeltlistv);
  return elementList;
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

  CLineModelElementList elementList =
      RuleStrongHigherThanWeak_fixture::makeSomeElementList(
          lineWeak1InitialAmplitude, lineWeak2InitialAmplitude,
          lineStrong1InitialAmplitude, lineStrong2InitialAmplitude);

  Float64 correctedAmpWeak1 =
      elementList[0]->m_ElementParam->m_FittedAmplitudes[0];
  Float64 correctedAmpWeak2 =
      elementList[0]->m_ElementParam->m_FittedAmplitudes[1];
  Float64 correctedAmpStrong1 =
      elementList[1]->m_ElementParam->m_FittedAmplitudes[0];
  Float64 correctedAmpStrong2 =
      elementList[1]->m_ElementParam->m_FittedAmplitudes[1];

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

  CLineModelElementList elementList =
      RuleStrongHigherThanWeak_fixture::makeSomeElementList(
          lineWeak1InitialAmplitude, lineWeak2InitialAmplitude,
          lineStrong1InitialAmplitude, lineStrong2InitialAmplitude);

  Float64 correctedAmpWeak1 =
      elementList[0]->m_ElementParam->m_FittedAmplitudes[0];
  Float64 correctedAmpWeak2 =
      elementList[0]->m_ElementParam->m_FittedAmplitudes[1];
  Float64 correctedAmpStrong1 =
      elementList[1]->m_ElementParam->m_FittedAmplitudes[0];
  Float64 correctedAmpStrong2 =
      elementList[1]->m_ElementParam->m_FittedAmplitudes[1];

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
