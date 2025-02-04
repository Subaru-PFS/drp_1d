
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

#include "RedshiftLibrary/linemodel/abstractfitter.h"
#include "RedshiftLibrary/linemodel/element.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/linemodel/individualfitter.h"
#include "RedshiftLibrary/linemodel/obsiterator.h"
#include "RedshiftLibrary/method/linemodelsolve.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/context.h"
#include <boost/test/unit_test.hpp>

// namespace NSEpic

using namespace NSEpic;

class CTestFitter : public CAbstractFitter {
public:
  using CAbstractFitter::CAbstractFitter;

private:
  friend struct AbstractFitterFixture;
  void doFit(Float64 redshift) override{};
};

struct AbstractFitterFixture {

  // Creates necessary objects
  CTLambdaRangePtrVector lambdaRanges;
  std::shared_ptr<Int32> curObsPtr = std::make_shared<Int32>(0);
  Int32 nb_spectra = 1;
  std::unique_ptr<CSpectraGlobalIndex> spcIndex =
      std::make_unique<CSpectraGlobalIndex>(nb_spectra);
  CLineMap restLineList = CLineMap();

  std::shared_ptr<CLMEltListVector> elementsVector =
      std::make_shared<CLMEltListVector>(
          lambdaRanges, *spcIndex, restLineList,
          ElementComposition::EmissionAbsorption);
  CCSpectrumVectorPtr inputSpcs;
  CSpcModelVectorPtr spectrumModels;

  CTestFitter initializeTestFitter(std::string method, std::string igmFit) {
    CAutoScope autoscope1(Context.m_ScopeStack, "l1");
    CAutoScope autoscope2(Context.m_ScopeStack, "l2", ScopeType::STAGE);
    CAutoScope autoscope3(Context.m_ScopeStack, method, ScopeType::METHOD);
    const std::string jsonString = createInputJson(method, igmFit);
    Context.LoadParameterStore(jsonString);
    CTestFitter testFitter(elementsVector, inputSpcs, lambdaRanges,
                           spectrumModels, restLineList, *spcIndex);
    return testFitter;
  }

  std::string createInputJson(std::string method, std::string profile) {
    std::string methodTargetWord = "someMethod";
    std::string profileTargetWord = "someProfile";

    std::string jsonString = {
        "{\"l1\" : {\"l2\": {\"someMethod\": {\"lineModel\": {"
        "\"lya\": {"
        "\"profile\": \"someProfile\","
        "\"asymProfile\": {"
        "\"asymFitMin\": 1.0,"
        "\"asymFitMax\": 2.0,"
        "\"asymFitStep\": 3.0,"
        "\"widthFitMin\": 4.0,"
        "\"widthFitMax\": 5.0,"
        "\"widthFitStep\": 6.0,"
        "\"deltaFitMin\": 7.0,"
        "\"deltaFitMax\": 8.0,"
        "\"deltaStepMax\": 9.0"
        "}}}}}}}"};

    size_t methodPosition = jsonString.find(methodTargetWord);
    jsonString.replace(methodPosition, methodTargetWord.length(), method);

    size_t profilePosition = jsonString.find(profileTargetWord);
    jsonString.replace(profilePosition, profileTargetWord.length(), profile);

    return jsonString;
  }

  void testLyaValues(CTestFitter testFitter, bool nanExpected = false) {
    if (nanExpected) {
      BOOST_CHECK(isnan(testFitter.m_opt_lya_fit_asym_min));
      BOOST_CHECK(isnan(testFitter.m_opt_lya_fit_asym_max));
      BOOST_CHECK(isnan(testFitter.m_opt_lya_fit_asym_step));
      BOOST_CHECK(isnan(testFitter.m_opt_lya_fit_width_min));
      BOOST_CHECK(isnan(testFitter.m_opt_lya_fit_width_max));
      BOOST_CHECK(isnan(testFitter.m_opt_lya_fit_width_step));
      BOOST_CHECK(isnan(testFitter.m_opt_lya_fit_delta_min));
      BOOST_CHECK(isnan(testFitter.m_opt_lya_fit_delta_max));
      BOOST_CHECK(isnan(testFitter.m_opt_lya_fit_delta_step));
    } else {
      BOOST_CHECK(testFitter.m_opt_lya_fit_asym_min == 1.0);
      BOOST_CHECK(testFitter.m_opt_lya_fit_asym_max == 2.0);
      BOOST_CHECK(testFitter.m_opt_lya_fit_asym_step == 3.0);
      BOOST_CHECK(testFitter.m_opt_lya_fit_width_min == 4.0);
      BOOST_CHECK(testFitter.m_opt_lya_fit_width_max == 5.0);
      BOOST_CHECK(testFitter.m_opt_lya_fit_width_step == 6.0);
      BOOST_CHECK(testFitter.m_opt_lya_fit_delta_min == 7.0);
      BOOST_CHECK(testFitter.m_opt_lya_fit_delta_max == 8.0);
      BOOST_CHECK(testFitter.m_opt_lya_fit_delta_step == 9.0);
    }
  }
};

BOOST_FIXTURE_TEST_SUITE(abstractFitter_test, AbstractFitterFixture)

BOOST_AUTO_TEST_CASE(lineModelSolveWithIgmFit_test) {
  // Checks that values are correctly set if method is lineModelSolve and
  // profile is asym
  std::string method = "lineModelSolve";
  CTestFitter testFitter = initializeTestFitter(method, "asym");
  testLyaValues(testFitter);
}

BOOST_AUTO_TEST_CASE(lineModelSolveWithoutIgmFit_test) {
  // Checks that values are set to NAN if profile is igm (even if method is
  // lineModelSolve)
  std::string method = "lineModelSolve";
  CTestFitter testFitter = initializeTestFitter(method, "igm");
  testLyaValues(testFitter, true);
}

BOOST_AUTO_TEST_CASE(lineMeasSolveWithIgmFit_test) {
  // Checks that values are correctly set if method is lineMeasSolve and profile
  // is asym
  std::string method = "lineMeasSolve";
  CTestFitter testFitter = initializeTestFitter(method, "asym");
  testLyaValues(testFitter);
}

BOOST_AUTO_TEST_CASE(lineMeasSolveWithoutIgmFit_test) {
  // Checks that values are set to NAN if profile is igm (even if method is
  // lineMeasSolve)
  std::string method = "lineMeasSolve";
  CTestFitter testFitter = initializeTestFitter(method, "igm");
  testLyaValues(testFitter, true);
}
BOOST_AUTO_TEST_CASE(otherMethod_test) {
  // Checks that values are correctly set to NAN if method is other than defined
  // ones even if profile is asym
  std::string method = "other";
  CTestFitter testFitter = initializeTestFitter(method, "asym");
  testLyaValues(testFitter, true);
}

BOOST_AUTO_TEST_SUITE_END()
