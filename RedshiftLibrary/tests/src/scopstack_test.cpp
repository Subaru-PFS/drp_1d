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

#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/processflow/scopestack.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(scopestack)

CScopeStack scope;

BOOST_AUTO_TEST_CASE(scopestack_ctor_test) {

  CScopeStack cscope0;
  BOOST_CHECK(cscope0.size() == 0);

  TScopeStack scope{"level_1"};
  CScopeStack cscope1(scope);
  BOOST_CHECK(cscope1.size() == 1);

  CScopeStack cscope2(TScopeStack{"level_1"});
  BOOST_CHECK(cscope2.size() == 1);

  CScopeStack cscope3{"level_1"};
  BOOST_CHECK(cscope3.size() == 1);
}

BOOST_AUTO_TEST_CASE(scopestack_info_test) {
  CScopeStack cscope0;
  BOOST_CHECK(cscope0.empty());
  BOOST_CHECK(cscope0.size() == 0);

  CScopeStack cscope1{"level_1"};
  BOOST_CHECK(!cscope1.empty());
  BOOST_CHECK(cscope1.size() == 1);

  CScopeStack cscope2{"level_1", "level_2"};
  BOOST_CHECK(!cscope2.empty());
  BOOST_CHECK(cscope2.size() == 2);
}

BOOST_AUTO_TEST_CASE(scopestack_accessors_test) {
  CScopeStack cscope0;
  BOOST_CHECK_THROW(cscope0.at(0), std::exception);

  CScopeStack const &cscope0_const = cscope0;
  BOOST_CHECK_THROW(cscope0_const.at(0), std::exception);

  CScopeStack cscope1{"level_1"};
  BOOST_CHECK(cscope1.front() == "level_1");
  BOOST_CHECK(cscope1.back() == "level_1");
  BOOST_CHECK(cscope1.at(0) == "level_1");
  BOOST_CHECK(cscope1[0] == "level_1");

  CScopeStack cscope2{"level_1", "level_2"};
  BOOST_CHECK(cscope2.front() == "level_1");
  BOOST_CHECK(cscope2.back() == "level_2");
  BOOST_CHECK(cscope2.at(0) == "level_1");
  BOOST_CHECK(cscope2.at(1) == "level_2");
  BOOST_CHECK(cscope2[0] == "level_1");
  BOOST_CHECK(cscope2[1] == "level_2");
}

BOOST_AUTO_TEST_CASE(scopestack_iterators_test) {
  CScopeStack cscope0;
  BOOST_CHECK(cscope0.begin() == cscope0.end());

  CScopeStack cscope1{"level_1"};
  BOOST_CHECK(*cscope1.begin() == "level_1");
  BOOST_CHECK(*(cscope1.end() - 1) == "level_1");
  size_t i = 0;
  for (auto scope : cscope1) {
    BOOST_CHECK(scope == cscope1[i++]);
  }

  CScopeStack cscope2{"level_1", "level_2"};
  BOOST_CHECK(*cscope2.begin() == "level_1");
  BOOST_CHECK(*(cscope2.end() - 1) == "level_2");
  i = 0;
  for (auto scope : cscope2) {
    BOOST_CHECK(scope == cscope2[i++]);
  }
}

BOOST_AUTO_TEST_CASE(scopestack_push_test) {
  CScopeStack cscope0;
  BOOST_CHECK_NO_THROW(cscope0.push_back("level_1"));
  BOOST_CHECK(cscope0.size() == 1);
  BOOST_CHECK(cscope0[0] == "level_1");

  std::string const level2 = "level_2";
  BOOST_CHECK_NO_THROW(cscope0.push_back(level2));
  BOOST_CHECK(cscope0.size() == 2);
  BOOST_CHECK(cscope0.back() == level2);

  BOOST_CHECK_NO_THROW(cscope0.pop_back());
  BOOST_CHECK(cscope0.size() == 1);
  BOOST_CHECK(cscope0.back() == "level_1");

  CScopeStack cscope1;
  BOOST_CHECK(!cscope1.has_type(ScopeType::UNDEFINED));
  BOOST_CHECK(!cscope1.has_type(ScopeType::STAGE));
  BOOST_CHECK(!cscope1.has_type(ScopeType::SPECTRUMMODEL));
  BOOST_CHECK(!cscope1.has_type(ScopeType::METHOD));
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::UNDEFINED),
                    std::exception);
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::STAGE), std::exception);
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::SPECTRUMMODEL),
                    std::exception);
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::METHOD), std::exception);
  BOOST_CHECK(cscope1.get_current_type() == ScopeType::UNDEFINED);

  BOOST_CHECK_NO_THROW(cscope1.push_back("stage1", ScopeType::STAGE));
  BOOST_CHECK(!cscope1.has_type(ScopeType::UNDEFINED));
  BOOST_CHECK(cscope1.has_type(ScopeType::STAGE));
  BOOST_CHECK(!cscope1.has_type(ScopeType::SPECTRUMMODEL));
  BOOST_CHECK(!cscope1.has_type(ScopeType::METHOD));
  BOOST_CHECK(cscope1.get_type_value(ScopeType::STAGE) == "stage1");
  BOOST_CHECK(cscope1.get_type_level(ScopeType::STAGE) == 0);
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::UNDEFINED),
                    std::exception);
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::SPECTRUMMODEL),
                    std::exception);
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::METHOD), std::exception);
  BOOST_CHECK(cscope1.get_current_type() == ScopeType::STAGE);

  BOOST_CHECK_NO_THROW(cscope1.pop_back()); // stage1
  BOOST_CHECK(!cscope1.has_type(ScopeType::UNDEFINED));
  BOOST_CHECK(!cscope1.has_type(ScopeType::STAGE));
  BOOST_CHECK(!cscope1.has_type(ScopeType::SPECTRUMMODEL));
  BOOST_CHECK(!cscope1.has_type(ScopeType::METHOD));
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::UNDEFINED),
                    std::exception);
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::STAGE), std::exception);
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::SPECTRUMMODEL),
                    std::exception);
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::METHOD), std::exception);
  BOOST_CHECK(cscope1.get_current_type() == ScopeType::UNDEFINED);

  BOOST_CHECK_NO_THROW(cscope1.push_back("stage1", ScopeType::STAGE));
  // cannot push again STAGE after STAGE
  BOOST_CHECK_THROW(cscope1.push_back("stage2", ScopeType::STAGE),
                    std::exception);
  // cannot push SPECTRUMMODEL after STAGE
  BOOST_CHECK_THROW(cscope1.push_back("model1", ScopeType::SPECTRUMMODEL),
                    std::exception);
  BOOST_CHECK_NO_THROW(cscope1.push_back("garbage"));
  // cannot push METHOD not following STAGE
  BOOST_CHECK_THROW(cscope1.push_back("method1", ScopeType::METHOD),
                    std::exception);
  BOOST_CHECK_NO_THROW(cscope1.pop_back()); // garbage1

  // push METHOD following STAGE
  BOOST_CHECK_NO_THROW(cscope1.push_back("method1", ScopeType::METHOD));
  BOOST_CHECK(!cscope1.has_type(ScopeType::UNDEFINED));
  BOOST_CHECK(cscope1.has_type(ScopeType::STAGE));
  BOOST_CHECK(!cscope1.has_type(ScopeType::SPECTRUMMODEL));
  BOOST_CHECK(cscope1.has_type(ScopeType::METHOD));
  BOOST_CHECK(cscope1.get_current_type() == ScopeType::METHOD);
  BOOST_CHECK(cscope1.get_type_value(ScopeType::STAGE) == "stage1");
  BOOST_CHECK(cscope1.get_type_level(ScopeType::STAGE) == 0);
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::UNDEFINED),
                    std::exception);
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::SPECTRUMMODEL),
                    std::exception);
  BOOST_CHECK(cscope1.get_type_value(ScopeType::METHOD) == "method1");
  BOOST_CHECK(cscope1.get_type_level(ScopeType::METHOD) == 1);
  // cannot push again METHOD after METHOD
  BOOST_CHECK_THROW(cscope1.push_back("method2", ScopeType::METHOD),
                    std::exception);

  BOOST_CHECK_NO_THROW(cscope1.pop_back()); // method1
  BOOST_CHECK(!cscope1.has_type(ScopeType::UNDEFINED));
  BOOST_CHECK(cscope1.has_type(ScopeType::STAGE));
  BOOST_CHECK(!cscope1.has_type(ScopeType::SPECTRUMMODEL));
  BOOST_CHECK(!cscope1.has_type(ScopeType::METHOD));
  BOOST_CHECK(cscope1.get_type_value(ScopeType::STAGE) == "stage1");
  BOOST_CHECK(cscope1.get_type_level(ScopeType::STAGE) == 0);
  BOOST_CHECK(cscope1.get_current_type() == ScopeType::STAGE);
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::UNDEFINED),
                    std::exception);
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::SPECTRUMMODEL),
                    std::exception);
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::METHOD), std::exception);

  BOOST_CHECK_NO_THROW(cscope1.pop_back()); // stage1
  BOOST_CHECK(!cscope1.has_type(ScopeType::UNDEFINED));
  BOOST_CHECK(!cscope1.has_type(ScopeType::STAGE));
  BOOST_CHECK(!cscope1.has_type(ScopeType::SPECTRUMMODEL));
  BOOST_CHECK(!cscope1.has_type(ScopeType::METHOD));
  BOOST_CHECK(cscope1.get_current_type() == ScopeType::UNDEFINED);
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::UNDEFINED),
                    std::exception);
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::STAGE), std::exception);
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::SPECTRUMMODEL),
                    std::exception);
  BOOST_CHECK_THROW(cscope1.get_type_value(ScopeType::METHOD), std::exception);

  CScopeStack cscope2;
  BOOST_CHECK_NO_THROW(cscope2.push_back("model1", ScopeType::SPECTRUMMODEL));
  BOOST_CHECK(!cscope2.has_type(ScopeType::UNDEFINED));
  BOOST_CHECK(!cscope2.has_type(ScopeType::STAGE));
  BOOST_CHECK(cscope2.has_type(ScopeType::SPECTRUMMODEL));
  BOOST_CHECK(!cscope2.has_type(ScopeType::METHOD));
  BOOST_CHECK(cscope2.get_current_type() == ScopeType::SPECTRUMMODEL);
  BOOST_CHECK_THROW(cscope2.get_type_value(ScopeType::STAGE), std::exception);
  BOOST_CHECK_THROW(cscope2.get_type_value(ScopeType::UNDEFINED),
                    std::exception);
  BOOST_CHECK(cscope2.get_type_value(ScopeType::SPECTRUMMODEL) == "model1");
  BOOST_CHECK(cscope2.get_type_level(ScopeType::SPECTRUMMODEL) == 0);
  BOOST_CHECK_THROW(cscope2.get_type_value(ScopeType::METHOD), std::exception);
  // cannot push again MODEL after MODEL
  BOOST_CHECK_THROW(cscope2.push_back("model2", ScopeType::SPECTRUMMODEL),
                    std::exception);
  // cannot push METHOD not following STAGE
  BOOST_CHECK_THROW(cscope2.push_back("method1", ScopeType::METHOD),
                    std::exception);

  BOOST_CHECK_NO_THROW(cscope2.push_back("stage1", ScopeType::STAGE));
  BOOST_CHECK(!cscope2.has_type(ScopeType::UNDEFINED));
  BOOST_CHECK(cscope2.has_type(ScopeType::STAGE));
  BOOST_CHECK(cscope2.has_type(ScopeType::SPECTRUMMODEL));
  BOOST_CHECK(!cscope2.has_type(ScopeType::METHOD));
  BOOST_CHECK(cscope2.get_current_type() == ScopeType::STAGE);
  BOOST_CHECK(cscope2.get_type_value(ScopeType::STAGE) == "stage1");
  BOOST_CHECK(cscope2.get_type_level(ScopeType::STAGE) == 1);
  BOOST_CHECK_THROW(cscope2.get_type_value(ScopeType::UNDEFINED),
                    std::exception);
  BOOST_CHECK(cscope2.get_type_value(ScopeType::SPECTRUMMODEL) == "model1");
  BOOST_CHECK(cscope2.get_type_level(ScopeType::SPECTRUMMODEL) == 0);
  BOOST_CHECK_THROW(cscope2.get_type_value(ScopeType::METHOD), std::exception);
}

BOOST_AUTO_TEST_SUITE_END()