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
#include "RedshiftLibrary/processflow/scopestore.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(ScopeStore)

void testConstructor(CScopeStore store_ref, CScopeStore store_test) {
  BOOST_CHECK(store_ref.GetCurrentScopeName() ==
              store_test.GetCurrentScopeName());
  BOOST_CHECK(store_ref.GetScopedName("last_level") ==
              store_test.GetScopedName("last_level"));
}

BOOST_AUTO_TEST_CASE(ScopeStore_test) {
  auto scopeStack = std::make_shared<CScopeStack>();

  BOOST_CHECK_NO_THROW(CScopeStore store = CScopeStore(scopeStack));

  CScopeStore store = CScopeStore(scopeStack);
  BOOST_CHECK(store.GetCurrentScopeName() == "");

  scopeStack->push_back("scope_1");

  // Test copy constructor
  CScopeStore store_2(store);
  testConstructor(store, store_2);

  // Test assignment copy constructor
  CScopeStore store_3 = store;
  testConstructor(store, store_3);

  // Test move constructor
  CScopeStore store_4(std::move(store));
  testConstructor(store, store_4);

  // Test assignment move constructor
  CScopeStore store_5 = std::move(store);
  testConstructor(store, store_5);

  // Methods
  scopeStack->push_back("scope_2", ScopeType::STAGE);
  BOOST_CHECK(store.GetCurrentScopeName() == "scope_1.scope_2");
  BOOST_CHECK(store.GetScopedName("last_level") ==
              "scope_1.scope_2.last_level");
  BOOST_CHECK(store.getCurrentScopeNameAt(0) == "");
  BOOST_CHECK(store.getCurrentScopeNameAt(1) == "scope_1");
  BOOST_CHECK(store.getCurrentScopeNameAt(2) == "scope_1.scope_2");
  BOOST_CHECK_THROW(store.getCurrentScopeNameAt(3), AmzException);
  BOOST_CHECK(store.GetScopedNameAt("new_last_level", 2) ==
              "scope_1.scope_2.new_last_level");
  BOOST_CHECK(store.GetScopedNameAt("new_last_level", ScopeType::STAGE) ==
              "scope_1.scope_2.new_last_level");
  BOOST_CHECK(store.getScopeDepth() == 2);
}

BOOST_AUTO_TEST_SUITE_END()