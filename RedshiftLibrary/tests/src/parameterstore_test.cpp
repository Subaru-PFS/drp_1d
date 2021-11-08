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
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"

#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <boost/range/combine.hpp>
#include <boost/foreach.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(ParameterStore)

BOOST_AUTO_TEST_CASE(ParameterStore1)
{
  TScopeStack scopeStack;

  CParameterStore store = CParameterStore(scopeStack);
  //CParameterStore store_(scopeStack);

  TFloat64List float64_list {1.0, 2.0, 3.14};
  TInt64List int64_list {42, 99, -100};
  TBoolList bool_list {true, true, false};
  TStringList string_list {"foo", "bar", "baz"};
  TFloat64Range float64range(0.01,0.011);

  store.Set( "TFloat64List", float64_list );
  store.Set( "TInt64List", int64_list );
  store.Set( "TBoolList", bool_list );
  store.Set( "TStringList", string_list );
  store.Set( "TFloat64Range", float64range );
  store.Set( "string", "string");
  store.Set( "Float64", 12.3);
  store.Set( "Int64", Int64(12));
  store.Set( "Bool", Bool(true));

  boost::filesystem::path _path = boost::filesystem::unique_path("file_%%%%%%%%%%");
  BOOST_CHECK_NO_THROW(store.Save(_path.c_str()));
  //BOOST_CHECK_NO_THROW(store_.Load(_path.c_str()));

  TFloat64List float64_list_ = store.GetList<Float64>("TFloat64List");
  TInt64List int64_list_ = store.GetList<Int64>("TInt64List");
  TBoolList bool_list_ = store.GetList<Bool>("TBoolList");
  TStringList string_list_ = store.GetList<std::string>("TStringList");
  TFloat64Range float64range_ = store.Get<TFloat64Range>("TFloat64Range");
  std::string string_ = store.Get<std::string>("string");
  Float64 float64_ = store.Get<Float64>("Float64");
  Int64 int64_ = store.Get<Int64>("Int64");
  Bool bool_ = store.Get<Bool>("Bool");

  Float64 f, f_;
  BOOST_FOREACH(boost::tie(f, f_), boost::combine(float64_list, float64_list_)) {
    BOOST_TEST_MESSAGE("comparing " << f << " with " << f_);
    BOOST_CHECK_CLOSE(f, f_, 1e-6 );
  }

  Int64 i, i_;
  BOOST_FOREACH(boost::tie(i, i_), boost::combine(int64_list, int64_list_)) {
    BOOST_TEST_MESSAGE("comparing " << i << " with " << i_);
    BOOST_CHECK( i == i_ );
  }

  Bool b, b_;
  BOOST_FOREACH(boost::tie(b, b_), boost::combine(bool_list, bool_list_)) {
    BOOST_TEST_MESSAGE("comparing " << b << " with " << b_);
    BOOST_CHECK( b == b_ );
  }

  std::string s, s_;
  BOOST_FOREACH(boost::tie(s, s_), boost::combine(string_list, string_list_)) {
    BOOST_TEST_MESSAGE("comparing " << s << " with " << s_);
    BOOST_CHECK( s == s_ );
  }

  BOOST_CHECK_CLOSE( float64range.GetBegin(), float64range_.GetBegin(), 1e-6);
  BOOST_CHECK_CLOSE( float64range.GetEnd(), float64range_.GetEnd(), 1e-6);

  BOOST_CHECK_CLOSE( float64_, 12.3, 1e-6 );
  BOOST_CHECK( int64_ == 12 );
  BOOST_CHECK( bool_ == true );

}

BOOST_AUTO_TEST_CASE(ParameterStore3)
{
  TScopeStack scopeStack;
  CParameterStore store(scopeStack);
/*
  try {
    store.Load("/this/file/should/not/exist");
    BOOST_FAIL("store.Load() should have failed");
  } catch (std::runtime_error&) {
    BOOST_CHECK(true);
  }*/

  try {
    store.Save("/this/file/should/not/exist");
    BOOST_FAIL("store.Save() should have failed");
  } catch (std::runtime_error&) {
    BOOST_CHECK(true);
  }

}


BOOST_AUTO_TEST_SUITE_END()
