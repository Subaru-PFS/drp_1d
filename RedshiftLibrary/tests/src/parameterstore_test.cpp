#include <RedshiftLibrary/processflow/parameterstore.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>

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
  CParameterStore store = CParameterStore();
  CParameterStore store_;

  TFloat64List float64_list {1.0, 2.0, 3.14};
  TInt64List int64_list {42, 99, -100};
  TBoolList bool_list {true, true, false};
  TStringList string_list {"foo", "bar", "baz"};
  TFloat64Range float64range(0.01,0.011);

  TFloat64List float64_list_;
  TInt64List int64_list_;
  TBoolList bool_list_;
  TStringList string_list_;
  TFloat64Range float64range_;
  std::string string_;
  Float64 float64_;
  Int64 int64_;
  Bool bool_;

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
  BOOST_CHECK_NO_THROW(store_.Load(_path.c_str()));

  store_.Get( "TFloat64List", float64_list_ );
  store_.Get( "TInt64List", int64_list_ );
  store_.Get( "TBoolList", bool_list_ );
  store_.Get( "TStringList", string_list_ );
  store_.Get( "TFloat64Range", float64range_ );
  store_.Get( "string", string_);
  store_.Get( "Float64", float64_);
  store_.Get( "Int64", int64_);
  store_.Get( "Bool", bool_);

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

BOOST_AUTO_TEST_CASE(ParameterStore2)
{
  CParameterStore store = CParameterStore();

  TFloat64List float64_list {1.0, 2.0, 3.14};
  TInt64List int64_list {42, 99, -100};
  TBoolList bool_list {true, true, false};
  TStringList string_list {"foo", "bar", "baz"};
  TFloat64Range float64range(0.01,0.011);

  TFloat64List float64_list_;
  TInt64List int64_list_;
  TBoolList bool_list_;
  TStringList string_list_;
  TFloat64Range float64range_;
  std::string string_;
  Float64 float64_;
  Int64 int64_;
  Bool bool_;

  store.Get( "TFloat64List", float64_list_, float64_list );
  store.Get( "TInt64List", int64_list_, int64_list );
  store.Get( "TBoolList", bool_list_, bool_list );
  store.Get( "TStringList", string_list_, string_list );
  store.Get( "TFloat64Range", float64range_, float64range );
  store.Get( "string", string_, std::string("foo"));
  store.Get( "Float64", float64_, Float64(12.3) );
  store.Get( "Int64", int64_, Int64(12));
  store.Get( "Bool", bool_, Bool(true));

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
  CParameterStore store;

  try {
    store.Load("/this/file/should/not/exist");
    BOOST_FAIL("store.Load() should have failed");
  } catch (std::runtime_error&) {
    BOOST_CHECK(true);
  }

  try {
    store.Save("/this/file/should/not/exist");
    BOOST_FAIL("store.Save() should have failed");
  } catch (std::runtime_error&) {
    BOOST_CHECK(true);
  }

}


BOOST_AUTO_TEST_SUITE_END()
