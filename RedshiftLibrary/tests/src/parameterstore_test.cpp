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

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/range/combine.hpp>
#include <boost/test/unit_test.hpp>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/processflow/parameterstore.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(ParameterStore)

// Set parameters
TFloat64List float64_list{1.0, 2.0, 3.14};
TInt64List int64_list{42, 99, -100};
TBoolList bool_list{true, true, false};
TStringList string_list{"foo", "bar", "baz"};
TFloat64Range float64range(0.01, 0.011);
std::string string_val = "string";
Float64 float_val = 12.3;
Int64 int_val = 12;
bool bool_val = true;

std::string jsonString =
    "{\"TFloat64List\" : [ 1, 2, 3.1400000000000001 ],\"TInt64List\" : [ 42, "
    "99, -100 ],\"TBoolList\" : [ true, true, false ],\"TStringList\" : [ "
    "\"foo\", \"bar\", \"baz\" ],\"TFloat64Range\" : [ 0.01, "
    "0.010999999999999999 ],\"string\" : \"string\",\"Float64\" : "
    "12.300000000000001,\"Int64\" : 12,\"bool\" : true}";

// set json file
boost::filesystem::path _path =
    boost::filesystem::unique_path("file_%%%%%%%%%%");

BOOST_AUTO_TEST_CASE(ParameterStore_GetTest) {
  std::shared_ptr<CScopeStack> scopeStack = std::make_shared<CScopeStack>();

  CParameterStore store = CParameterStore(scopeStack);

  store.Set("TFloat64List", float64_list);
  store.Set("TInt64List", int64_list);
  store.Set("TBoolList", bool_list);
  store.Set("TStringList", string_list);
  store.Set("TFloat64Range", float64range);
  store.Set("string", string_val);
  store.Set("Float64", float_val);
  store.Set("Int64", int_val);
  store.Set("bool", bool_val);

  BOOST_CHECK_NO_THROW(store.Save(_path.c_str()));
  // BOOST_CHECK_NO_THROW(store_.Load(_path.c_str()));

  TFloat64List float64_list_ = store.GetList<Float64>("TFloat64List");
  TInt64List int64_list_ = store.GetList<Int64>("TInt64List");
  TBoolList bool_list_ = store.GetList<bool>("TBoolList");
  TStringList string_list_ = store.GetList<std::string>("TStringList");
  TFloat64Range float64range_ = store.Get<TFloat64Range>("TFloat64Range");
  std::string string_ = store.Get<std::string>("string");
  Float64 float64_ = store.Get<Float64>("Float64");
  Int64 int64_ = store.Get<Int64>("Int64");
  bool bool_ = store.Get<bool>("bool");

  Float64 f, f_;
  BOOST_FOREACH (boost::tie(f, f_),
                 boost::combine(float64_list, float64_list_)) {
    BOOST_TEST_MESSAGE("comparing " << f << " with " << f_);
    BOOST_CHECK_CLOSE(f, f_, 1e-6);
  }

  Int64 i, i_;
  BOOST_FOREACH (boost::tie(i, i_), boost::combine(int64_list, int64_list_)) {
    BOOST_TEST_MESSAGE("comparing " << i << " with " << i_);
    BOOST_CHECK(i == i_);
  }

  bool b, b_;
  BOOST_FOREACH (boost::tie(b, b_), boost::combine(bool_list, bool_list_)) {
    BOOST_TEST_MESSAGE("comparing " << b << " with " << b_);
    BOOST_CHECK(b == b_);
  }

  std::string s, s_;
  BOOST_FOREACH (boost::tie(s, s_), boost::combine(string_list, string_list_)) {
    BOOST_TEST_MESSAGE("comparing " << s << " with " << s_);
    BOOST_CHECK(s == s_);
  }

  BOOST_CHECK_CLOSE(float64range.GetBegin(), float64range_.GetBegin(), 1e-6);
  BOOST_CHECK_CLOSE(float64range.GetEnd(), float64range_.GetEnd(), 1e-6);

  BOOST_CHECK_CLOSE(float64_, 12.3, 1e-6);
  BOOST_CHECK(int64_ == 12);
  BOOST_CHECK(bool_ == true);
  BOOST_CHECK(string_ == "string");

  // test getScoped
  float64_ = store.GetScoped<Float64>("Float64");
  BOOST_CHECK_CLOSE(float64_, 12.3, 1e-6);

  // test getScopedAt
  float64_ = store.GetScopedAt<Float64>("Float64", 0);
  BOOST_CHECK_CLOSE(float64_, 12.3, 1e-6);

  scopeStack->push_back("model1", ScopeType::SPECTRUMMODEL);
  store.Set(store.GetScopedName("Float64"), float_val);
  float64_ = store.GetScopedAt<Float64>("Float64", ScopeType::SPECTRUMMODEL);
  BOOST_CHECK_CLOSE(float64_, 12.3, 1e-6);
  scopeStack->pop_back();

  // test getList
  float64_list_ = store.GetList<Float64>("TFloat64List");
  BOOST_FOREACH (boost::tie(f, f_),
                 boost::combine(float64_list, float64_list_)) {
    BOOST_TEST_MESSAGE("comparing " << f << " with " << f_);
    BOOST_CHECK_CLOSE(f, f_, 1e-6);
  }

  // test getListScoped
  float64_list_ = store.GetListScoped<Float64>("TFloat64List");
  BOOST_FOREACH (boost::tie(f, f_),
                 boost::combine(float64_list, float64_list_)) {
    BOOST_TEST_MESSAGE("comparing " << f << " with " << f_);
    BOOST_CHECK_CLOSE(f, f_, 1e-6);
  }

  // test throw
  BOOST_CHECK_THROW(store.Get<Float64>("Float64_"), AmzException);
  BOOST_CHECK_THROW(store.GetScoped<Float64>("Float64_"), AmzException);
  BOOST_CHECK_THROW(store.GetList<Float64>("TFloat64List_"), AmzException);
  BOOST_CHECK_THROW(store.GetListScoped<Float64>("TFloat64List_"),
                    AmzException);
}

BOOST_AUTO_TEST_CASE(ParameterStore_HasTest) {
  std::shared_ptr<CScopeStack> scopeStack = std::make_shared<CScopeStack>();
  CParameterStore store(scopeStack);
  store.Set("bool", bool_val);

  // test property
  bool check_property = store.Has<bool>("bool");
  BOOST_CHECK(check_property == true);

  check_property = store.Has<bool>("bool_");
  BOOST_CHECK(check_property == false);

  // test scope
  bool check_scope = store.HasScoped<bool>("bool");
  BOOST_CHECK(check_scope == true);

  check_scope = store.HasScoped<bool>("bool_");
  BOOST_CHECK(check_scope == false);
}

BOOST_AUTO_TEST_CASE(ParameterStore_ReadJsonTest) {
  // Test reading json string
  std::shared_ptr<CScopeStack> scopeStack = std::make_shared<CScopeStack>();
  CParameterStore store(scopeStack);

  store.FromString(jsonString);

  TFloat64List float64_list_ = store.GetList<Float64>("TFloat64List");
  TInt64List int64_list_ = store.GetList<Int64>("TInt64List");
  TBoolList bool_list_ = store.GetList<bool>("TBoolList");
  TStringList string_list_ = store.GetList<std::string>("TStringList");
  TFloat64Range float64range_ = store.Get<TFloat64Range>("TFloat64Range");
  std::string string_ = store.Get<std::string>("string");
  Float64 float64_ = store.Get<Float64>("Float64");
  Int64 int64_ = store.Get<Int64>("Int64");
  bool bool_ = store.Get<bool>("bool");

  Float64 f, f_;
  BOOST_FOREACH (boost::tie(f, f_),
                 boost::combine(float64_list, float64_list_)) {
    BOOST_TEST_MESSAGE("comparing " << f << " with " << f_);
    BOOST_CHECK_CLOSE(f, f_, 1e-6);
  }

  Int64 i, i_;
  BOOST_FOREACH (boost::tie(i, i_), boost::combine(int64_list, int64_list_)) {
    BOOST_TEST_MESSAGE("comparing " << i << " with " << i_);
    BOOST_CHECK(i == i_);
  }

  bool b, b_;
  BOOST_FOREACH (boost::tie(b, b_), boost::combine(bool_list, bool_list_)) {
    BOOST_TEST_MESSAGE("comparing " << b << " with " << b_);
    BOOST_CHECK(b == b_);
  }

  std::string s, s_;
  BOOST_FOREACH (boost::tie(s, s_), boost::combine(string_list, string_list_)) {
    BOOST_TEST_MESSAGE("comparing " << s << " with " << s_);
    BOOST_CHECK(s == s_);
  }

  BOOST_CHECK_CLOSE(float64range.GetBegin(), float64range_.GetBegin(), 1e-6);
  BOOST_CHECK_CLOSE(float64range.GetEnd(), float64range_.GetEnd(), 1e-6);

  BOOST_CHECK_CLOSE(float64_, 12.3, 1e-6);
  BOOST_CHECK(int64_ == 12);
  BOOST_CHECK(bool_ == true);
  BOOST_CHECK(string_ == "string");
}

BOOST_AUTO_TEST_CASE(ParameterStore_SpecificTest) {
  std::shared_ptr<CScopeStack> scopeStack = std::make_shared<CScopeStack>();
  CParameterStore store(scopeStack);
  std::string object_type = "galaxy";
  std::string object_type2 = "qso";
  bool check_property;

  // HasTplIsmExtinction
  store.Set(object_type + ".redshiftSolver.method",
            std::string("templateFittingSolve"));
  store.Set(object_type + ".redshiftSolver.templateFittingSolve.ismFit", true);
  check_property = store.HasTplIsmExtinction(object_type);
  BOOST_CHECK(check_property == true);

  store.Set(object_type + ".redshiftSolver.templateFittingSolve.ismFit", false);
  check_property = store.HasTplIsmExtinction(object_type);
  BOOST_CHECK(check_property == false);

  store.Set(object_type + ".redshiftSolver.method",
            std::string("tplCombinationSolve"));
  store.Set(object_type + ".redshiftSolver.tplCombinationSolve.ismFit", true);
  check_property = store.HasTplIsmExtinction(object_type);
  BOOST_CHECK(check_property == true);

  store.Set(object_type + ".redshiftSolver.tplCombinationSolve.ismFit", false);
  check_property = store.HasTplIsmExtinction(object_type);
  BOOST_CHECK(check_property == false);

  store.Set(object_type + ".redshiftSolver.method",
            std::string("lineModelSolve"));
  store.Set(object_type +
                ".redshiftSolver.lineModelSolve.lineModel.continuumFit.ismFit",
            true);
  check_property = store.HasTplIsmExtinction(object_type);
  BOOST_CHECK(check_property == true);

  store.Set(object_type +
                ".redshiftSolver.lineModelSolve.lineModel.continuumFit.ismFit",
            false);
  check_property = store.HasTplIsmExtinction(object_type);
  BOOST_CHECK(check_property == false);

  // HasTplIgmExtinction
  store.Set(object_type + ".redshiftSolver.method",
            std::string("templateFittingSolve"));
  store.Set(object_type + ".redshiftSolver.templateFittingSolve.igmFit", true);
  check_property = store.HasTplIgmExtinction(object_type);
  BOOST_CHECK(check_property == true);

  store.Set(object_type + ".redshiftSolver.templateFittingSolve.igmFit", false);
  check_property = store.HasTplIgmExtinction(object_type);
  BOOST_CHECK(check_property == false);

  store.Set(object_type + ".redshiftSolver.method",
            std::string("tplCombinationSolve"));
  store.Set(object_type + ".redshiftSolver.tplCombinationSolve.igmFit", true);
  // Here
  check_property = store.HasTplIgmExtinction(object_type);
  BOOST_CHECK(check_property == true);

  store.Set(object_type + ".redshiftSolver.tplCombinationSolve.igmFit", false);
  check_property = store.HasTplIgmExtinction(object_type);
  BOOST_CHECK(check_property == false);

  store.Set(object_type + ".redshiftSolver.method",
            std::string("lineModelSolve"));
  store.Set(object_type +
                ".redshiftSolver.lineModelSolve.lineModel.continuumFit.igmFit",
            true);
  check_property = store.HasTplIgmExtinction(object_type);
  BOOST_CHECK(check_property == true);

  store.Set(object_type +
                ".redshiftSolver.lineModelSolve.lineModel.continuumFit.igmFit",
            false);
  check_property = store.HasTplIgmExtinction(object_type);
  BOOST_CHECK(check_property == false);

  // HasFFTProcessing
  store.Set(object_type + ".redshiftSolver.method",
            std::string("templateFittingSolve"));
  store.Set(object_type + ".redshiftSolver.templateFittingSolve.fftProcessing",
            true);
  check_property = store.HasFFTProcessing(object_type);
  BOOST_CHECK(check_property == true);

  store.Set(object_type + ".redshiftSolver.templateFittingSolve.fftProcessing",
            false);
  store.Set(object_type + ".redshiftSolver.method",
            std::string("lineModelSolve"));
  store.Set(
      object_type +
          ".redshiftSolver.lineModelSolve.lineModel.continuumFit.fftProcessing",
      true);
  check_property = store.HasFFTProcessing(object_type);
  BOOST_CHECK(check_property == true);

  store.Set(
      object_type +
          ".redshiftSolver.lineModelSolve.lineModel.continuumFit.fftProcessing",
      false);
  check_property = store.HasFFTProcessing(object_type);
  BOOST_CHECK(check_property == false);

  // hasToLogRebin
  store.Set(
      object_type +
          ".redshiftSolver.lineModelSolve.lineModel.continuumFit.fftProcessing",
      true);
  store.Set(
      object_type2 +
          ".redshiftSolver.lineModelSolve.lineModel.continuumFit.fftProcessing",
      true);
  std::map<std::string, bool> fft_processing;
  bool m_use_LogLambaSpectrum =
      store.hasToLogRebin({object_type, object_type2}, fft_processing);
  BOOST_CHECK(fft_processing[object_type] = true);
  BOOST_CHECK(fft_processing[object_type2] = true);
  BOOST_CHECK(m_use_LogLambaSpectrum = true);

  // getMinZStepForFFTProcessing
  store.Set(object_type + ".redshiftStep", 0.0005);
  store.Set(object_type2 + ".redshiftStep", 0.0003);
  Float64 logGridStep = store.getMinZStepForFFTProcessing(fft_processing);
  BOOST_CHECK(logGridStep == 0.0003);

  // HasToOrthogonalizeTemplates
  store.Set(object_type + ".redshiftSolver.method",
            std::string("lineModelSolve"));
  store.Set(object_type +
                ".redshiftSolver.lineModelSolve.lineModel.continuumComponent",
            std::string("tplFit"));
  check_property = store.HasToOrthogonalizeTemplates(object_type);
  BOOST_CHECK(check_property == true);

  store.Set(object_type +
                ".redshiftSolver.lineModelSolve.lineModel.continuumComponent",
            std::string("tplFitAuto"));
  check_property = store.HasToOrthogonalizeTemplates(object_type);
  BOOST_CHECK(check_property == true);

  store.Set(object_type +
                ".redshiftSolver.lineModelSolve.lineModel.continuumComponent",
            std::string("tplFitAuto__"));
  check_property = store.HasToOrthogonalizeTemplates(object_type);
  BOOST_CHECK(check_property == false);

  // EnableTemplateOrthogonalization
  check_property = store.EnableTemplateOrthogonalization(object_type);
  BOOST_CHECK(check_property == false);

  store.Set(object_type +
                ".redshiftSolver.lineModelSolve.lineModel.continuumComponent",
            std::string("tplFitAuto"));
  store.Set(object_type + ".redshiftSolver.lineModelSolve.lineModel."
                          "continuumFit.ignoreLineSupport",
            false);
  check_property = store.EnableTemplateOrthogonalization(object_type);
  BOOST_CHECK(check_property == true);

  store.Set(object_type + ".redshiftSolver.lineModelSolve.lineModel."
                          "continuumFit.ignoreLineSupport",
            true);
  check_property = store.EnableTemplateOrthogonalization(object_type);
  BOOST_CHECK(check_property == false);
}

BOOST_AUTO_TEST_CASE(ParameterStore_saveTest) {
  std::shared_ptr<CScopeStack> scopeStack = std::make_shared<CScopeStack>();
  CParameterStore store(scopeStack);

  try {
    store.Save("/this/file/should/not/exist");
    BOOST_FAIL("store.Save() should have failed");
  } catch (std::runtime_error &) {
    BOOST_CHECK(true);
  }
}

BOOST_AUTO_TEST_SUITE_END()
