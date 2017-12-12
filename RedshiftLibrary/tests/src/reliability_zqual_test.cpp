#include <RedshiftLibrary/reliability/zqualresult.h>
#include <RedshiftLibrary/reliability/zqual.h>
#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/processflow/resultstore.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>
#include <RedshiftLibrary/common/datatypes.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <boost/unordered_map.hpp>

#include <cstring>
#include <iostream>
#include <sstream>

#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(ReliabilityzQual)

BOOST_AUTO_TEST_CASE(QualzResult_Save)
{
  CQualzResult result;
  COperatorResultStore result_store;
  CParameterStore parameter_store;
  CDataStore data_store(result_store, parameter_store);
  std::stringstream stream;
  CDataStore::CAutoScope result_scope(data_store, "zReliability/result");
  //std::shared_ptr<const CQualzResult> SolveResult = (std::shared_ptr<const CQualzResult>) new CQualzResult();

  // data_store.StoreScopedGlobalResult("zpredict", (std::shared_ptr<const CQualzResult>) &result);

  // result_store.m_predLabel = "foo";
  // result_store.m_Labels = TStringList({"bar", "baz", "blah" });

  // result.m_posterior = gsl_vector_alloc(3);
  // for (size_t i=0; i<3; ++i) {
  //   gsl_vector_set(result.m_posterior, i, i * 1.1);
  // }

  // result.m_score = gsl_vector_alloc(5);
  // for (size_t i=0; i<5; ++i) {
  //   gsl_vector_set(result.m_score, i, i * 0.1);
  // }

  //result.Save(data_store, stream);

  BOOST_TEST_MESSAGE("FIXME: don't know how to test this");
}

BOOST_AUTO_TEST_SUITE_END()
