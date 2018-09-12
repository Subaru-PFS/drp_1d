#include <RedshiftLibrary/reliability/pdfzPredictResult.h>
#include <RedshiftLibrary/reliability/pdfzFeatureResult.h>
#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/processflow/resultstore.h>
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

BOOST_AUTO_TEST_SUITE(ReliabilityzPDF)

BOOST_AUTO_TEST_CASE(PdfzPredictResult_Save)
{
  CPdfzPredictResult result;
  COperatorResultStore result_store;
  CParameterStore parameter_store;
  CDataStore dummy_store(result_store, parameter_store);
  stringstream stream(ios::in|ios::out);

  result.m_predLabel = "foo";
  result.m_Labels = TStringList({"bar", "baz", "blah" });

  result.m_posterior = gsl_vector_alloc(3);
  for (size_t i=0; i<3; ++i) {
    gsl_vector_set(result.m_posterior, i, i * 1.1);
  }

  result.m_score = gsl_vector_alloc(5);
  for (size_t i=0; i<5; ++i) {
    gsl_vector_set(result.m_score, i, i * 0.1);
  }

  result.Save(dummy_store, stream);

  stringstream expected("#Predicted_Label(zReliability)  \tfoo\n\
#Class \t PosteriorProbability_Prediction\t\n\
bar\t0\t\n\
baz\t1.1\t\n\
blah\t2.2\t\n\
\n#Learner \t Score_Prediction \t\n\
Learner_1\t0\t\n\
Learner_2\t0.1\t\n\
Learner_3\t0.2\t\n\
Learner_4\t0.3\t\n\
Learner_5\t0.4\t\n");
   BOOST_CHECK(stream.str() == expected.str());

  // BOOST_TEST_MESSAGE("\nresult:\n" << stream.str());
  // BOOST_TEST_MESSAGE("\nexpected:\n" << expected.str());
}


BOOST_AUTO_TEST_CASE(PdfzPredictResult_SaveLine)
{
  CPdfzPredictResult result;
  COperatorResultStore result_store;
  CParameterStore parameter_store;
  CDataStore data_store(result_store, parameter_store);
  stringstream stream(ios::in|ios::out);

  result.m_predLabel = "foo";
  result.m_predProba = 12.34;
  data_store.SetSpectrumName("blahblah");

  result.SaveLine(data_store, stream);

  stringstream expected("#zReliability_Pred \tfoo\tPosteriorProba_Pred \t12.34\tblahblah\n");
  BOOST_CHECK(stream.str() == expected.str());

  // BOOST_TEST_MESSAGE("\nresult:\n" << stream.str());
}

BOOST_AUTO_TEST_CASE(PdfzFeatureResult_Save)
{
  CPdfzFeatureResult result;
  COperatorResultStore result_store;
  CParameterStore parameter_store;
  CDataStore dummy_store(result_store, parameter_store);
  stringstream stream(ios::in|ios::out);

  result.mapzfeatures["foo"] = 11.11;
  result.mapzfeatures["bar"] = 12.21;
  result.mapzfeatures["baz"] = 13.31;

  result.Save(dummy_store, stream);

  stringstream expected("#zPDF_descriptors \tValue\n\
baz\t13.31\n\
bar\t12.21\n\
foo\t11.11\n");

  BOOST_CHECK(stream.str() == expected.str());
  //BOOST_TEST_MESSAGE("\nresult:\n" << stream.str());
  //BOOST_TEST_MESSAGE("\nexpected:\n" << expected.str());

  result.SaveLine(dummy_store, stream);
}

BOOST_AUTO_TEST_SUITE_END()

