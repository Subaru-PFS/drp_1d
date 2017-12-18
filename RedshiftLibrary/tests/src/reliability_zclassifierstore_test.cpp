#include <RedshiftLibrary/reliability/zqualresult.h>
#include <RedshiftLibrary/reliability/zclassifierstore.h>
#include <RedshiftLibrary/reliability/zqual.h>
#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/processflow/resultstore.h>
#include <RedshiftLibrary/common/datatypes.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <boost/unordered_map.hpp>
#include <boost/filesystem.hpp>

#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>

#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(ZClassifierStore)

const Int32 ROWS = 8+1, FEATURES=2, CLASSES=5, LEARNERS=5;

void create_test_file_params(boost::filesystem::path path)
{
  FILE *f;
  gsl_vector* params;
  boost::filesystem::path filename;

  filename = path/"zClassifier_params_ova.dat";
  f = fopen(filename.c_str(), "w");
  params = gsl_vector_alloc(ROWS);
  params->data[0] = FEATURES;
  params->data[1] = CLASSES;
  params->data[2] = LEARNERS;
  for (Int32 i=3; i<ROWS; ++i) {
    gsl_vector_set(params, i,  (i-2)*2);
    // so that temp_sv[i] = (i+1)*2. ie : sv_{i}_vector.dat has (i+1)*2 rows
  }

  gsl_vector_fprintf(f, params,"%f");
  gsl_vector_free( params );
  fclose(f);
}

void create_test_file_weigths(boost::filesystem::path path)
{
  FILE *f;
  gsl_vector* params;
  boost::filesystem::path filename;

  filename = path/"zClassifier_learnersWeight.dat";
  f = fopen(filename.c_str(), "w");
  params = gsl_vector_alloc(LEARNERS + 1);
  for (Int32 i=0; i<LEARNERS+1; ++i) {
    gsl_vector_set(params, i, (Float64) i*0.123);
  }

  gsl_vector_fprintf(f, params, "%f");
  gsl_vector_free( params );
  fclose(f);
}

void create_test_file_codingmatrix(boost::filesystem::path path)
{
  FILE *f;
  gsl_matrix* coding_matrix;
  boost::filesystem::path filename;

  filename = path/"zClassifier_codingMatrix.dat";
  f = fopen(filename.c_str(), "w");
  coding_matrix = gsl_matrix_alloc(CLASSES, LEARNERS);

  gsl_matrix_set_zero(coding_matrix);
  for (Int32 i=0; i<std::min(CLASSES, LEARNERS); ++i) {
    for (Int32 j=0; j<=i; ++j) {
      gsl_matrix_set(coding_matrix, j, i, std::pow(-1, i+j));
    }
  }

  gsl_matrix_fprintf(f, coding_matrix, "%f");
  gsl_matrix_free( coding_matrix );
  fclose(f);
}

void create_test_file_learnersparams(boost::filesystem::path path)
{
  FILE *f;
  gsl_matrix* learners_matrix;
  boost::filesystem::path filename;

  filename = path/"zClassifier_learnersParams.dat";
  f = fopen(filename.c_str(), "w");
  learners_matrix = gsl_matrix_alloc(3, LEARNERS);
  for (Int32 i=0; i<3; ++i) {
    for (Int32 j=0; j<LEARNERS; ++j) {
      gsl_matrix_set(learners_matrix, i, j, (Float64) i*0.1 + j*2.2);
    }
  }

  gsl_matrix_fprintf(f, learners_matrix, "%f");
  gsl_matrix_free( learners_matrix );
  fclose(f);
}

void create_test_file_learners_vectors(boost::filesystem::path path)
{
  FILE *f;
  gsl_matrix* sv_vectors;
  boost::filesystem::path filename;

  for (Int32 i=0; i<LEARNERS; ++i) {
    boost::filesystem::path pathname = path/("sv" + std::to_string(i+1));
    boost::filesystem::create_directory(pathname);
    filename = pathname/("sv" + std::to_string(i+1) + "_vectors.dat");
    f = fopen(filename.c_str(), "w");
    sv_vectors = gsl_matrix_alloc( (i+1)*2, FEATURES + 2);
    // (i+1)*2 rows, as created in create_test_file_params
    gsl_matrix_set_identity(sv_vectors);
    gsl_matrix_fprintf(f, sv_vectors, "%f");
    gsl_matrix_free( sv_vectors );
    fclose(f);
  }
}

void create_test_file_learners_params(boost::filesystem::path path)
{
  FILE *f;
  gsl_matrix* sv_params;
  boost::filesystem::path filename;

  for (Int32 i=0; i<LEARNERS; ++i) {
    boost::filesystem::path pathname = path/("sv" + std::to_string(i+1));
    boost::filesystem::create_directory(pathname);
    filename = pathname/("sv" + std::to_string(i+1) + "_params.dat");
    f = fopen(filename.c_str(), "w");
    sv_params = gsl_matrix_alloc(FEATURES, 3);

    gsl_matrix_set_all(sv_params, 2.0);
    gsl_matrix_fprintf(f, sv_params, "%f");
    gsl_matrix_free( sv_params );
    fclose(f);
  }
}

void delete_test_files()
{
  boost::filesystem::remove("/tmp/zClassifier_params_ova.dat");
  boost::filesystem::remove("/tmp/zClassifier_learnersWeight.dat");
  boost::filesystem::remove("/tmp/zClassifier_codingMatrix.dat");
  boost::filesystem::remove("/tmp/zClassifier_learnersParams.dat");
  for (Int32 i=0; i<LEARNERS; ++i) {
    boost::filesystem::path pathname = "/tmp/sv" + std::to_string(i+1);
    boost::filesystem::remove_all(pathname);
  }
}

BOOST_AUTO_TEST_CASE(ZClassifierStore1)
{
  CClassifierStore store;

  delete_test_files();

  BOOST_TEST_MESSAGE("The following error regarding zClassifier_params_ova.dat is normal:");
  BOOST_CHECK(store.Load_params("/this/path/should/not/exist") == false);
  create_test_file_params("/tmp");
  BOOST_TEST_MESSAGE("The following error regarding zClassifier_learnersWeight.dat is normal:");
  BOOST_CHECK(store.Load_params("/tmp") == false);
  create_test_file_weigths("/tmp");
  BOOST_CHECK(store.Load_params("/tmp") == true);

  BOOST_CHECK( store.GetNbClasses() == CLASSES );
  BOOST_CHECK( store.GetNbLearners() == LEARNERS );
  BOOST_CHECK( store.GetNbFeatures() == FEATURES );

  delete_test_files();
}

BOOST_AUTO_TEST_CASE(ZClassifierStore2)
{
  CClassifierStore store;

  delete_test_files();

  create_test_file_params("/tmp");
  create_test_file_weigths("/tmp");
  BOOST_TEST_MESSAGE("The following error regarding zClassifier_codingMatrix.dat is normal:");
  BOOST_CHECK(store.Load("/tmp") == false);
  create_test_file_codingmatrix("/tmp");
  BOOST_TEST_MESSAGE("The following error regarding Classifier_learnersParams.dat is normal:");
  BOOST_CHECK(store.Load("/tmp") == false);
  create_test_file_learnersparams("/tmp");
  BOOST_TEST_MESSAGE("The following error regarding sv1/sv1_vectors.dat is normal:");
  BOOST_CHECK(store.Load("/tmp") == false);
  create_test_file_learners_vectors("/tmp");
  BOOST_TEST_MESSAGE("The following error regarding sv1/sv1_params.dat is normal:");
  BOOST_CHECK(store.Load("/tmp") == false);
  create_test_file_learners_params("/tmp");
  BOOST_CHECK(store.Load("/tmp") == true);

  delete_test_files();

}

BOOST_AUTO_TEST_CASE(Qualz_Solve)
{
  CQualz zqual;
  COperatorResultStore result_store;
  CParameterStore parameter_store;
  CDataStore data_store(result_store, parameter_store);
  CDataStore::CAutoScope result_scope(data_store, "linemodelsolve");

  CClassifierStore classifier_store;
  TFloat64Range redshift_range(0.8,1.9);
  Float64 redshift_step=0.1;
  auto postmargZResult = std::shared_ptr<CPdfMargZLogResult>(new CPdfMargZLogResult());

  postmargZResult->Redshifts = TFloat64List({0.3, 0.4, 0.5, 0.6, 0.7});
  postmargZResult->valProbaLog = TFloat64List({0.1, 0.2, 0.3, 0.4, 0.5});
  postmargZResult->countTPL = postmargZResult->Redshifts.size();

  data_store.StoreScopedGlobalResult( "zPDF/logposterior.logMargP_Z_data", postmargZResult);

  delete_test_files();
  create_test_file_params("/tmp");
  create_test_file_weigths("/tmp");
  create_test_file_codingmatrix("/tmp");
  create_test_file_learnersparams("/tmp");
  create_test_file_learners_vectors("/tmp");
  create_test_file_learners_params("/tmp");
  BOOST_CHECK(classifier_store.Load("/tmp") == true);

  delete_test_files();

  zqual.Compute(data_store, classifier_store, redshift_range, redshift_step);

}


BOOST_AUTO_TEST_SUITE_END()
