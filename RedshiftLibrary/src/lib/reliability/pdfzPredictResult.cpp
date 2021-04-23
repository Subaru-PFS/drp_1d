#include <RedshiftLibrary/reliability/pdfzPredictResult.h>

#include <RedshiftLibrary/processflow/datastore.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>

#include <iostream>
#include <iomanip>

using namespace std;
using namespace NSEpic;


CPdfzPredictResult::CPdfzPredictResult()
{

}

CPdfzPredictResult::~CPdfzPredictResult()
{
    if(m_score != nullptr) gsl_vector_free(m_score);
    if(m_posterior != nullptr) gsl_vector_free(m_posterior);
}


