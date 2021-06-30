#ifndef _REDSHIFT_RELIABILITY_PDFZPREDICTRESULT_
#define _REDSHIFT_RELIABILITY_PDFZPREDICTRESULT_

#include "RedshiftLibrary/processflow/result.h"
#include "RedshiftLibrary/operator/operator.h"
#include "RedshiftLibrary/common/datatypes.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <boost/unordered_map.hpp>

using namespace std;
namespace NSEpic
{

class CPdfzPredictResult : public COperatorResult
{

public:

	CPdfzPredictResult();
	virtual ~CPdfzPredictResult();


	TStringList m_Labels;

	gsl_vector* m_score = nullptr;
	gsl_vector* m_posterior = nullptr;

	Float64 m_predProba;
	std::string m_predLabel;

};


}

#endif
