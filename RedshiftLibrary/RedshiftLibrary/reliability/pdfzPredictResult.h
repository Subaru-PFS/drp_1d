#ifndef _REDSHIFT_OPERATOR_PDFZPREDICTRESULT_
#define _REDSHIFT_OPERATOR_PDFZPREDICTRESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/operator/operator.h>
#include <RedshiftLibrary/common/datatypes.h>

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

	void Save( const CDataStore& store, std::ostream& stream ) const;
	void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }

	TStringList m_Labels;

	gsl_vector* m_score;
	gsl_vector* m_posterior;

	Float64 m_predProba;
	std::string m_predLabel;

};


}

#endif
