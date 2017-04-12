#ifndef _REDSHIFT_RELIABILITY_ZQUAL_
#define _REDSHIFT_RELIABILITY_ZQUAL_

#define EPS_TOL 2.2204e-16

#include <epic/core/common/datatypes.h>
#include <epic/redshift/operator/operator.h>
#include <epic/redshift/processflow/result.h>

#include <epic/core/log/log.h>
#include <epic/core/debug/assert.h>
#include <epic/redshift/processflow/datastore.h>
#include <epic/redshift/processflow/resultstore.h>

#include <vector>
#include <gsl/gsl_matrix.h>
#include <boost/unordered_map.hpp>

#include <epic/redshift/reliability/pdfzFeatureResult.h>
#include <epic/redshift/operator/pdfMargZLogResult.h>

#include <epic/redshift/reliability/zqualresult.h>
#include <epic/redshift/reliability/zclassifierstore.h>

namespace NSEpic
{
class CDataStore;

class CQualz
{
public:

	CQualz();
	~CQualz();

	std::shared_ptr<const CQualzResult> Compute ( CDataStore& resultStore, CClassifierStore& classifierStore,
			const TFloat64Range &redshiftRange, Float64& redshiftStep );

	Bool disp_details = false;   // display temporary results for each learnerwhen computing the score
	Bool disp_time    = true;    // display time of full process (zFeature + zProject)
	Bool m_doTEST   = 0;    // matlab check

	boost::posix_time::time_duration T0, T1;


private:

	Bool Solve( CDataStore& resultStore, CClassifierStore& classifierStore, const TFloat64Range& redshiftRange,
			Float64& redshiftStep );

	Void DisplayPrediction();
	Bool CheckPDF ( const TFloat64List& zpdf);

	// extract descriptors from the zPDF
	Bool ExtractFeaturesPDF ( CDataStore& resultStore, const TFloat64Range& redshiftRange, Float64& redshiftStep );

	// predict a reliability label using const classifier
	Void ProjectPDF ( CClassifierStore& classifierStore );
	Void GetPosteriorPred ( CClassifierStore& classifierStore );
	Void GetScorePred ( CClassifierStore& classifierStore );
	Void GetLabelPred ( CClassifierStore& classifierStore );

	// method to derive number of significant peaks/modes in the zPDF
	Int32 GetNPeaksKM ( TPointList& peaks, Float64& zmap_estimate, Float64& lbins );
	Float64 GetDistanceKM ( gsl_vector *v1, gsl_vector *v2, std::string method );

	// KL minimization to compute posterior class probablities for each prediction
	gsl_vector* GetArgminKL ( CClassifierStore& classifierStore, gsl_vector* r, gsl_vector* w_learner, gsl_vector* p0, Float64& distance );
	gsl_vector* GetNumDenKL( CClassifierStore& classifierStore, gsl_vector* r, gsl_vector* r_estim , gsl_vector* pold );
	Float64 GetDistanceKL ( gsl_vector* r, gsl_vector* r_estim, gsl_vector* w_learner );
	Void GetSigmoid ( Float64& sc, TFloat64List& paramsSig, Float64& result );
	Bool GetLSQnonNegKL (gsl_matrix* c, gsl_vector* d, gsl_vector* result );
	gsl_vector* GetProductKL ( gsl_vector* x, gsl_matrix* sv );
	gsl_matrix* GetTimesKL ( const gsl_matrix* m, gsl_vector* v );
	gsl_vector* GetSumKL ( gsl_matrix* m, Bool opt_row );
	Void GetXc ( gsl_vector* mu, gsl_vector* sigma );
	//Float64 GetNormKL ( gsl_matrix* m );



protected:

	Bool nanVector = false;			// sometimes the flux and noise spectra are irrelevant, which produces a false PDF...

	TFloat64List m_zfeatures;
	vector<string> id_descriptors;

	Int32 m_idpredLabel;
	Float64 m_predProba;
	std::string m_predLabel;

	gsl_vector* m_Xc;
	gsl_vector* m_score;
	gsl_vector* m_posterior;
	gsl_vector* m_binaryPred; 	// prediction y_0 for each learner {-1; 1}

};

}
#endif
