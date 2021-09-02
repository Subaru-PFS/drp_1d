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
#ifndef _REDSHIFT_RELIABILITY_ZQUAL_
#define _REDSHIFT_RELIABILITY_ZQUAL_

#define EPS_TOL 2.2204e-16

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/operator/operator.h"
#include "RedshiftLibrary/processflow/result.h"

#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/debug/assert.h"
#include "RedshiftLibrary/processflow/datastore.h"
#include "RedshiftLibrary/processflow/resultstore.h"

#include <vector>
#include <gsl/gsl_matrix.h>
#include <boost/unordered_map.hpp>

#include "RedshiftLibrary/reliability/pdfzFeatureResult.h"
#include "RedshiftLibrary/operator/pdfMargZLogResult.h"

#include "RedshiftLibrary/reliability/zqualresult.h"
#include "RedshiftLibrary/reliability/zclassifierstore.h"

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
    Bool disp_time    = false;    // display time of full process (zFeature + zProject)
	Bool m_doTEST   = 0;         // matlab check

	boost::posix_time::time_duration T0, T1;


private:

	Bool Solve( CDataStore& resultStore, CClassifierStore& classifierStore, const TFloat64Range& redshiftRange,
			Float64& redshiftStep );

    void DisplayPrediction();
	Bool CheckPDF ( const TFloat64List& zpdf);

	// extract descriptors from the zPDF
	Bool ExtractFeaturesPDF ( CDataStore& resultStore, const TFloat64Range& redshiftRange, Float64& redshiftStep );

	// predict a reliability label using const classifier
	void ProjectPDF ( CClassifierStore& classifierStore );
	void GetPosteriorPred ( CClassifierStore& classifierStore );
	void GetScorePred ( CClassifierStore& classifierStore );
	void GetLabelPred ( CClassifierStore& classifierStore );

	// method to derive number of significant peaks/modes in the zPDF
	Int32 GetNPeaksKM ( TPointList& peaks, Float64& zmap_estimate, Float64& lbins );
	Float64 GetDistanceKM ( gsl_vector *v1, gsl_vector *v2, std::string method );

	// KL minimization to compute posterior class probablities for each prediction
	gsl_vector* GetArgminKL ( CClassifierStore& classifierStore, gsl_vector* r, gsl_vector* w_learner, gsl_vector* p0, Float64& distance );
	gsl_vector* GetNumDenKL( CClassifierStore& classifierStore, gsl_vector* r, gsl_vector* r_estim , gsl_vector* pold );
	Float64 GetDistanceKL ( gsl_vector* r, gsl_vector* r_estim, gsl_vector* w_learner );
	void GetSigmoid ( Float64& sc, TFloat64List& paramsSig, Float64& result );
	Bool GetLSQnonNegKL (gsl_matrix* c, gsl_vector* d, gsl_vector* result );
	gsl_vector* GetProductKL ( gsl_vector* x, gsl_matrix* sv );
	gsl_matrix* GetTimesKL ( const gsl_matrix* m, gsl_vector* v );
	gsl_vector* GetSumKL ( gsl_matrix* m, Bool opt_row );
	void GetXc ( gsl_vector* mu, gsl_vector* sigma );
	//Float64 GetNormKL ( gsl_matrix* m );



protected:

	Bool nanVector = false;			// sometimes the flux and noise spectra are irrelevant, which produces a false PDF...

	TFloat64List m_zfeatures;
	vector<string> id_descriptors;

	Int32 m_idpredLabel = 0;
	Float64 m_predProba = 0;
	std::string m_predLabel = "";

	gsl_vector* m_Xc = NULL;
	gsl_vector* m_score = NULL;
	gsl_vector* m_posterior = NULL;

};

}
#endif
