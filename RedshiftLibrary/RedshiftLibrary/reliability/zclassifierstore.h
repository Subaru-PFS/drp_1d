#ifndef _REDSHIFT_RELIABILITY_ZCLASSIFIERSTORE_
#define _REDSHIFT_RELIABILITY_ZCLASSIFIERSTORE_

#include <epic/core/common/datatypes.h>

#include <vector>
#include <gsl/gsl_matrix.h>
#include <boost/unordered_map.hpp>

namespace NSEpic
{


class CClassifierStore
{

private:

	class CLearner
	{
	public :

		CLearner();
		~CLearner();

		Bool Load ( gsl_matrix* params, gsl_matrix* vectors );

		Int32 m_nbSVectors;		// = M, number of SVector
		Int32 m_nbDesriptors;		// =  P, number of desriptors in an observation x
		Float64 m_SVbias;		// float

		gsl_vector* m_SVmu;		// size = [ 1 x P ]
		gsl_vector* m_SVsigma;		// size = [ 1 x P ]
		gsl_vector* m_SVbeta;		// size = [ P x 1 ]

		gsl_vector* m_SValpha;		// size = [ M x 1 ]
		gsl_vector* m_SVectorLabels;		// size =  [M x 1 ]
		gsl_matrix* m_SVectors;		// size = [ M x P ]

		TFloat64List m_SVsigmoiid ;		// [slope A, intercept C]  . SIGMOIID (x) = 1 / ( 1 + exp[ A*x + C ] )

	};

public:

	CClassifierStore();
	~CClassifierStore();

	Bool Load ( const char* dirPath );
	Bool Load_params( const char* dirPath );
    Bool m_isInitialized;

	Void DisplayQ ( const gsl_matrix* m );

	typedef boost::unordered_map<const Int32, std::shared_ptr<CLearner>> MapLearners;
	const MapLearners& GetLearners() const;
	Void SetLearners (MapLearners& learner );

	const Int32 GetNbClasses() const;
	const Int32 GetNbLearners() const;
	const Int32 GetNbFeatures() const;
	const std::string GetTypeCoding() const;
	const std::string GetTypeClassifier() const;
	const gsl_vector* GetLearnersWeight() const;
	const std::string GetLabel( Int32 idLabel ) const;
	const TStringList GetLabels() const;

	const gsl_matrix* GetCodingMatrix() const;
	const gsl_matrix* GetCodingMatrixPos() const;
	const gsl_matrix* GetCodingMatrixNeg() const;

	Void SetNbClasses( Int32 nbclasses );
	Void SetNbLearners( Int32 nblearners );
	Void SetNbFeatures( Int32 nbfeatures );
	Void SetTypeCoding( std::string typecoding );
	Void SetTypeClassifier( std::string typeclassifier );

	Void SetLearnerWeight( gsl_vector* w );
	Void SetCodingMatrix( gsl_matrix* m );
	Void SetCodingMatrixPos();
	Void SetCodingMatrixNeg();

	TFloat64List temp_sv;
	gsl_matrix* params_L;

protected:
	//typedef boost::unordered_map<const Int32 , int > MapLearners;
	MapLearners m_learners;

	CLearner* m_learner;

	TStringList m_Labels;

	Int32 m_nbClasses;
	Int32 m_nbLearners;
	Int32 m_nbFeatures;
	std::string m_typeCoding;
	std::string m_typeClassifier;

	gsl_vector* m_learnersWeight;

	gsl_matrix* m_codingMatrix;
	gsl_matrix* m_codingMatrixPos;
	gsl_matrix* m_codingMatrixNeg;

};


}
#endif

