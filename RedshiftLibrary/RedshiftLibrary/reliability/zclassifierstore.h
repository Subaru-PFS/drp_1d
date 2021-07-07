#ifndef _REDSHIFT_RELIABILITY_ZCLASSIFIERSTORE_
#define _REDSHIFT_RELIABILITY_ZCLASSIFIERSTORE_

#include "RedshiftLibrary/common/datatypes.h"

#include <boost/unordered_map.hpp>
#include <gsl/gsl_matrix.h>
#include <vector>

namespace NSEpic {

class CClassifierStore
{

  private:
    class CLearner
    {
      public:
        CLearner();
        ~CLearner();

        Bool LoadParams(gsl_matrix *params);
        Bool LoadVectors(gsl_matrix *vectors);

        Int32 m_nbSVectors = 0;   // = M, number of SVector
        Int32 m_nbDesriptors = 0; // =  P, number of desriptors in an observation x
        Float64 m_SVbias;     // float

        gsl_vector *m_SVmu = NULL;    // size = [ 1 x P ]
        gsl_vector *m_SVsigma = NULL; // size = [ 1 x P ]
        gsl_vector *m_SVbeta= NULL;  // size = [ P x 1 ]

        gsl_vector *m_SValpha = NULL;       // size = [ M x 1 ]
        gsl_vector *m_SVectorLabels = NULL; // size =  [M x 1 ]
        gsl_matrix *m_SVectors = NULL;      // size = [ M x P ]

        TFloat64List m_SVsigmoiid; // [slope A, intercept C]  . SIGMOIID (x) = 1
                                   // / ( 1 + exp[ A*x + C ] )
    };

  public:
    CClassifierStore();
    ~CClassifierStore();

    Bool Load(const char *dirPath);
    Float64 Load_version(const char *dirPath);
    Bool Load_params(const char *dirPath);
    Bool m_isInitialized;
    Float64 m_file_format_version;

    void DisplayQ(const gsl_matrix *m);

    typedef boost::unordered_map<const Int32, std::shared_ptr<CLearner>>
        MapLearners;
    const MapLearners &GetLearners() const;
    void SetLearners(MapLearners &learner);

    const Int32 GetNbClasses() const;
    const Int32 GetNbLearners() const;
    const Int32 GetNbFeatures() const;
    const std::string GetTypeCoding() const;
    const std::string GetTypeClassifier() const;
    const Int32 GetOptionClassifier() const;

    const gsl_vector *GetLearnersWeight() const;
    const std::string GetLabel(Int32 idLabel) const;
    const TStringList GetLabels() const;

    const gsl_matrix *GetCodingMatrix() const;
    const gsl_matrix *GetCodingMatrixPos() const;
    const gsl_matrix *GetCodingMatrixNeg() const;

    void SetNbClasses(Int32 nbclasses);
    void SetNbLearners(Int32 nblearners);
    void SetNbFeatures(Int32 nbfeatures);
    void SetTypeCoding(std::string typecoding);
    void SetTypeClassifier(std::string typeclassifier);

    void SetLearnerWeight(gsl_vector *w);
    void SetCodingMatrix(gsl_matrix *m);
    void SetCodingMatrixPos();
    void SetCodingMatrixNeg();

    TFloat64List temp_sv;
    gsl_matrix *params_L = NULL;

  protected:
    // typedef boost::unordered_map<const Int32 , int > MapLearners;
    MapLearners m_learners;

    CLearner *m_learner;

    TStringList m_Labels;

    Int32 m_nbClasses = 0;
    Int32 m_nbLearners = 0;
    Int32 m_nbFeatures = 0;
    std::string m_typeCoding;
    std::string m_typeClassifier;
    Int32 m_classifier_option;

    gsl_vector *m_learnersWeight = NULL;

    gsl_matrix *m_codingMatrix = NULL;
    gsl_matrix *m_codingMatrixPos = NULL;
    gsl_matrix *m_codingMatrixNeg = NULL;
};

} // namespace NSEpic
#endif
